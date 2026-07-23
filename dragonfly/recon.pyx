import os
import time
import numpy as np
import h5py
from mpi4py import MPI

cimport numpy as np
cimport openmp
from libc.stdlib cimport malloc, calloc, free
from libc.stdint cimport uint8_t
from libc.math cimport isnan, sqrt
from .iterate cimport Iterate
from . cimport recon as c_recon
from . cimport iterate as c_iterate
from . cimport model as c_model
from . cimport params as c_params
from . cimport quaternion as c_quat
from . cimport emcfile as c_dataset

cdef class EMCRecon():
    '''Main class for EMC reconstruction.

    Performs the Expand-Maximize-Compress algorithm for single-particle
    diffraction imaging using MPI and OpenMP.

    Args:
        num_threads (int): Number of OpenMP threads. Default -1 (auto).

    Example:
        >>> recon = EMCRecon(num_threads=8)
        >>> itr = Iterate('config.ini')
        >>> itr.params.num_iter = 100
        >>> recon.set_iterate(itr)
        >>> for itr.params.iteration in range(1, 101):
        ...     recon.run_iteration()
    '''

    def __init__(self, num_threads=-1):
        '''Initialize EMCRecon object.'''
        if num_threads <= 0:
            self.num_threads = openmp.omp_get_max_threads()
        else:
            self.num_threads = num_threads
        openmp.omp_set_num_threads(self.num_threads)

        self.mdata = <c_recon.max_data*> calloc(1, sizeof(c_recon.max_data))
        self.mdata.within_openmp = 0

    def set_iterate(self, Iterate itr):
        '''Set the Iterate object for reconstruction.

        Args:
            itr (Iterate): Iterate object with data, model, quat, and params.
        '''
        if itr.iter.dset == NULL:
            print('Set data for iterate first')
            return
        if itr.iter.mod == NULL:
            print('Set model for iterate first')
            return
        if itr.iter.quat == NULL:
            print('Set quat for iterate first')
            return
        if itr.iter.par == NULL:
            print('Set params for iterate first')
            return

        itr.iter.quat.num_rot_p = itr.quat.divide(MPI.COMM_WORLD.rank, MPI.COMM_WORLD.size, itr.iter.par.num_modes, itr.iter.par.nonrot_modes)
        itr.iter.par.rank = MPI.COMM_WORLD.rank
        itr.iter.par.num_proc = MPI.COMM_WORLD.size
        self.mdata.iter = itr.iter

        if itr.stochastic_active():
            MPI.COMM_WORLD.Bcast([<uint8_t[:itr.iter.tot_num_data]>itr.iter.blacklist, MPI.BYTE], 0)

    def calc_frame_counts(self):
        '''Recalculate eager per-frame counts and detector means.'''
        self.iter.calc_frame_counts()

    def calc_sum_fact(self):
        '''Recalculate eager per-frame log-factorial terms.'''
        self.iter.calc_sum_fact()

    def calc_powder(self):
        '''Recalculate eager detector powder images.'''
        self.iter.calc_powder()

    def run_iteration(self, dynamic=False):
        '''Run a single EMC iteration.

        Args:
            dynamic (bool): Use dynamic output. Default False.
        '''
        cdef double likelihood, beta_mean
        cdef double stime = time.time()
        cdef c_iterate.iterate *itr = self.mdata.iter
        cdef bint refreshed

        if itr == NULL:
            print('Set iterate first')
            return

        if itr.par.rtype == c_params.RECON3D:
            c_recon.slice_gen = &c_model.slice_gen3d
            c_recon.slice_merge = &c_model.slice_merge3d
        elif itr.par.rtype == c_params.RECON2D:
            c_recon.slice_gen = &c_model.slice_gen2d
            c_recon.slice_merge = &c_model.slice_merge2d

        cdef long vol = itr.mod.num_modes * itr.mod.vol
        MPI.COMM_WORLD.Bcast([<double[:vol]>itr.mod.model1, MPI.DOUBLE], 0)
        refreshed = self.iter.update_blacklist(finalize=False)
        if self.iter.stochastic_active() and refreshed:
            MPI.COMM_WORLD.Bcast([<uint8_t[:itr.tot_num_data]>itr.blacklist, MPI.BYTE], 0)
        self.iter._finalize_blacklist()
        if itr.par.iteration == 1:
            self.write_log_file_header()

        beta_mean = self.update_beta()
        likelihood = c_recon.maximize(self.mdata)

        if itr.par.rank == 0:
            self.update_model()
        if itr.par.need_scaling == 1 and itr.mod.mtype == c_model.MODEL_3D:
            self.iter.normalize_scale()
        if itr.par.rank == 0:
            self.save_output(dynamic)
            self.update_log_file(likelihood, beta_mean, time.time()-stime)

        if itr.par.verbosity > 0 and itr.par.rank == 0:
            print('Finished iteration', itr.par.iteration, '(%e)' % itr.rms_change)

        MPI.COMM_WORLD.Bcast([<double[:1]>&itr.rms_change, MPI.DOUBLE], 0)
        if isnan(itr.rms_change):
            if itr.par.rank == 0:
                print('rms_change = NAN')
            return False
        return True

    def update_model(self):
        '''Update model after iteration.'''
        cdef long x
        cdef double diff, change
        cdef c_model.model *mod = self.mdata.iter.mod
        cdef c_params.params *param = self.mdata.iter.par
        cdef c_quat.quaternion *quat = self.mdata.iter.quat

        for x in range(mod.num_modes * mod.vol):
            if mod.inter_weight[x] > 0.:
                mod.model2[x] /= mod.inter_weight[x]

        if param.rtype == c_params.RECONRZ or (param.rtype == c_params.RECON2D and param.friedel_sym):
            c_model.symmetrize_friedel2d(mod.model2, mod.inter_weight, mod.num_modes, mod.size)
        elif param.rtype == c_params.RECON3D and quat.icosahedral_flag:
            for x in range(mod.num_modes):
                c_model.symmetrize_icosahedral(&mod.model2[x*mod.vol], &mod.inter_weight[x*mod.vol], mod.size)
        elif param.rtype == c_params.RECON3D and quat.octahedral_flag:
            for x in range(mod.num_modes):
                c_model.symmetrize_octahedral(&mod.model2[x*mod.vol], &mod.inter_weight[x*mod.vol], mod.size)
        elif param.rtype == c_params.RECON3D:
            for x in range(mod.num_modes):
                c_model.symmetrize_friedel(&mod.model2[x*mod.vol], &mod.inter_weight[x*mod.vol], mod.size)
        
        if param.axial_sym > 1:
            for x in range(mod.num_modes):
                c_model.symmetrize_axial(&mod.model2[x*mod.vol], &mod.inter_weight[x*mod.vol], mod.size, param.axial_sym)

        change = 0
        for x in range(mod.num_modes * mod.vol):
            diff = mod.model2[x] - mod.model1[x]
            change += diff**2
            if param.alpha > 0.:
                mod.model1[x] = param.alpha * mod.model1[x] + (1-param.alpha) * mod.model2[x]
            else:
                mod.model1[x] = mod.model2[x]

        self.mdata.iter.rms_change = sqrt(change / mod.num_modes / mod.vol)

    def update_beta(self):
        '''Update beta values for current iteration.

        Returns:
            float: Mean beta value.
        '''
        cdef int d
        cdef double factor, beta_mean = 0.
        cdef c_iterate.iterate *itr = self.mdata.iter

        factor = itr.par.beta_jump ** ((itr.par.iteration-1) // itr.par.beta_period)
        factor *= itr.par.beta_factor

        for d in range(itr.tot_num_data):
            if itr.blacklist[d] == 1:
                continue
            itr.beta[d] = itr.beta_start[d] * factor
            if itr.beta[d] > 1.:
                itr.beta[d] = 1.
            beta_mean += itr.beta[d]
        beta_mean /= (itr.tot_num_data - itr.num_blacklist)

        return beta_mean

    def free(self):
        '''Free allocated reconstruction memory.'''
        if self.mdata == NULL or self.mdata.iter == NULL:
            return
        c_recon.free_max_data(self.mdata)
        free(self.mdata)
        self.mdata = NULL

    def write_log_file_header(self):
        '''Write header to log file.'''
        cdef c_iterate.iterate *itr = self.mdata.iter
        cdef c_model.model *mod = itr.mod
        cdef c_params.params *param = itr.par
        cdef c_quat.quaternion *quat = itr.quat
        cdef c_dataset.dataset *frames = itr.dset

        with open((<bytes> param.log_fname).decode(), 'w') as fp:
            fp.write('Cryptotomography with the EMC algorithm using MPI+OpenMP\n\n')
            fp.write('Data parameters:\n')
            if itr.num_blacklist == 0:
                fp.write('\tnum_data = %d\n\tmean_count = %f\n\n' % (itr.tot_num_data, itr.mean_count[0]))
            else:
                fp.write('\tnum_data = %d/%d\n\tmean_count = %f\n\n' % (itr.tot_num_data-itr.num_blacklist, itr.tot_num_data, itr.mean_count[0]))
            fp.write('System size:\n')
            fp.write('\tnum_rot = %d\n\tnum_pix = %d/%d\n\t' % (quat.num_rot, (self.iter.dets[0].raw_mask==0).sum(), itr.det.num_pix))
            if param.rtype == c_params.RECON3D:
                fp.write('system_volume = %d X %ld X %ld X %ld\n\n' % (mod.num_modes, mod.size, mod.size, mod.size))
            elif param.rtype == c_params.RECON2D or param.rtype == c_params.RECONRZ:
                fp.write('system_volume = %d X %ld X %ld\n\n' % (mod.num_modes, mod.size, mod.size))
            fp.write('Reconstruction parameters:\n')
            fp.write('\tnum_threads = %d\n\tnum_proc = %d\n\talpha = %.6f\n\tbeta = %.6f\n\tneed_scaling = %s' % (
                    self.num_threads,
                    param.num_proc,
                    param.alpha,
                    itr.beta_start[0],
                    'yes' if param.need_scaling == 1 else 'no'))
            fp.write('\n\nIter\ttime\trms_change\tinfo_rate\tlog-likelihood\tnum_rot\tbeta\n')
            fp.close()

    def update_log_file(self, double likelihood, double beta, double iter_time):
        '''Append iteration results to log file.

        Args:
            likelihood (float): Log-likelihood value.
            beta (float): Beta value.
            iter_time (float): Iteration time in seconds.
        '''
        cdef c_iterate.iterate *itr = self.mdata.iter
        cdef c_params.params *param = itr.par
        cdef double prev_beta

        log_fname = (<bytes>param.log_fname).decode()
        if (param.iteration != param.start_iter and
                (param.iteration - 1) % param.beta_period != 0 and os.path.exists(log_fname)):
            with open(log_fname, 'r') as fp:
                lines = fp.readlines()
                if len(lines) > 0:
                    fields = lines[len(lines)-1].split()
                    if len(fields) == 7 and fields[0].isdigit():
                        prev_beta = float(fields[6])
                        beta = prev_beta

        fp = open(log_fname, 'a')
        fp.write('%d\t' % param.iteration)
        fp.write('%4.2f\t' % iter_time)
        fp.write('%1.4e\t%f\t%.6e\t%-7d\t%f\n' % (itr.rms_change, itr.mutual_info, likelihood, itr.quat.num_rot, beta))
        fp.close()

    def save_output(self, dynamic=False):
        '''Save iteration output to HDF5 file.

        Args:
            dynamic (bool): Use dynamic output filename. Default False.
        '''
        itr = self.iter
        param = itr.params
        if not os.path.exists(param.output_folder):
            os.makedirs(param.output_folder)
        if dynamic:
            out_fname = '%s/output_dynamic.h5' % param.output_folder
        else:
            out_fname = '%s/output_%.3d.h5' % (param.output_folder, param.iteration)

        with h5py.File(out_fname, 'w') as f:
            for attr in dir(itr.params):
                if not (attr[0] == '_' or callable(getattr(itr.params, attr))):
                    f['params/' + attr] = getattr(itr.params, attr)
            f['params/det_fnames'] = [det.fname for det in itr.dets]
            f['params/emc_fnames'] = [emc.fname for emc in itr.data]
            f['params/mtype'] = itr.model.mtype
            f['params/num_rot'] = itr.quat.num_rot
            if itr.params.rtype == 'RECON3D':
                f['params/num_div'] = itr.quat.num_div

            f['orientations'] = self.rmax
            f['mutual_info'] = self.info
            f['likelihood'] = self.likelihood
            if self.mdata.iter.mod.num_modes > 1:
                f['occupancies'] = self.quat_norm

            f['intens'] = itr.model.model1
            f['inter_weight'] = itr.model.inter_weight
            if param.need_scaling != 0:
                f['scale'] = itr.scale
            if itr.stochastic_active():
                f['blacklist'] = itr.blacklist
                if param.coverage_bias != 0.:
                    f['last_selected'] = itr.last_selected

            f['quaternions'] = itr.quat.quats

            if param.save_prob != 0:
                f['probabilities/num_rot'] = [self.mdata.iter.quat.num_rot]
                num_data = itr.tot_num_data

                dtype = h5py.special_dtype(vlen='i4')
                place_dset = f.create_dataset('probabilities/place', (num_data,), dtype=dtype)
                place_dset[:] = self.place_prob

                dtype = h5py.special_dtype(vlen='f8')
                prob_dset = f.create_dataset('probabilities/prob', (num_data,), dtype=dtype)
                prob_dset[:] = self.prob

    def save_initial_model(self):
        '''Save initial model guess to output_000.h5.'''
        itr = self.iter
        param = itr.params
        if not os.path.exists(param.output_folder):
            os.makedirs(param.output_folder)
        out_fname = '%s/output_000.h5' % param.output_folder

        with h5py.File(out_fname, 'w') as f:
            f['intens'] = itr.model.model1
            if param.need_scaling != 0:
                f['scale'] = itr.scale

    @property
    def num_threads(self):
        '''Number of OpenMP threads.'''
        return self.num_threads
    @num_threads.setter
    def num_threads(self, int val):
        '''Set number of OpenMP threads.'''
        if val <= 0:
            self.num_threads = openmp.omp_get_max_threads()
        else:
            self.num_threads = val

    @property
    def iter(self):
        '''Associated Iterate object.'''
        if self.mdata.iter == NULL:
            return
        retval = Iterate()
        retval.free()
        retval.iter = self.mdata.iter
        return retval

    @property
    def refinement(self):
        '''Refinement mode flag.'''
        return self.mdata.refinement
    @refinement.setter
    def refinement(self, int val):
        '''Set refinement flag.'''
        self.mdata.refinement = val
    @property
    def within_openmp(self):
        '''OpenMP threading flag.'''
        return self.mdata.within_openmp
    @within_openmp.setter
    def within_openmp(self, int val):
        '''Set OpenMP flag.'''
        self.mdata.within_openmp = val

    @property
    def model(self):
        '''Per-thread model (OpenMP private).'''
        return np.asarray(<double[:self.mdata.iter.mod.vol]>self.mdata.model).reshape(3*(self.mdata.iter.mod.size,)) if self.mdata.model != NULL else None
    @property
    def weight(self):
        '''Per-thread weights (OpenMP private).'''
        return np.asarray(<double[:self.mdata.iter.mod.vol]>self.mdata.weight).reshape(3*(self.mdata.iter.mod.size,)) if self.mdata.weight != NULL else None
    @property
    def all_views(self):
        '''All detector views (OpenMP private).'''
        cdef int num_pix
        if self.mdata.all_views == NULL:
            return
        else:
            retval = []
            for i in range(self.mdata.iter.num_det):
                num_pix = self.mdata.iter.det[i].num_pix
                retval.append(<double[:num_pix]> self.mdata.all_views[i])
            return retval
    @property
    def psum_r(self):
        '''Per-rotation sums (OpenMP private).'''
        return np.asarray(<double[:self.mdata.iter.quat.num_rot_p]>self.mdata.psum_r) if self.mdata.psum_r != NULL else None
    @property
    def psum_d(self):
        '''Per-data sums (OpenMP private).'''
        return np.asarray(<double[:self.mdata.iter.tot_num_data]>self.mdata.psum_d) if self.mdata.psum_d != NULL else None

    @property
    def max_exp(self):
        '''Max exposure values (shared among threads).'''
        return np.asarray(<double[:self.mdata.iter.tot_num_data]>self.mdata.max_exp) if self.mdata.max_exp != NULL else None
    @property
    def p_norm(self):
        '''Probability normalization (shared among threads).'''
        return np.asarray(<double[:self.mdata.iter.tot_num_data]>self.mdata.p_norm) if self.mdata.p_norm != NULL else None
    @property
    def u(self):
        '''Intermediate values (shared among threads).'''
        if self.mdata.u == NULL:
            return
        retval = []
        for i in range(self.mdata.iter.num_det):
            retval.append(np.asarray(<double[:self.mdata.iter.quat.num_rot_p]>self.mdata.u[i]))
        return retval
    @property
    def offset_prob(self):
        '''Probability offsets (shared among threads).'''
        return np.asarray(<int[:self.num_threads*self.mdata.iter.tot_num_data]>self.mdata.offset_prob).reshape(self.num_threads, -1) if self.mdata.offset_prob != NULL else None

    @property
    def max_exp_p(self):
        '''Per-thread max exposure (both private and shared).'''
        return np.asarray(<double[:self.mdata.iter.tot_num_data]>self.mdata.max_exp_p) if self.mdata.max_exp_p != NULL else None
    @property
    def info(self):
        '''Per-frame mutual information.'''
        return np.asarray(<double[:self.mdata.iter.tot_num_data]>self.mdata.info) if self.mdata.info != NULL else None
    @property
    def likelihood(self):
        '''Per-frame log-likelihood.'''
        return np.asarray(<double[:self.mdata.iter.tot_num_data]>self.mdata.likelihood) if self.mdata.likelihood != NULL else None
    @property
    def rmax(self):
        '''Most likely rotation index per frame.'''
        return np.asarray(<int[:self.mdata.iter.tot_num_data]>self.mdata.rmax) if self.mdata.rmax != NULL else None
    @property
    def quat_norm(self):
        '''Quaternion normalization factors.'''
        return np.asarray(<double[:self.mdata.iter.tot_num_data*self.mdata.iter.mod.num_modes]>self.mdata.quat_norm).reshape(-1, self.mdata.iter.mod.num_modes) if self.mdata.quat_norm != NULL else None
    @property
    def num_prob(self):
        '''Number of probabilities per frame.'''
        return np.asarray(<int[:self.mdata.iter.tot_num_data]>self.mdata.num_prob) if self.mdata.num_prob != NULL else None
    @property
    def prob(self):
        '''Probability arrays per frame.'''
        return [np.asarray(<double[:self.mdata.num_prob[d]]>self.mdata.prob[d]) if self.mdata.num_prob[d] > 0 else None for d in range(self.mdata.iter.tot_num_data)] if self.mdata.prob != NULL else None
    @property
    def place_prob(self):
        '''Probability place indices per frame.'''
        return [np.asarray(<int[:self.mdata.num_prob[d]]>self.mdata.place_prob[d]) for d in range(self.mdata.iter.tot_num_data)] if self.mdata.place_prob != NULL else None


def main():
    '''Command-line entry point for EMC reconstruction.'''
    import argparse

    parser = argparse.ArgumentParser(description='Cryptotomography with the EMC algorithm using MPI+OpenMP')
    parser.add_argument('niter', help='Number of iterations', type=int)
    parser.add_argument('-c', '--config_file', help='Path to configuration file. Default: config.ini', default='config.ini')
    parser.add_argument('-r', '--resume', help='Whether to resume from last iteration', action='store_true')
    parser.add_argument('-t', '--threads', help='Number of OpenMP threads per MPI process', type=int, default=-1)
    args = parser.parse_args()

    recon = EMCRecon(args.threads)
    stime = time.time()

    itr = Iterate(args.config_file, resume=args.resume)
    itr.params.num_iter = args.niter
    recon.set_iterate(itr)

    st = itr.params.start_iter
    en = itr.params.start_iter + itr.params.num_iter
    
    # Save initial model guess before first iteration
    if st == 1 and itr.params.rank == 0:
        recon.save_initial_model()
    if itr.params.rank == 0:
        print('Finished setup (%f)' % (time.time() - stime))
    
    for itr.params.iteration in range(st, en):
        if not recon.run_iteration():
            break
    if itr.params.verbosity > 0 and itr.params.rank == 0:
        print('Finished all iterations')
