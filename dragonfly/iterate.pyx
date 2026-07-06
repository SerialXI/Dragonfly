import sys
import os.path as op
import re
import numpy as np
import h5py
import configparser
from .utils.py_src.read_config import MyConfigParser

from libc.stdlib cimport malloc, calloc, free, atoi
from libc.string cimport memcpy
cimport numpy as np
from . cimport iterate as c_iterate
from . cimport params as c_params
from . cimport emcfile as c_emcfile
from .iterate cimport Iterate
from .detector cimport CDetector, detector
from .model cimport Model
from .emcfile cimport CDataset, dataset
from .quaternion cimport Quaternion
from .params cimport EMCParams

cdef class Iterate:
    '''Class managing state for each EMC iteration.

    Handles loading configuration, data, model, and quaternion data for
    iterative reconstruction.

    Args:
        config_fname (str): Path to configuration file. Default ''.
        section_name (str): Section name. Default 'emc'.
        resume (bool): Resume from last iteration. Default False.
    '''

    def __init__(self, config_fname='', section_name='emc', resume=False):
        '''Initialize Iterate object.'''
        self.iter = <c_iterate.iterate*> calloc(1, sizeof(c_iterate.iterate))
        self.iter.par = <c_params.params*> calloc(1, sizeof(c_params.params))
        if config_fname != '':
            self.from_config(config_fname, section_name, resume=resume)

    def from_config(self, config_fname, section_name='emc', resume=False):
        '''Load configuration and initialize from file.

        Args:
            config_fname (str): Path to configuration file.
            section_name (str): Section name. Default 'emc'.
            resume (bool): Resume from last iteration. Default False.
        '''
        param = EMCParams()
        param.from_config(config_fname, section_name)
        self.set_params(param)
        rtypes = {str(c_params.RECON3D): '3d',
                  str(c_params.RECON2D): '2d',
                  str(c_params.RECONRZ): 'rz'}
        rtype = rtypes[str(self.iter.par.rtype)]
        print('Doing %s recon from %s' % (rtype.upper(), config_fname))

        config_folder = op.dirname(config_fname)
        config = MyConfigParser()
        config.read(config_fname)

        det_fname = config.get(section_name, 'in_detector_file', fallback=None)
        det_flist = config.get(section_name, 'in_detector_list', fallback=None)
        if det_fname is not None and det_flist is not None:
            raise ValueError('Both in_detector_file and in_detector_list specified. Pick one.')
        elif det_fname is not None:
            if ':::' in det_fname:
                sname, oname = re.split(':::', det_fname)
                det_fname = config.get(sname, oname)
            dets = [CDetector(op.join(config_folder, det_fname), norm=True, rtype=rtype)]
        elif det_flist is not None:
            if ':::' in det_flist:
                sname, oname = re.split(':::', det_flist)
                det_flist = config.get(sname, oname)
            fptr = open(det_flist, 'r')
            dets = [CDetector(op.join(config_folder, line.strip()), norm=True, rtype=rtype) for line in fptr.readlines()]
            fptr.close()
        else:
            raise ValueError('Need either in_detector_file or in_detector_list.')

        ph_fname = config.get(section_name, 'in_photons_file', fallback=None)
        ph_flist = config.get(section_name, 'in_photons_list', fallback=None)
        if ph_fname is not None and ph_flist is not None:
            raise ValueError('Both in_photons_file and in_photons_list specified. Pick one.')
        elif ph_fname is not None:
            if len(dets) > 1:
                print('WARNING: Multiple detectors but only one photons file. Using first detector')
            if ':::' in ph_fname:
                sname, oname = re.split(':::', ph_fname)
                ph_fname = config.get(sname, oname)
            frames = CDataset(op.join(config_folder, ph_fname), dets[0], lazy=bool(param.lazy_data))
        elif ph_flist is not None:
            if ':::' in ph_flist:
                sname, oname = re.split(':::', ph_flist)
                ph_flist = config.get(sname, oname)

            fptr = open(ph_flist, 'r')
            fnames = [op.join(config_folder, line.strip()) for line in fptr.readlines()]
            fptr.close()

            frames = CDataset(fnames[0], dets[0], lazy=bool(param.lazy_data))
            if len(dets) > 1:
                [frames.append(CDataset(op.join(config_folder, fnames[i]), dets[i], lazy=bool(param.lazy_data))) for i in range(1, len(fnames))]
            else:
                [frames.append(CDataset(op.join(config_folder, fnames[i]), dets[0], lazy=bool(param.lazy_data))) for i in range(1, len(fnames))]
        else:
            raise ValueError('Need either in_photons_file or in_photons_list.')

        background_file = config.get_filename(section_name, 'background_file', fallback=None)
        background_list = config.get_filename(section_name, 'background_list', fallback=None)
        if background_file is not None and background_list is not None:
            raise ValueError('Both background_file and background_list specified. Pick one.')

        if background_list is not None:
            bg_list_dir = op.dirname(background_list)
            with open(background_list, 'r') as f:
                bg_lines = [line.strip() for line in f if line.strip()]
            if len(bg_lines) != len(dets):
                raise ValueError('Background list has %d entries but there are %d unique detectors' %
                                 (len(bg_lines), len(dets)))
            for i, bg_fname in enumerate(bg_lines):
                dets[i].parse_background(op.join(bg_list_dir, bg_fname))
        elif background_file is not None:
            if len(dets) > 1:
                print('WARNING: Multiple detectors and single background file. '
                      'Assuming same background for all detectors')
            for det in dets:
                det.parse_background(background_file)

        self.set_data(frames)

        quat = Quaternion()
        quat.from_config(config_fname, section_name)
        self.set_quat(quat)

        try:
            model_file = config.get_filename(section_name, 'start_model_file')
        except configparser.NoOptionError:
            model_file = ''

        try:
            scale_file = config.get_filename(section_name, 'scale_file')
        except configparser.NoOptionError:
            scale_file = ''

        try:
            bgscale_file = config.get_filename(section_name, 'bgscale_file')
        except configparser.NoOptionError:
            bgscale_file = ''
        if resume:
            try:
                fp = open(param.log_fname, 'r')
                lines = fp.readlines()
                fp.close()
            except FileNotFoundError:
                print('No log file found to resume reconstruction. Starting from scratch.')
            else:
                fields = lines[len(lines)-1].split() if len(lines) > 0 else []
                try:
                    last_iter = int(fields[0])
                except (IndexError, ValueError):
                    print('No completed iterations found in log file. Starting from scratch and overwriting log file.')
                else:
                    self.params.start_iter = last_iter + 1
                    model_file = op.join(param.output_folder, 'output_%.3d.h5' % last_iter)
                    scale_file = op.join(param.output_folder, 'output_%.3d.h5' % last_iter)
                    print('Resuming from iteration', self.params.start_iter)

        if self.iter.par.need_scaling == 1:
            self.parse_scale(scale_file)
            self.parse_scale(bgscale_file, bg=True)
        sel_string = config.get(section_name, 'selection', fallback=None)
        try:
            blacklist_file = config.get_filename(section_name, 'blacklist_file')
        except configparser.NoOptionError:
            blacklist_file = ''
        self.parse_blacklist(blacklist_file, sel_string)
        self.refresh_data_stats()
        if (self.iter.par.need_scaling == 1 and self.iter.par.data_fraction < 1.
                and scale_file == ''):
            self.initialize_stochastic_scale()
        if self.iter.par.data_fraction < 1. and self.iter.par.coverage_bias != 0.:
            self.initialize_last_selected(scale_file)

        if self.iter.par.lazy_data:
            self.params.iteration = self.params.start_iter
            if not self.load_previous_blacklist(model_file):
                self.update_blacklist(force=True)
            self.load_active_data(update_powder=True)

        qmax = max([det.qmax() for det in dets])
        model = Model(self.calculate_size(qmax), self.iter.par.num_modes, rtype)
        model_mean = self.mean_count[0] / (self.dets[0].raw_mask==0).sum() * 2.
        model.allocate(model_file, model_mean, fixed_seed=bool(param.fixed_seed))
        self.set_model(model)

        if self.iter.par.lazy_data == 0:
            if self.iter.par.beta_config < 0.:
                self.calc_beta()
            else:
                self.calc_beta(self.iter.par.beta_config)

    def set_model(self, Model model):
        '''Set the model.

        Args:
            model (Model): Model object.
        '''
        self.iter.mod = model.mod

    def update_data(self, CDataset in_dset):
        '''Update the data with a new dataset.

        Args:
            in_dset (CDataset): New dataset to use.
        '''
        cdef int total = 0
        cdef dataset *curr = in_dset.dset

        self.iter.num_dfiles = 0
        while curr != NULL:
            curr.num_offset = total
            total += curr.num_data
            self.iter.num_dfiles += 1
            curr = curr.next
        self.iter.tot_num_data = total
        print('New tot_num_data =', self.iter.tot_num_data)

        self.iter.dset = in_dset.dset

        free(self.iter.det_mapping)
        self.iter.det_mapping = NULL
        self.iter.det_mapping = <int*> calloc(self.iter.num_dfiles, sizeof(int))

        free(self.iter.blacklist)
        self.iter.blacklist = <uint8_t*> calloc(self.iter.tot_num_data, sizeof(uint8_t))
        free(self.iter.static_blacklist)
        self.iter.static_blacklist = <uint8_t*> calloc(self.iter.tot_num_data, sizeof(uint8_t))
        self.iter.num_blacklist = 0

        free(self.iter.fcounts)
        self.iter.fcounts = NULL
        if self.iter.par.lazy_data == 0:
            self.calc_frame_counts()
            self._blacklist_zeros()

        free(self.iter.last_selected)
        self.iter.last_selected = NULL

        free(self.iter.sum_fact)
        self.iter.sum_fact = NULL
        if self.iter.par.lazy_data == 0:
            self.calc_sum_fact()

        if self.iter.par.need_scaling == 1:
            free(self.iter.scale)
            self.iter.scale = NULL
            self.parse_scale('')

            free(self.iter.bgscale)
            self.parse_scale('', bg=True)

        free(self.iter.beta)
        free(self.iter.beta_start)
        if self.iter.par.lazy_data == 0:
            self.calc_beta()

        if self.iter.par.lazy_data == 0:
            self.calc_powder()

    def set_data(self, CDataset in_dset):
        '''Set the data with a new dataset.

        Args:
            in_dset (CDataset): Dataset to use.
        '''
        cdef int total = 0
        cdef dataset *curr = in_dset.dset

        if self.iter.tot_num_data != 0:
            raise ValueError('Please use update_data() to change dataset')

        while curr != NULL:
            curr.num_offset = total
            total += curr.num_data
            curr = curr.next
        self.iter.tot_num_data = total

        self.iter.dset = in_dset.dset

        self._gen_detlist()

        self.iter.blacklist = <uint8_t*> calloc(self.iter.tot_num_data, sizeof(uint8_t))
        self.iter.static_blacklist = <uint8_t*> calloc(self.iter.tot_num_data, sizeof(uint8_t))
        self.iter.num_blacklist = 0
        if self.iter.par.lazy_data == 0:
            self.calc_frame_counts()
            self.calc_sum_fact()
            self.calc_powder()

    def refresh_data_stats(self, update_powder=True):
        '''Recalculate eager data statistics after changing the blacklist.'''
        if self.iter.par.lazy_data == 0:
            self.calc_frame_counts()
            self._blacklist_zeros()
            self.calc_sum_fact()
            if update_powder:
                self.calc_powder()

    def calc_frame_counts(self):
        '''Calculate per-frame counts.'''
        c_iterate.calc_frame_counts(self.iter)

    def calc_frame_counts_partial(self, int rank, int num_proc):
        '''Calculate this rank's partial per-frame counts.'''
        cdef int[:] num_data

        num_data = np.zeros(self.iter.num_det, dtype=np.int32)
        c_iterate.calc_frame_counts_partial(self.iter, rank, num_proc, &num_data[0])
        return np.asarray(num_data)

    def _blacklist_zeros(self):
        '''Add active zero-photon frames to the static blacklist.'''
        cdef int dset = 0
        cdef int d, curr_d, detn
        cdef int num_added = 0
        cdef int[:] active_count
        cdef int[:] zero_count
        cdef dataset *curr = self.iter.dset

        if self.iter.fcounts == NULL:
            raise AttributeError('Need frame counts before blacklisting zero-count frames')
        if self.iter.mean_count == NULL:
            raise AttributeError('Need mean counts before blacklisting zero-count frames')
        if self.iter.blacklist == NULL or self.iter.static_blacklist == NULL:
            raise AttributeError('Blacklist has not been initialized')

        active_count = np.zeros(self.iter.num_det, dtype=np.int32)
        zero_count = np.zeros(self.iter.num_det, dtype=np.int32)

        while curr != NULL:
            detn = self.iter.det_mapping[dset]
            for curr_d in range(curr.num_data):
                d = curr.num_offset + curr_d
                if self.iter.blacklist[d] == 0:
                    active_count[detn] += 1
                    if self.iter.fcounts[d] == 0:
                        self.iter.blacklist[d] = 1
                        self.iter.static_blacklist[d] = 1
                        zero_count[detn] += 1
                        num_added += 1
            dset += 1
            curr = curr.next

        if num_added == 0:
            return 0

        for detn in range(self.iter.num_det):
            if zero_count[detn] > 0 and zero_count[detn] < active_count[detn]:
                self.iter.mean_count[detn] *= (<double> active_count[detn] /
                                               (active_count[detn] - zero_count[detn]))

        if self.iter.par.rank == 0:
            sys.stderr.write('WARNING: Auto-blacklisted %d zero-photon frame(s)\n' % num_added)
        self._finalize_blacklist()
        return num_added

    def calc_sum_fact(self):
        '''Calculate per-frame log-factorial terms.'''
        c_iterate.calc_sum_fact(self.iter)

    def calc_sum_fact_partial(self, int rank, int num_proc):
        '''Calculate this rank's partial per-frame log-factorial terms.'''
        c_iterate.calc_sum_fact_partial(self.iter, rank, num_proc)

    def calc_powder(self):
        '''Calculate detector powder images.'''
        c_iterate.calc_powder(self.iter)

    def calc_powder_partial(self, int rank, int num_proc):
        '''Calculate this rank's partial detector powder images.'''
        cdef int[:] nframes

        nframes = np.zeros(self.iter.num_det, dtype=np.int32)
        c_iterate.calc_powder_partial(self.iter, rank, num_proc, &nframes[0])
        return np.asarray(nframes)

    def set_quat(self, Quaternion quaternion):
        '''Set the quaternion orientations.

        Args:
            quaternion (Quaternion): Quaternion object.
        '''
        self.iter.quat = quaternion.quat

    def set_params(self, EMCParams param):
        '''Set the parameters.

        Args:
            param (EMCParams): Parameters object.
        '''
        self.iter.par = param.par

    def load_previous_blacklist(self, fname):
        '''Load held dynamic blacklist when resuming within a data fraction period.'''
        cdef uint8_t[:] blist
        cdef int d

        if (not self.stochastic_active() or
                (self.iter.par.iteration - 1) % self.iter.par.data_fraction_period == 0):
            return False
        if not (op.isfile(fname) and h5py.is_hdf5(fname)):
            return False

        with h5py.File(fname, 'r') as f:
            if 'blacklist' not in f:
                return False
            arr = f['blacklist'][:].astype('u1')
        if arr.shape[0] != self.iter.tot_num_data:
            raise ValueError('Mismatched number of frames in saved blacklist')
        blist = arr
        memcpy(self.iter.blacklist, &blist[0], self.iter.tot_num_data*sizeof(uint8_t))
        self.iter.num_blacklist = 0
        for d in range(self.iter.tot_num_data):
            if self.iter.blacklist[d] != 0:
                self.iter.num_blacklist += 1
        return True

    def load_active_data(self, update_powder=False):
        '''Load currently active sparse binary frames for lazy_data mode.'''
        if self.iter.par.lazy_data == 0:
            return
        self.load_active_frames()

        self.calc_frame_counts()
        self._blacklist_zeros()
        self.calc_sum_fact()
        if update_powder:
            self.calc_powder()
        if self.iter.beta != NULL:
            free(self.iter.beta)
            self.iter.beta = NULL
        if self.iter.beta_start != NULL:
            free(self.iter.beta_start)
            self.iter.beta_start = NULL
        if self.iter.par.beta_config < 0.:
            self.calc_beta()
        else:
            self.calc_beta(self.iter.par.beta_config)

    def load_active_frames(self):
        '''Load currently active sparse binary frames without recalculating stats.'''
        cdef dataset *curr = self.iter.dset

        if self.iter.par.lazy_data == 0:
            return
        while curr != NULL:
            if c_emcfile.load_active_frames(curr, self.iter.blacklist) != 0:
                raise ValueError('Could not load active frames')
            curr = curr.next

    def parse_scale(self, fname, bg=False):
        '''Parse or initialize scale factors.

        Args:
            fname (str): Path to scale file or empty string for defaults.
            bg (bool): Parse background scale. Default False.
        '''
        if not op.isfile(fname):
            if self.iter.tot_num_data == 0:
                raise AttributeError('Need tot_num_data to initialize scale factors')
            scale = np.ones(self.iter.tot_num_data)
        else:
            with h5py.File(fname, 'r') as f:
                scale = f['scale'][:]

        cdef double[:] scale_view = scale

        if self.iter.tot_num_data == 0:
            self.iter.tot_num_data = scale_view.shape[0]
        elif self.iter.tot_num_data < scale_view.shape[0]:
            print('More scale factors than required. Taking only first', self.iter.tot_num_data)
        elif self.iter.tot_num_data > scale_view.shape[0]:
            print('Insufficient scale factors in file. Setting rest to unity')
            scale = np.append(scale, np.ones(self.iter.tot_num_data - scale_view.shape[0]))
            scale_view = scale

        if bg:
            self.iter.bgscale = <double*> malloc(self.iter.tot_num_data * sizeof(double))
            memcpy(self.iter.bgscale, &scale_view[0], self.iter.tot_num_data * sizeof(double))
        else:
            self.iter.scale = <double*> malloc(self.iter.tot_num_data * sizeof(double))
            memcpy(self.iter.scale, &scale_view[0], self.iter.tot_num_data * sizeof(double))

    def parse_blacklist(self, fname, sel_string=None, refresh=False):
        '''Generate or update blacklist.

        Blacklist file contains one number (0 or 1) per frame indicating
        whether the frame is blacklisted (1) or good (0).

        Args:
            fname (str): Path to blacklist file.
            sel_string (str): Selection string ('odd_only' or 'even_only'). Default None.
            refresh (bool): Refresh the blacklist. Default False.
        '''
        if self.iter.tot_num_data == 0:
            raise AttributeError('Need to define tot_num_data before generating blacklist')

        cdef int d
        cdef uint8_t curr
        cdef uint8_t[:] arr
        self.iter.num_blacklist = 0
        if refresh:
            if self.iter.blacklist != NULL:
                free(self.iter.blacklist)
            self.iter.blacklist = <uint8_t*> calloc(self.iter.tot_num_data, sizeof(uint8_t))
            if self.iter.static_blacklist != NULL:
                free(self.iter.static_blacklist)
            self.iter.static_blacklist = <uint8_t*> calloc(self.iter.tot_num_data, sizeof(uint8_t))
        else:
            if self.iter.blacklist == NULL:
                self.iter.blacklist = <uint8_t*> calloc(self.iter.tot_num_data, sizeof(uint8_t))
            if self.iter.static_blacklist == NULL:
                self.iter.static_blacklist = <uint8_t*> calloc(self.iter.tot_num_data, sizeof(uint8_t))

        if op.isfile(fname):
            arr = np.loadtxt(fname, dtype='u1')
            if arr.shape[0] != self.iter.tot_num_data:
                raise ValueError('Mismatched number of frames in blacklist file')
            memcpy(self.iter.static_blacklist, &arr[0], self.iter.tot_num_data*sizeof(uint8_t))

        if sel_string is None:
            pass
        elif sel_string == 'odd_only':
            curr = 0
            for d in range(self.iter.tot_num_data):
                if self.iter.static_blacklist[d] == 0:
                    self.iter.static_blacklist[d] = curr
                    curr = 1 - curr
        elif sel_string == 'even_only':
            curr = 1
            for d in range(self.iter.tot_num_data):
                if self.iter.static_blacklist[d] == 0:
                    self.iter.static_blacklist[d] = curr
                    curr = 1 - curr
        else:
            raise ValueError('Unrecognized sel_string, %s' % sel_string)

        memcpy(self.iter.blacklist, self.iter.static_blacklist, self.iter.tot_num_data*sizeof(uint8_t))
        for d in range(self.iter.tot_num_data):
            if self.iter.blacklist[d] != 0:
                self.iter.num_blacklist += 1

    def initialize_stochastic_scale(self):
        '''Mark fresh stochastic-EMC scale factors as not yet estimated.'''
        cdef int d
        if self.iter.scale == NULL:
            return
        if self.iter.static_blacklist == NULL:
            raise AttributeError('Static blacklist has not been initialized')

        for d in range(self.iter.tot_num_data):
            self.iter.scale[d] = -1.

    def initialize_last_selected(self, fname=''):
        '''Initialize or load stochastic coverage sampling history.'''
        cdef int[:] last_view
        cdef int d

        if self.iter.last_selected != NULL:
            free(self.iter.last_selected)
        self.iter.last_selected = <int*> malloc(self.iter.tot_num_data * sizeof(int))

        if op.isfile(fname) and h5py.is_hdf5(fname):
            with h5py.File(fname, 'r') as f:
                if 'last_selected' in f:
                    last = f['last_selected'][:].astype('i4')
                else:
                    last = np.full(self.iter.tot_num_data, -1, dtype='i4')
        else:
            last = np.full(self.iter.tot_num_data, -1, dtype='i4')

        if last.shape[0] != self.iter.tot_num_data:
            raise ValueError('Mismatched number of frames in last_selected')
        last_view = last
        memcpy(self.iter.last_selected, &last_view[0], self.iter.tot_num_data*sizeof(int))

    def _selection_prob(self, good):
        '''Calculate per-frame stochastic selection probabilities.'''
        if self.iter.par.coverage_bias == 0.:
            return self.iter.par.data_fraction

        if self.iter.last_selected == NULL:
            self.initialize_last_selected()

        last_selected = np.asarray(<int[:self.iter.tot_num_data]>self.iter.last_selected)
        max_age = 5. / self.iter.par.data_fraction
        ages = (self.iter.par.iteration - last_selected[good]).astype('f8')
        ages[last_selected[good] < 0] = max_age
        ages = np.minimum(ages, max_age)
        weights = 1. + self.iter.par.coverage_bias * ages / max_age
        target = self.iter.par.data_fraction * good.size

        lo, hi = 0., 1.
        while np.minimum(1., hi * weights).sum() < target:
            hi *= 2.
        for _ in range(32):
            mid = 0.5 * (lo + hi)
            if np.minimum(1., mid * weights).sum() < target:
                lo = mid
            else:
                hi = mid
        return np.minimum(1., hi * weights)

    def stochastic_active(self):
        '''Whether stochastic frame subsets are active for the current iteration.'''
        return (self.iter.par.data_fraction < 1. and
                (self.iter.par.data_fraction_until <= 0 or
                 self.iter.par.iteration <= self.iter.par.data_fraction_until))

    def _update_last_selected(self, uint8_t[:] active_view):
        '''Update rank-0 coverage sampling history.'''
        if (not self.stochastic_active() or self.iter.par.coverage_bias == 0.
                or self.iter.par.rank != 0):
            return

        if self.iter.last_selected == NULL:
            self.initialize_last_selected()
        last_selected = np.asarray(<int[:self.iter.tot_num_data]>self.iter.last_selected)
        last_selected[np.asarray(active_view) == 0] = self.iter.par.iteration

    def update_blacklist(self, force=False, finalize=True):
        '''Update active blacklist for the current iteration.'''
        cdef bint refresh
        cdef bint stochastic
        cdef uint8_t[:] active_view
        cdef uint8_t[:] static_view

        if self.iter.blacklist == NULL or self.iter.static_blacklist == NULL:
            raise AttributeError('Blacklist has not been initialized')

        active_view = <uint8_t[:self.iter.tot_num_data]>self.iter.blacklist
        static_view = <uint8_t[:self.iter.tot_num_data]>self.iter.static_blacklist
        stochastic = self.stochastic_active()
        refresh = (force or not stochastic or
                   (self.iter.par.iteration - 1) % self.iter.par.data_fraction_period == 0)
        if refresh:
            memcpy(self.iter.blacklist, self.iter.static_blacklist, self.iter.tot_num_data*sizeof(uint8_t))

        if stochastic and refresh:
            if self.iter.par.rank == 0:
                active = np.asarray(active_view)
                static = np.asarray(static_view)
                good = np.flatnonzero(static == 0)
                if good.size > 0:
                    prob = self._selection_prob(good)
                    selected = np.random.random(good.size) < prob
                    active[good[~selected]] = 1

        if finalize:
            self._finalize_blacklist()
        return refresh

    def _finalize_blacklist(self):
        '''Finalize local blacklist-derived counters after any distribution.'''
        cdef int d
        cdef uint8_t[:] active_view

        active_view = <uint8_t[:self.iter.tot_num_data]>self.iter.blacklist
        self.iter.num_blacklist = 0
        for d in range(self.iter.tot_num_data):
            if self.iter.blacklist[d] != 0:
                self.iter.num_blacklist += 1
            elif self.iter.scale != NULL and self.iter.scale[d] < 0.:
                self.iter.scale[d] = 1.
        if self.iter.num_blacklist == self.iter.tot_num_data:
            raise ValueError('All frames are blacklisted')
        self._update_last_selected(active_view)

    def normalize_scale(self):
        '''Normalize scale factors.'''
        cdef long x, d
        cdef double mean_scale

        blist = self.blacklist
        if blist is None:
            mean_scale = self.scale.mean()
        else:
            mean_scale = self.scale[blist==0].mean()

        if self.iter.mod != NULL and self.iter.mod.model1 != NULL:
            for x in range(self.iter.mod.vol):
                self.iter.mod.model1[x] *= mean_scale

        if blist is None:
            for d in range(self.iter.tot_num_data):
                self.iter.scale[d] /= mean_scale
        elif self.iter.par.data_fraction < 1.:
            for d in range(self.iter.tot_num_data):
                if self.iter.static_blacklist[d] == 0 and self.iter.scale[d] > 0.:
                    self.iter.scale[d] /= mean_scale
        else:
            for d in range(self.iter.tot_num_data):
                if self.iter.blacklist[d] == 0:
                    self.iter.scale[d] /= mean_scale

        self.iter.rms_change *= mean_scale

    def calc_beta(self, setval=None):
        '''Calculate beta values for each frame.

        Args:
            setval (float): Set all betas to this value. Default None.
        '''
        if self.iter.dset == NULL:
            raise AttributeError('Add data first before calculating beta for each frame')

        cdef double start

        if setval is not None:
            start = setval
        else:
            start = -1.

        c_iterate.calc_beta(start, self.iter)

    def free(self):
        '''Free allocated iteration memory.'''
        if self.iter == NULL:
            return

        if self.iter.fcounts != NULL: free(self.iter.fcounts)
        if self.iter.scale != NULL: free(self.iter.scale)
        if self.iter.bgscale != NULL: free(self.iter.bgscale)
        if self.iter.beta != NULL: free(self.iter.beta)
        if self.iter.beta_start != NULL: free(self.iter.beta_start)
        if self.iter.sum_fact != NULL: free(self.iter.sum_fact)
        if self.iter.blacklist != NULL: free(self.iter.blacklist)
        if self.iter.static_blacklist != NULL: free(self.iter.static_blacklist)
        if self.iter.last_selected != NULL: free(self.iter.last_selected)
        if self.iter.quat_mapping != NULL: free(self.iter.quat_mapping)
        if self.iter.num_rel_quat != NULL: free(self.iter.num_rel_quat)
        if self.iter.det_mapping != NULL: free(self.iter.det_mapping)
        if self.iter.rescale != NULL: free(self.iter.rescale)
        if self.iter.mean_count != NULL: free(self.iter.mean_count)

        if self.iter.rel_quat != NULL:
            for d in range(iter.tot_num_data):
                free(self.iter.rel_quat[d])
            free(self.iter.rel_quat)
        if self.iter.rel_prob != NULL:
            for d in range(iter.tot_num_data):
                free(self.iter.rel_prob[d])
            free(self.iter.rel_prob)
        self.iter = NULL

    def _gen_detlist(self):
        '''Generate list of unique detectors.'''
        cdef int i, d, ind
        cdef dataset *curr = self.iter.dset

        fnames = []
        self.iter.num_dfiles = 0
        while curr != NULL:
            fnames.append(curr.det.fname)
            self.iter.num_dfiles += 1
            curr = curr.next
        _, indices, mapping = np.unique(fnames, return_index=True, return_inverse=True)
        self.iter.num_det = len(indices)
        print(self.iter.num_det, 'unique detectors')

        self.iter.det = <detector*> calloc(len(indices), sizeof(detector))
        for i, ind in enumerate(indices):
            curr = self.iter.dset
            for d in range(ind):
                curr = curr.next
            memcpy(&self.iter.det[i], curr.det, sizeof(detector))

        self.iter.det_mapping = <int*> malloc(self.iter.num_dfiles * sizeof(int))
        for i, d in enumerate(mapping):
            self.iter.det_mapping[i] = d

    @staticmethod
    def calculate_size(qmax):
        '''Calculate model size from qmax.

        Args:
            qmax (float): Maximum q value.

        Returns:
            int: Model grid size.
        '''
        return int(2 * np.ceil(qmax) + 3)

    @property
    def tot_num_data(self):
        '''Total number of frames.'''
        return self.iter.tot_num_data
    @tot_num_data.setter
    def tot_num_data(self, val):
        '''Set total number of frames.'''
        self.iter.tot_num_data = val
    @property
    def fcounts(self):
        '''Frame photon counts.'''
        return np.asarray(<int[:self.tot_num_data]>self.iter.fcounts) if self.iter.fcounts != NULL else None
    @property
    def scale(self):
        '''Per-frame scale factors.'''
        return np.asarray(<double[:self.tot_num_data]>self.iter.scale) if self.iter.scale != NULL else None
    @scale.setter
    def scale(self, arr):
        '''Set scale factors.'''
        if len(arr.shape) != 1 or arr.dtype != 'f8':
            raise ValueError('scale must be 1D array of float64 dtype')
        if arr.shape[0] != self.iter.tot_num_data:
            raise ValueError('tot_num_data mismatch. One scale factor per frame (%d vs %d)'%(arr.shape[0], self.det.num_pix))

        if self.iter.scale == NULL:
            self.iter.scale = <double*> malloc(arr.size * sizeof(double))
        for i in range(arr.size):
            self.iter.scale[i] = arr[i]
    @property
    def bgscale(self):
        '''Per-frame background scale factors.'''
        return np.asarray(<double[:self.tot_num_data]>self.iter.bgscale) if self.iter.bgscale != NULL else None
    @property
    def beta(self):
        '''Per-frame beta values.'''
        return np.asarray(<double[:self.tot_num_data]>self.iter.beta) if self.iter.beta != NULL else None
    @property
    def beta_start(self):
        '''Initial per-frame beta values.'''
        return np.asarray(<double[:self.tot_num_data]>self.iter.beta_start) if self.iter.beta_start != NULL else None
    @property
    def blacklist(self):
        '''Frame blacklist mask (1=blacklisted, 0=good).'''
        return np.asarray(<uint8_t[:self.tot_num_data]>self.iter.blacklist) if self.iter.blacklist != NULL else None
    @property
    def static_blacklist(self):
        '''Static frame blacklist mask before stochastic sampling.'''
        return np.asarray(<uint8_t[:self.tot_num_data]>self.iter.static_blacklist) if self.iter.static_blacklist != NULL else None
    @property
    def last_selected(self):
        '''Last iteration where each frame was selected for stochastic EMC.'''
        return np.asarray(<int[:self.tot_num_data]>self.iter.last_selected) if self.iter.last_selected != NULL else None
    @property
    def det_mapping(self):
        '''Mapping from frames to detectors.'''
        return np.asarray(<int[:self.num_dfiles]>self.iter.det_mapping) if self.iter.det_mapping != NULL else None
    @property
    def mean_count(self):
        '''Mean photon count per detector.'''
        return np.asarray(<double[:self.num_det]>self.iter.mean_count) if self.iter.mean_count != NULL else None
    @property
    def rescale(self):
        '''Per-detector rescale factors.'''
        return np.asarray(<double[:self.num_det]>self.iter.rescale) if self.iter.rescale != NULL else None
    @property
    def likelihood(self):
        '''Current log-likelihood value.'''
        return self.iter.likelihood
    @property
    def mutual_info(self):
        '''Current mutual information value.'''
        return self.iter.mutual_info
    @property
    def rms_change(self):
        '''RMS change between iterations.'''
        return self.iter.rms_change

    @property
    def num_det(self):
        '''Number of unique detectors.'''
        return self.iter.num_det
    @property
    def num_dfiles(self):
        '''Number of data files.'''
        return self.iter.num_dfiles
    @property
    def dets(self):
        '''List of detector objects.'''
        retval = [None for d in range(self.iter.num_det)]
        for d in range(self.iter.num_det):
            curr = CDetector()
            curr.free()
            curr.det = &self.iter.det[d]
            retval[d] = curr
        return retval
    @property
    def model(self):
        '''Current model.'''
        if self.iter.mod == NULL:
            return
        retval = Model()
        retval.free()
        retval.mod = self.iter.mod
        return retval
    @property
    def quat(self):
        '''Quaternion orientations.'''
        if self.iter.quat == NULL:
            return
        retval = Quaternion()
        retval.free()
        retval.quat = self.iter.quat
        return retval
    @property
    def data(self):
        '''Dataset.'''
        if self.iter.dset == NULL:
            return
        retval = CDataset()
        retval.free()
        retval.dset = self.iter.dset
        return retval
    @property
    def params(self):
        '''Parameters.'''
        if self.iter.par == NULL:
            return
        retval = EMCParams()
        retval.free()
        retval.par = self.iter.par
        return retval
