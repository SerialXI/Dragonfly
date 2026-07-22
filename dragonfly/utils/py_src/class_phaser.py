import sys
import numpy as np
from scipy import ndimage

class ClassPhaser():
    def __init__(self, intens, num_supp=None, maxfrac=None,
                 positivity=True):
        self.fobs = np.empty_like(intens)
        self.fobs[intens>=0] = np.sqrt(intens[intens>=0])
        self.fobs[intens<0] = -1
        self.rel_qpix = (self.fobs >= 0)
        self.positivity = positivity

        if maxfrac is None and num_supp is None:
            raise ValueError('Need either num_supp or maxfrac for shrinkwrap')
        elif maxfrac is None:
            self.num_supp = num_supp
            self.maxfrac = None
        elif num_supp is None:
            self.maxfrac = maxfrac
            self.num_supp = None
        else:
            raise ValueError('Cannot use both num_supp and maxfrac')

    def phase(self, algorithms, num_runs=1, progress_callback=None):
        if num_runs < 1:
            raise ValueError('Need at least one phasing run')

        if num_runs == 1:
            self._phase_once(algorithms, progress_callback=progress_callback)
            return

        num_valid = 0
        total = None
        support_total = None
        for run in range(num_runs):
            self._phase_once(algorithms, progress_callback=progress_callback,
                             run_num=run+1, num_runs=num_runs)
            if np.isnan(self.current).sum() > 0:
                sys.stderr.write('NaNs in density, ignoring run %d\n' % (run+1))
                continue
            if self.current.sum() == 0.:
                sys.stderr.write('Zero-valued density, ignoring run %d\n' % (run+1))
                continue

            current = self.current
            support = self.support.astype('f8')
            if total is not None:
                current, support = self._align(total, current, support)

            if total is None:
                total = current.copy()
                support_total = support.copy()
            else:
                total += current
                support_total += support
            num_valid += 1

        if num_valid == 0:
            raise ValueError('No valid phasing runs')

        self.current = total / num_valid
        self.support = support_total / num_valid
        self.num_valid_runs = num_valid

    def _phase_once(self, algorithms, progress_callback=None, run_num=None, num_runs=None):
        self.current = np.random.random(self.fobs.shape)
        for i, algo in enumerate(algorithms):
            self.run_iteration(algo)
            status = '%d/%d: %s' % (i+1, len(algorithms), algo)
            if run_num is not None and num_runs is not None:
                status = 'Run %d/%d, %s' % (run_num, num_runs, status)
            if progress_callback is not None:
                progress_callback(status)
            else:
                sys.stderr.write('\rIteration %s' % status)
        if progress_callback is None:
            sys.stderr.write('\n')

        # Shift support center to center of array
        x, y = np.indices(self.current.shape)
        cen = self.current.shape[-1] // 2
        shift = int(cen - (x*self.support).sum() / self.support.sum()), int(cen - (y*self.support).sum() / self.support.sum())
        self.current = np.roll(self.current, shift, axis=(0,1))
        self.support = np.roll(self.support, shift, axis=(0,1))
        shift = int(cen - (x*self.support).sum() / self.support.sum()), int(cen - (y*self.support).sum() / self.support.sum())
        self.current = np.roll(self.current, shift, axis=(0,1))
        self.support = np.roll(self.support, shift, axis=(0,1))

    def _align(self, total, current, support):
        fcurrent = np.fft.fftn(np.abs(current))
        ftotal = np.fft.fftn(np.abs(total))
        corr = np.abs(np.fft.ifftn(ftotal * np.conj(fcurrent)))
        icorr = np.abs(np.fft.ifftn(ftotal * fcurrent))

        if corr.max() > icorr.max():
            pos = np.unravel_index(corr.argmax(), corr.shape)
            return (np.roll(current, pos, axis=(0, 1)),
                    np.roll(support, pos, axis=(0, 1)))

        pos = tuple([1+x for x in np.unravel_index(icorr.argmax(), icorr.shape)])
        return (np.roll(current[::-1, ::-1], pos, axis=(0, 1)),
                np.roll(support[::-1, ::-1], pos, axis=(0, 1)))
    
    def run_iteration(self, algo='DM'):
        if algo == 'ER':
            self.current = self.proj_fourier(self.proj_direct(self.current))
        elif algo == 'DM':
            p1 = self.proj_fourier(self.current)
            p2 = self.proj_direct(2 * p1 - self.current)
            self.current += p2 - p1
        else:
            raise ValueError('Unknown algorithm name: %s' % algo)

    def proj_fourier(self, dens):
        fdens = np.fft.fftshift(np.fft.fftn(dens))
        sel = self.rel_qpix & (np.abs(fdens) > 0)
        fdens[sel] *= self.fobs[sel] / np.abs(fdens[sel])
        return np.real(np.fft.ifftn(np.fft.ifftshift(fdens)))

    def proj_direct(self, dens):
        out_dens = np.copy(dens)

        # Shrinkwrap / Volume constraint
        smdens = ndimage.gaussian_filter(dens, 2)
        if self.num_supp is not None:
            thresh = np.sort(smdens.ravel())[-self.num_supp]
        else:
            thresh = smdens.max() * self.maxfrac
        self.support = smdens > thresh
        out_dens[~self.support] = 0

        if self.positivity:
            out_dens[out_dens < 0] = 0

        return out_dens

def preproc(intens, rad, minrad, maxrad, filtsize=30):
    out = intens.copy()
    out[rad<=minrad] = 1e20
    out[rad>=maxrad] = 1e20
    out -= ndimage.minimum_filter(out, filtsize)
    out[rad<=minrad] = -1
    #out[rad<=minrad] = 0
    out[rad>=maxrad] = 0
    return out
