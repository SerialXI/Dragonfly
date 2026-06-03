import sys
import os.path as op
import time
import argparse
import logging
import configparser

import numpy as np
import dragonfly
from cython.parallel import parallel, prange

from .py_src import read_config, py_utils
from . cimport make_data as c_make_data
from dragonfly.model cimport Model
from dragonfly.detector cimport CDetector

def rand_quat():
    while True:
        qvals = np.random.random(4) - 0.5
        qnorm = np.linalg.norm(qvals)
        if qnorm < 0.25:
            break
    return qvals / qnorm

cdef make_data_cdef(config_fname, yes=False, verbose=False):
    config = read_config.MyConfigParser()
    config.read(config_fname)

    num_data = config.getint('make_data', 'num_data')
    fluence = config.getfloat('make_data', 'fluence', fallback=-1.)
    mean_count = config.getfloat('make_data', 'mean_count', fallback=-1.)
    do_gamma = config.getboolean('make_data', 'do_gamma', fallback=False)
    intens_fname = config.get_filename('make_data', 'in_intensity_file')
    det_fname = config.get_filename('make_data', 'in_detector_file')
    try:
        background_fname = config.get_filename('make_data', 'background_file')
    except configparser.NoOptionError:
        background_fname = None
    out_fname = config.get_filename('make_data', 'out_photons_file').encode('UTF-8')
    likelihood_fname = config.get_filename('make_data', 'out_likelihood_file', fallback='').encode('UTF-8')
    scale_fname = config.get_filename('make_data', 'out_scale_file', fallback='').encode('UTF-8')

    if not (yes or py_utils.check_to_overwrite(out_fname.decode())):
        return

    timer = py_utils.MyTimer()
    model = Model(1)
    try:
        model.allocate(intens_fname)
    except ValueError as err:
        size = err.args[1]
        model.free()
        model = Model(size)
        model.allocate(intens_fname)

    det = CDetector(det_fname, norm=False)
    
    if background_fname:
        det.parse_background(background_fname)

    hdf5_output = True
    if op.splitext(out_fname)[1].decode('utf-8') == '.emc':
        hdf5_output = False
    timer.reset('Ready to generate data', report=verbose)

    mean_count = c_make_data.rescale_intens(fluence, mean_count, model.mod, det.det)
    if mean_count < 0:
        return
    timer.reset('Rescaled model', report=verbose)

    cdef char* c_out_fname = out_fname
    cdef char* c_likelihood_fname = likelihood_fname
    cdef char* c_scale_fname = scale_fname
    mean_count = c_make_data.gen_and_save_dataset(num_data, mean_count, int(do_gamma),
                                                  c_out_fname, int(hdf5_output),
                                                  c_likelihood_fname, c_scale_fname,
                                                  model.mod, det.det)
    timer.reset('Generated and saved data', report=verbose)
    timer.report_time_since_beginning()

    print('Generated %d frames with %f photons/frame' % (num_data, mean_count))
    print('Time taken = %f s' % (time.time() - timer._time_start))

def make_data(config_fname, yes=False, verbose=False):
    make_data_cdef(config_fname, yes=yes, verbose=verbose)

def main():
    '''Parse command line arguments and generate electron density volume with config file'''
    parser = argparse.ArgumentParser(description='Make data file from 3D intensities')
    parser.add_argument('-c', '--config_fname', default='config.ini',
                        help='Path to config file (Default: config.ini)')
    parser.add_argument('-v', '--verbose', action='store_true', help='Verbose logging')
    parser.add_argument('-y', '--yes', action='store_true',
                        help='Say yes to all prompts')
    args = parser.parse_args()

    logging.basicConfig(filename='simdata.log', level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s')
    logging.info('\n\nStarting.... make_data')
    logging.info(' '.join(sys.argv))

    make_data_cdef(args.config_fname, yes=args.yes, verbose=args.verbose)

if __name__ == '__main__':
    main()
