#!/usr/bin/env python

'''Module to generate detector file from configuration parameters.

This module creates detector geometry files containing pixel coordinates,
correction factors, and masks. The detector is computed based on
experimental geometry parameters in the configuration file.

Functions:
    make_detector: Generate detector file from configuration.
    main: Command-line interface.
'''

import sys
import logging
import argparse

import numpy as np
import dragonfly

from .py_src import read_config
from .py_src import py_utils

def make_detector(config_fname, yes=False, verbose=False):
    '''Generate detector file from parameters in config file.

    Args:
        config_fname (str): Path to configuration file.
        yes (bool): Skip confirmation prompts. Default False.
        verbose (bool): Enable verbose logging. Default False.

    Returns:
        None: Outputs detector file specified in config.
    '''
    config = read_config.MyConfigParser()
    config.read(config_fname)

    log_file = config.get_filename('make_detector', 'log_file', fallback='simdata.log')

    logging.basicConfig(filename=log_file, level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s')
    logging.info('\n\nStarting make_detector....')
    logging.info(' '.join(sys.argv))



    det_fname = config.get_filename('make_detector', 'out_detector_file')

    if not (yes or py_utils.check_to_overwrite(det_fname)):
        return

    timer = py_utils.MyTimer()
    pm = config.get_detector_config(show=verbose)
    read_config.compute_q_params(pm['detd'], pm['dets_x'],
                                 pm['dets_y'], pm['pixsize'],
                                 pm['wavelength'], pm['ewald_rad'], show=verbose)
    timer.reset('Reading experiment parameters', report=verbose)

    det = dragonfly.Detector()
    x, y = np.meshgrid(np.arange(pm['dets_x'])-pm['detc_x'],
                       np.arange(pm['dets_y'])-pm['detc_y'],
                       indexing='ij')
    det.cx = x.flatten()
    det.cy = y.flatten()
    det.detd = pm['detd'] / pm['pixsize']
    det.ewald_rad = pm['ewald_rad']
    det.calc_from_coords(pol=pm['polarization'])

    rad = np.sqrt(det.cx**2 + det.cy**2)
    if pm['mask_fname'] is None:
        det.raw_mask = np.zeros(det.corr.shape, dtype='u1')
        det.raw_mask[rad > min(pm['detc_x'], pm['detc_y'])] = 1
        det.raw_mask[rad < pm['stoprad']] = 2
    else:
        det.raw_mask = np.fromfile(pm['mask_fname'], '=u1')
        det.raw_mask[(rad > min(pm['detc_x'], pm['detc_y'])) & (det.raw_mask == 0)] = 1
    timer.reset('Creating detector', report=verbose)

    det.write(det_fname)
    timer.reset('Writing detector to %s' % det_fname, report=verbose)

    timer.report_time_since_beginning()

def main():
    '''Parses command line arguments and config file to generate detector file'''
    parser = argparse.ArgumentParser(description='Make detector file')
    parser.add_argument('-c', '--config_fname', default='config.ini',
                        help='Path to config file (Default: config.ini)')
    parser.add_argument('-v', '--verbose', action='store_true', help='Verbose logging')
    parser.add_argument('-y', '--yes', action='store_true',
                        help='Say yes to all prompts')
    args = parser.parse_args()

    make_detector(args.config_fname, yes=args.yes, verbose=args.verbose)

if __name__ == '__main__':
    main()
