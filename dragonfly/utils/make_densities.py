#!/usr/bin/env python

'''Make electron density volume from PDB file and configuration parameters.

This module generates 3D electron density maps from PDB structure files
using atomic scattering factors. The density is computed on a cubic grid
and optionally filtered to remove high-frequency artifacts.

Functions:
    make_dens: Generate density map from configuration file.
    main: Command-line interface.
'''

import sys
import argparse
import logging

from .py_src import read_config
from .py_src import process_pdb
from .py_src import py_utils

def make_dens(config_fname, yes=False, verbose=False):
    '''Generate density map from parameters in config file.

    Args:
        config_fname (str): Path to configuration file.
        yes (bool): Skip confirmation prompts. Default False.
        verbose (bool): Enable verbose logging. Default False.

    Returns:
        None: Outputs binary density file specified in config.

    Raises:
        :py:exc:`configparser.NoOptionError`: If required config options are missing.
    '''
    config = read_config.MyConfigParser()
    config.read(config_fname)
    log_file = config.get_filename('make_densities', 'log_file', fallback='simdata.log')

    logging.basicConfig(filename=log_file, level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s')
    logging.info('\n\nStarting.... make_densities')
    logging.info(' '.join(sys.argv))

    pdb_code = None
    pdb_file = config.get_filename('make_densities', 'in_pdb_file', fallback=None)
    if pdb_file is None:
        pdb_code = config.get_filename('make_densities', 'pdb_code')
        pdb_file = 'aux/%s.pdb' % pdb_code.upper()

    num_threads = config.getint('make_densities', 'num_threads', fallback=4)

    aux_dir = config.get_filename('make_densities', 'scatt_dir')
    dens_fname = config.get_filename('make_densities', 'out_density_file')

    if not (yes or py_utils.check_to_overwrite(dens_fname)):
        return

    timer = py_utils.MyTimer()
    pm = config.get_detector_config(show=verbose)
    q_pm = read_config.compute_q_params(pm['detd'], pm['dets_x'],
                                        pm['dets_y'], pm['pixsize'],
                                        pm['wavelength'], pm['ewald_rad'], show=verbose)
    timer.reset('Reading experiment parameters', report=verbose)

    if pdb_code is not None:
        process_pdb.fetch_pdb(pdb_code)
    all_atoms = process_pdb.process(pdb_file, aux_dir, pm['wavelength'])
    timer.reset('Reading PDB', report=verbose)

    den = process_pdb.atoms_to_density_map(all_atoms, q_pm['half_p_res'])
    lp_den = process_pdb.low_pass_filter_density_map(den, threads=num_threads)
    timer.reset('Creating density map', report=verbose)

    py_utils.write_density(dens_fname, lp_den, binary=True)
    timer.reset('Writing densities to file', report=verbose)

    timer.report_time_since_beginning()

def main():
    '''Parse command line arguments and generate electron density volume with config file'''
    parser = argparse.ArgumentParser(description='Make density map from PDB')
    parser.add_argument('-c', '--config_fname', default='config.ini',
                        help='Path to config file (Default: config.ini)')
    parser.add_argument('-v', '--verbose', action='store_true', help='Verbose logging')
    parser.add_argument('-y', '--yes', action='store_true',
                        help='Say yes to all prompts')
    args = parser.parse_args()


    make_dens(args.config_fname, yes=args.yes, verbose=args.verbose)

if __name__ == '__main__':
    main()
