'''Miscellaneous utility functions called by one or more python modules'''

from __future__ import print_function
from builtins import input
import time
import argparse
import os
import sys
import logging
from six.moves import configparser
import numpy as np
from . import detector, reademc, readvol

class MyTimer(object):
    '''Class to report elapsed time for logging'''
    def __init__(self):
        self._time0 = time.time()
        self._time_start = self._time0

    def reset(self):
        '''Update time'''
        self._time0 = time.time()

    def reset_global(self):
        '''Update global start time'''
        self._time_start = time.time()

    def reset_and_report(self, msg):
        '''Update time and log difference since last reset'''
        time1 = time.time()
        logging.info("%-30s:%5.5f seconds", msg, time1-self._time0)
        self._time0 = time1

    def report_time_since_beginning(self):
        '''Log time difference since start'''
        logging.info("="*20)
        logging.info("%-30s:%5.5f seconds", "Since beginning", time.time() - self._time_start)

class MyArgparser(argparse.ArgumentParser):
    '''Modified ArgumentParser to add default arguments common across all utilities'''
    def __init__(self, description=""):
        argparse.ArgumentParser.__init__(self, description=description,
                                         formatter_class=argparse.RawTextHelpFormatter)
        self.add_argument("-c", "--config_file", dest="config_file",
                          help="config file (default config.ini)")
        self.add_argument("-v", "--verbose", dest="vb", action="store_true", default=False)
        self.add_argument("-m", "--main_dir", dest="main_dir",
                          help="relative path to main repository directory\n"
                               "(where data aux utils are stored)")
        self.add_argument("-y", "--yes", help="say yes to all question prompts",
                          action="store_true")

    def special_parse_args(self):
        '''Parse command line arguments
        Log if config file and main directory used are defaults
        '''
        args = self.parse_args()
        if not args.config_file:
            args.config_file = "config.ini"
            logging.info("Config file not specified. Using %s", args.config_file)
        if not args.main_dir:
            args.main_dir = os.path.split(os.path.abspath(args.config_file))[0]
            logging.info("Main directory not specified. Using %s", args.main_dir)
        return args

def write_density(in_den_file, in_den, binary=True):
    '''Write density volume to file (binary or text format)'''
    if binary:
        in_den.astype('float64').tofile(in_den_file)
    else:
        with open(in_den_file, "w") as fptr:
            for line0 in in_den:
                for line1 in line0:
                    tmp = ' '.join(line1.astype('str'))
                    fptr.write(tmp + '\n')

def read_density(in_den_file):
    '''Read density volume from file (binary)'''
    den = np.fromfile(in_den_file, dtype="float64")
    vol = len(den)
    size = int(np.round(np.power(vol, 1./3.)))
    out_den = den.reshape(3*(size,))
    return out_den

def check_to_overwrite(fname):
    '''Check if file exists and prompt before overwriting
    By default, the file is overwritten
    '''
    overwrite = True
    yes_val = set(['yes', 'y', '', 'yup', 'ya'])
    no_val = set(['no', 'n', 'nope', 'nay', 'not'])
    if os.path.isfile(fname):
        sys.stdout.write("%s is present. Overwrite? [Y or Return/N]: " % fname)
        choice = input().lower()
        if choice in yes_val:
            overwrite = True
            print("Overwriting " + fname)
            logging.info("Overwriting %s", fname)
        elif choice in no_val:
            overwrite = False
            print("Not overwriting " + fname)
            logging.info("Not overwriting %s", fname)
        else:
            sys.stdout.write("Please respond with 'yes' or 'no'")
            overwrite = False
    return overwrite

def confirm_oversampling(ratio):
    '''Print message for user if oversampling ratio is too high
    Prompt for continuation
    '''
    proceed = True
    done = False
    yes_val = set(['yes', 'y', '', 'yup', 'ya'])
    no_val = set(['no', 'n', 'nope', 'nay', 'not'])
    print('Oversampling ratio = %.2f is a little high. This is inefficient.' % ratio)
    print('Please see https://dragonfly-spi.readthedocs.io/en/latest/background/oversampling.html for tips')
    sys.stdout.write('Continue anyway? [Y or Return/N]: ')
    while not done:
        choice = input().lower()
        if choice in yes_val:
            proceed = True
            done = True
            logging.info("Continuing with minimum oversampling ratio = %.2f", ratio)
        elif choice in no_val:
            proceed = False
            done = True
            logging.info("Minimum oversampling ratio of %.2f too high. Quitting", ratio)
        else:
            sys.stdout.write("Please respond with 'yes' or 'no': ")
            proceed = False
    return proceed

def increment_quat_file_sensibly(config_fname, incr):
    '''Increments num_div in config file by incr'''
    config = configparser.ConfigParser()
    config.read(config_fname)

    quat_num_div = int(config.get("emc", "num_div"))
    logging.info("Setting quaternion from %s to %s", str(quat_num_div), str(quat_num_div+incr))
    quat_num_div += incr
    config.set("emc", "num_div", quat_num_div)

    with open(config_fname, "w") as fptr:
        config.write(fptr)

def gen_det_and_emc(gui, classifier=False, mask=False):
    '''Creates EMCReader and Detector instances for GUIs'''
    if len(set(gui.det_list)) == 1:
        geom_list = [detector.Detector(gui.det_list[0], gui.detd, gui.ewald_rad, mask_flag=mask)]
        geom_mapping = None
    else:
        if classifier:
            print('The Classifier GUI will likely have problems with multiple geometries')
            print('We recommend classifying patterns with a common geometry')
        uniq = sorted(set(gui.det_list))
        geom_list = [detector.Detector(fname, gui.detd, gui.ewald_rad, mask_flag=mask)
                     for fname in uniq]
        geom_mapping = [uniq.index(fname) for fname in gui.det_list]
    gui.geom = geom_list[0]
    gui.emc_reader = reademc.EMCReader(gui.photons_list, geom_list, geom_mapping)

class DummyGeom(object):
    pass

def gen_stack(gui):
    size = gui.stack_size
    gui.emc_reader = readvol.VolReader(gui.photons_list[0], size)
    gui.geom = DummyGeom()
    ix, iy = np.indices((size,size))
    gui.geom.cx = ix.ravel().astype('f8') - size//2
    gui.geom.cy = iy.ravel().astype('f8') - size//2
    gui.geom.unassembled_mask = np.ones(size*size)

def increment_beta_sensibly(config_fname, incr):
    config = ConfigParser.ConfigParser()
    config.read(config_fname)

    beta = float(config.get("emc", "beta"))
    msg = ["Setting beta from", str(beta), "to", str(beta*incr)]
    logging.info(' '.join(msg))
    beta *= incr
    config.set("emc", "beta", beta)
    
    with open("config.ini", "w") as fptr:
        config.write(fptr)
