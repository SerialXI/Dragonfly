import os.path as op
import configparser

from libc.stdlib cimport malloc, calloc, free
from libc.string cimport strcpy
from . cimport params as c_params
from .params cimport EMCParams, params

cdef class EMCParams:
    '''Class storing reconstruction parameters.

    Holds all configuration parameters for EMC reconstruction including
    iteration settings, algorithm parameters, and symmetry options.
    '''

    def __init__(self):
        '''Initialize EMCParams with default values.'''
        self.par = <c_params.params*> calloc(1, sizeof(c_params.params))

        self.par.rank = 0
        self.par.num_proc = 1
        self.par.start_iter = 1
        self.par.known_scale = 0

        self.par.log_fname = <char*> malloc(1024)
        strcpy(self.par.log_fname, b'EMC.log')
        self.par.output_folder= <char*> malloc(1024)
        strcpy(self.par.output_folder , b'data/')
        self.par.rtype = c_params.RECON3D

        self.par.need_scaling = 0
        self.par.alpha = 0.
        self.par.data_fraction = 1.
        self.par.beta_factor = -1.
        self.par.radius = 0.
        self.par.num_modes = 1
        self.par.nonrot_modes = 0
        self.par.rot_per_mode = 0

        self.par.beta_jump = 1.
        self.par.beta_period = 100
        self.par.data_fraction_period = 1
        self.par.radius_jump = 0.
        self.par.radius_period = 100

        self.par.oversampling = 10.
        self.par.sigmasq = 0.
        self.par.friedel_sym = 0
        self.par.axial_sym = 1
        self.par.save_prob = 0
        self.par.update_scale = 1
        self.par.refine = 0
        self.par.coarse_div = 0
        self.par.fine_div = 0
        self.par.fixed_seed = 0

    def free(self):
        '''Free allocated parameter memory.'''
        if self.par == NULL:
            return

        if self.par.log_fname != NULL:
            free(self.par.log_fname)
        if self.par.output_folder != NULL:
            free(self.par.output_folder)
        free(self.par)
        self.par = NULL

    def from_config(self, config_fname, section_name='emc'):
        '''Load parameters from configuration file.

        Args:
            config_fname (str): Path to configuration file.
            section_name (str): Section name. Default 'emc'.
        '''
        from .utils.py_src.read_config import MyConfigParser
        config = MyConfigParser()
        config.read(config_fname)

        try:
            output_folder = config.get_filename(section_name, 'output_folder', fallback='data/')
        except configparser.NoOptionError:
            output_folder = 'data/'
        strcpy(self.par.output_folder, output_folder.encode())

        try:
            log_fname = config.get_filename(section_name, 'log_file', fallback='EMC.log')
        except configparser.NoOptionError:
            log_fname = 'EMC.log'
        strcpy(self.par.log_fname, log_fname.encode())
        rtype = config.get(section_name, 'recon_type', fallback='3d')
        if rtype == '3d':
            self.par.rtype = c_params.RECON3D
        elif rtype == '2d':
            self.par.rtype = c_params.RECON2D
        elif rtype == 'rz':
            self.par.rtype = c_params.RECONRZ
        else:
            raise ValueError('Unknown recon_type %s'%rtype)
        
        self.par.verbosity = config.getint(section_name, 'verbosity', fallback=9)
        self.par.need_scaling = config.getint(section_name, 'need_scaling', fallback=0)
        self.par.alpha = config.getfloat(section_name, 'alpha', fallback=0.)
        self.par.data_fraction = config.getfloat(section_name, 'data_fraction', fallback=1.)
        if self.par.data_fraction <= 0. or self.par.data_fraction > 1.:
            raise ValueError('data_fraction must be > 0 and <= 1')
        self.par.data_fraction_period = config.getint(section_name, 'data_fraction_period', fallback=1)
        if self.par.data_fraction_period < 1:
            raise ValueError('data_fraction_period must be >= 1')
        self.par.beta_factor = config.getfloat(section_name, 'beta_factor', fallback=1.)
        self.par.radius = config.getfloat(section_name, 'radius', fallback=0.)
        self.par.num_modes = config.getint(section_name, 'num_modes', fallback=1)
        self.par.nonrot_modes = config.getint(section_name, 'num_nonrot_modes', fallback=0)
        self.par.rot_per_mode = config.getint(section_name, 'num_rot', fallback=0)

        beta_schedule = tuple(config.get(section_name, 'beta_schedule', fallback='1. 100').split())
        self.par.beta_jump = float(beta_schedule[0])
        self.par.beta_period = int(beta_schedule[1])
        radius_schedule= tuple(config.get(section_name, 'radius_schedule', fallback='0. 100').split())
        self.par.radius_jump = float(radius_schedule[0])
        self.par.radius_period = int(radius_schedule[1])

        self.par.oversampling = config.getfloat(section_name, 'oversampling', fallback=10.)
        self.par.sigmasq = config.getfloat(section_name, 'gaussian_sigma', fallback=0.)**2
        self.par.friedel_sym = config.getint(section_name, 'friedel_sym', fallback=0)
        self.par.axial_sym = config.getint(section_name, 'axial_sym', fallback=1)
        self.par.save_prob = config.getint(section_name, 'save_prob', fallback=0)
        self.par.update_scale = config.getint(section_name, 'update_scale', fallback=1)
        num_divs = config.get(section_name, 'num_div', fallback='0').split()
        if len(num_divs) > 1:
            self.par.refine = 1
            self.par.fine_div = int(num_divs[0])
            self.par.coarse_div = int(num_divs[1])
            print('Doing refinement from num_div = %d -> %d\n' % (self.par.coarse_div, self.par.fine_div))
        
        self.par.fixed_seed = config.getint(section_name, 'fixed_seed', fallback=0)

    @property
    def rank(self):
        '''MPI rank.'''
        return self.par.rank
    @rank.setter
    def rank(self, int val):
        '''Set MPI rank.'''
        self.par.rank = val
    @property
    def num_proc(self):
        '''Number of MPI processes.'''
        return self.par.num_proc
    @num_proc.setter
    def num_proc(self, int val):
        '''Set number of MPI processes.'''
        self.par.num_proc = val
    @property
    def iteration(self):
        '''Current iteration number.'''
        return self.par.iteration
    @iteration.setter
    def iteration(self, int val):
        '''Set current iteration number.'''
        self.par.iteration = val
    @property
    def start_iter(self):
        '''Starting iteration number.'''
        return self.par.start_iter
    @start_iter.setter
    def start_iter(self, int val):
        '''Set starting iteration number.'''
        self.par.start_iter = val
    @property
    def num_iter(self):
        '''Number of iterations.'''
        return self.par.num_iter
    @num_iter.setter
    def num_iter(self, int val):
        '''Set number of iterations.'''
        self.par.num_iter = val
    @property
    def beta_period(self):
        '''Iterations per beta period.'''
        return self.par.beta_period
    @beta_period.setter
    def beta_period(self, int val):
        '''Set iterations per beta period.'''
        self.par.beta_period = val
    @property
    def data_fraction_period(self):
        '''Iterations to keep each stochastic data subset.'''
        return self.par.data_fraction_period
    @data_fraction_period.setter
    def data_fraction_period(self, int val):
        '''Set stochastic data subset period.'''
        if val < 1:
            raise ValueError('data_fraction_period must be >= 1')
        self.par.data_fraction_period = val
    @property
    def need_scaling(self):
        '''Whether frame scaling is enabled.'''
        return self.par.need_scaling
    @need_scaling.setter
    def need_scaling(self, int val):
        '''Set need_scaling flag.'''
        self.par.need_scaling = val
    @property
    def known_scale(self):
        '''Whether scale factors are known.'''
        return self.par.known_scale
    @known_scale.setter
    def known_scale(self, int val):
        '''Set known_scale flag.'''
        self.par.known_scale = val
    @property
    def update_scale(self):
        '''Whether to update scale factors during iteration.'''
        return self.par.update_scale
    @update_scale.setter
    def update_scale(self, int val):
        '''Set update_scale flag.'''
        self.par.update_scale = val
    @property
    def save_prob(self):
        '''Whether to save probabilities.'''
        return self.par.save_prob
    @save_prob.setter
    def save_prob(self, int val):
        '''Set save_prob flag.'''
        self.par.save_prob = val
    @property
    def verbosity(self):
        '''Verbosity level.'''
        return self.par.verbosity
    @verbosity.setter
    def verbosity(self, int val):
        '''Set verbosity level.'''
        self.par.verbosity = val
    @property
    def friedel_sym(self):
        '''Whether Friedel symmetry is applied.'''
        return self.par.friedel_sym
    @friedel_sym.setter
    def friedel_sym(self, int val):
        '''Set friedel_sym flag.'''
        self.par.friedel_sym = val
    @property
    def axial_sym(self):
        '''N-fold axial symmetry order.'''
        return self.par.axial_sym
    @axial_sym.setter
    def axial_sym(self, int val):
        '''Set axial symmetry order.'''
        self.par.axial_sym = val
    @property
    def refine(self):
        '''Whether refinement mode is enabled.'''
        return self.par.refine
    @refine.setter
    def refine(self, int val):
        '''Set refine flag.'''
        self.par.refine = val
    @property
    def coarse_div(self):
        '''Coarse division for refinement.'''
        return self.par.coarse_div
    @coarse_div.setter
    def coarse_div(self, int val):
        '''Set coarse division.'''
        self.par.coarse_div = val
    @property
    def fine_div(self):
        '''Fine division for refinement.'''
        return self.par.fine_div
    @fine_div.setter
    def fine_div(self, int val):
        '''Set fine division.'''
        self.par.fine_div = val
    @property
    def radius_period(self):
        '''Iterations per radius period.'''
        return self.par.radius_period
    @radius_period.setter
    def radius_period(self, int val):
        '''Set radius period.'''
        self.par.radius_period = val
    @property
    def num_modes(self):
        '''Number of modes.'''
        return self.par.num_modes
    @num_modes.setter
    def num_modes(self, int val):
        '''Set number of modes.'''
        self.par.num_modes = val
    @property
    def rot_per_mode(self):
        '''Rotations per mode.'''
        return self.par.rot_per_mode
    @rot_per_mode.setter
    def rot_per_mode(self, int val):
        '''Set rotations per mode.'''
        self.par.rot_per_mode = val
    @property
    def nonrot_modes(self):
        '''Number of non-rotating modes.'''
        return self.par.nonrot_modes
    @nonrot_modes.setter
    def nonrot_modes(self, int val):
        '''Set non-rotating modes.'''
        self.par.nonrot_modes = val

    @property
    def rtype(self):
        '''Returns the reconstruction type (RECON3D, RECON2D, or RECONRZ).'''
        return ['RECON3D', 'RECON2D', 'RECONRZ'][self.par.rtype]
    @rtype.setter
    def rtype(self, val): 
        '''Set reconstruction type.'''
        if val == 'RECON3D':
            self.par.rtype = RECON3D
        elif val == 'RECON2D':
            self.par.rtype = RECON2D
        elif val == 'RECONRZ':
            self.par.rtype = RECONRZ
        else:
            raise ValueError('Unknown recon_type %s'%val)
            
    @property
    def alpha(self):
        '''Model smoothing factor.'''
        return self.par.alpha
    @alpha.setter
    def alpha(self, double val):
        '''Set alpha smoothing factor.'''
        self.par.alpha = val
    @property
    def data_fraction(self):
        '''Fraction of non-blacklisted frames to use per iteration.'''
        return self.par.data_fraction
    @data_fraction.setter
    def data_fraction(self, double val):
        '''Set stochastic EMC data fraction.'''
        if val <= 0. or val > 1.:
            raise ValueError('data_fraction must be > 0 and <= 1')
        self.par.data_fraction = val
    @property
    def beta_jump(self):
        '''Beta value jump per period.'''
        return self.par.beta_jump
    @beta_jump.setter
    def beta_jump(self, double val):
        '''Set beta jump.'''
        self.par.beta_jump = val
    @property
    def beta_factor(self):
        '''Beta scaling factor.'''
        return self.par.beta_factor
    @beta_factor.setter
    def beta_factor(self, double val):
        '''Set beta factor.'''
        self.par.beta_factor = val
    @property
    def radius(self):
        '''Reconstruction radius.'''
        return self.par.radius
    @radius.setter
    def radius(self, double val):
        '''Set radius.'''
        self.par.radius = val
    @property
    def radius_jump(self):
        '''Radius jump per period.'''
        return self.par.radius_jump
    @radius_jump.setter
    def radius_jump(self, double val):
        '''Set radius jump.'''
        self.par.radius_jump = val
    @property
    def oversampling(self):
        '''Pixel oversampling factor.'''
        return self.par.oversampling
    @oversampling.setter
    def oversampling(self, double val):
        '''Set oversampling factor.'''
        self.par.oversampling = val
    @property
    def sigmasq(self):
        '''Squared Gaussian sigma.'''
        return self.par.sigmasq
    @sigmasq.setter
    def sigmasq(self, double val):
        '''Set sigmasq.'''
        self.par.sigmasq = val
    
    @property
    def log_fname(self):
        '''Path to log file.'''
        return (<bytes> self.par.log_fname).decode() if self.par.log_fname != NULL else None
    @log_fname.setter
    def log_fname(self, val):
        '''Set log file path.'''
        if self.par.log_fname != NULL:
            free(self.par.log_fname)
        self.par.log_fname = <char*> malloc(len(val) + 1)
        strcpy(self.par.log_fname, bytes(val, 'utf-8'))
    @property
    def output_folder(self):
        '''Path to output folder.'''
        return (<bytes> self.par.output_folder).decode() if self.par.output_folder != NULL else None
    @output_folder.setter
    def output_folder(self, val):
        '''Set output folder path.'''
        if self.par.output_folder != NULL:
            free(self.par.output_folder)
        self.par.output_folder = <char*> malloc(len(val) + 1)
        strcpy(self.par.output_folder, bytes(val, 'utf-8'))
    @property
    def fixed_seed(self):
        '''Fixed seed for testing.'''
        return self.par.fixed_seed
    @fixed_seed.setter
    def fixed_seed(self, int val):
        '''Set fixed seed flag.'''
        self.par.fixed_seed = val
