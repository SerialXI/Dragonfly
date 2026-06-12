'''Module to generate and manipulate uniformly sampled 3D quaternions'''

from __future__ import print_function
import sys
import os.path as op
from socket import gethostname
import itertools
from configparser import ConfigParser

import numpy as np
import h5py
from scipy.spatial import distance

NUM_VERT = 120
NUM_EDGE = 720
NUM_FACE = 1200
NUM_CELL = 600
NNN = 12

cimport numpy as np
from libc.stdlib cimport malloc, calloc, free
from . cimport quaternion as c_quat
from .quaternion cimport Quaternion


def _parse_pair(config, section_name, option, default, cast):
    '''Parse one or two values from a config option.'''
    vals = config.get(section_name, option, fallback=default).split()
    if len(vals) == 1:
        return cast(vals[0]), cast(vals[0])
    if len(vals) == 2:
        return cast(vals[0]), cast(vals[1])
    raise ValueError('%s must contain one or two values' % option)

cdef class Quaternion:
    '''Class for generating and manipulating orientations for SO(3).

    Supports 3D quaternion generation with various point group symmetries,
    as well as 2D in-plane rotation generation.

    Args:
        num_div (int): Number of divisions for 3D sampling. Default 0.
        num_rot (int): Number of rotations for 2D sampling. Default 0.
        point_group (str): Point group symmetry ('1', '2', 'S4', 'A5'). Default ''.

    Example:
        >>> quat = Quaternion()
        >>> quat.generate_3d(10)
        >>> quat = Quaternion(num_rot=36, point_group='1')
    '''

    def __init__(self, int num_div=0, int num_rot=0, point_group=''):
        '''Initialize Quaternion object.'''
        self.quat = <c_quat.quaternion*> malloc(sizeof(c_quat.quaternion))
        self.quat.octahedral_flag = (point_group == 'S4')
        self.quat.icosahedral_flag = (point_group == 'A5')
        self.quat.num_rot_p = 0
        self.quat.quats = NULL
        self.reduced = False

        if num_div > 0 and num_rot > 0:
            raise ValueError('Cannot specify both num_div and num_rot')
        elif num_div > 0:
            self.generate_3d(num_div)
        elif num_rot > 0:
            self.generate_2d(num_rot, point_group)

    def from_config(self, config_fname, section_name='emc'):
        '''Load orientations from configuration file.

        Args:
            config_fname (str): Path to configuration file.
            section_name (str): Section name. Default 'emc'.
        '''
        config_folder = op.dirname(config_fname)
        config = ConfigParser()
        config.read(config_fname)

        rtype = config.get(section_name, 'recon_type', fallback='3d').lower()

        if rtype == '2d':
            num_rot = config.getint(section_name, 'num_rot')
            friedel_sym = config.getboolean(section_name, 'friedel_sym', fallback=False)
            point_group = '2' if friedel_sym else '1'
            shift_max = _parse_pair(config, section_name, 'shift_max', '0.', float)
            num_shifts = _parse_pair(config, section_name, 'num_shifts', '1', int)
            self.generate_2d(num_rot, point_group, shift_max=shift_max, num_shifts=num_shifts)
        elif rtype =='3d':
            num_div = config.get(section_name, 'num_div', fallback='0').split()
            self.generate_3d(int(num_div[0]))
            if len(num_div) == 2:
                self.voronoi_subset(int(num_div[1]))

            point_group = config.get(section_name, 'point_group', fallback='1')
            if point_group == 'A5':
                print('Reducing to A5 sub-group')
                self.reduce_icosahedral()
                self.quat.icosahedral_flag = True
            elif point_group == 'S4':
                print('Reducing to S4 sub-group')
                self.reduce_octahedral()
                self.quat.octahedral_flag = True
            elif point_group != '1':
                raise ValueError('Unknown point group: ' + point_group)

    def generate_3d(self, int num_div):
        '''Generate quaternions for SO(3) with uniform sampling.

        Args:
            num_div (int): Number of rotational subdivisions.
        '''
        self.quat.num_div = num_div
        self.reduced = False
        c_quat.quat_gen(num_div, self.quat)

    def generate_2d(self, int num_rot, point_group, shift_max=(0., 0.), num_shifts=(1, 1)):
        '''Generate quaternions for in-plane rotations.

        Args:
            num_rot (int): Number of rotational samples.
            point_group (str): N-fold rotational symmetry ('1', '2', etc.).
            shift_max (tuple): Maximum absolute x/y beam-center shifts.
            num_shifts (tuple): Number of x/y beam-center shifts.
        '''
        if point_group == '':
            point_group = '1'
        max_angle = 2 * np.pi / int(point_group)

        if len(shift_max) != 2 or len(num_shifts) != 2:
            raise ValueError('shift_max and num_shifts must have two expanded values')
        shift_max_x, shift_max_y = float(shift_max[0]), float(shift_max[1])
        num_shift_x, num_shift_y = int(num_shifts[0]), int(num_shifts[1])
        if shift_max_x < 0. or shift_max_y < 0.:
            raise ValueError('shift_max must be >= 0')
        if num_shift_x < 1 or num_shift_y < 1:
            raise ValueError('num_shifts must be >= 1')
        if (num_shift_x > 1 and shift_max_x == 0.) or (num_shift_y > 1 and shift_max_y == 0.):
            raise ValueError('shift_max must be > 0 for dimensions with num_shifts > 1')

        shifts_x = np.array([0.]) if num_shift_x == 1 else np.linspace(-shift_max_x, shift_max_x, num_shift_x)
        shifts_y = np.array([0.]) if num_shift_y == 1 else np.linspace(-shift_max_y, shift_max_y, num_shift_y)
        angles = np.arange(num_rot) * max_angle / num_rot
        angle_grid, shift_x_grid, shift_y_grid = np.meshgrid(
            angles, shifts_x, shifts_y, indexing='ij')
        q = np.zeros((angle_grid.size, 5))
        q[:,0] = angle_grid.ravel()
        q[:,1] = shift_x_grid.ravel()
        q[:,2] = shift_y_grid.ravel()
        q[:,4] = 1. / q.shape[0]
        self.quats = q

    def save(self, fname):
        '''Save quaternions to HDF5 file.

        Args:
            fname (str): Output file path.
        '''
        if self.quat.quats == NULL:
            raise AttributeError('Generate quaternion first before saving')

        with h5py.File(fname, 'w') as fptr:
            fptr['quaternions'] = self.quats

    def parse(self, fname):
        '''Load quaternions from HDF5 file.

        Args:
            fname (str): Input file path.
        '''
        cdef int i

        with h5py.File(fname, 'r') as fptr:
            quats = fptr['quaternions'][:]

        if quats.shape[1] == 4:
            quats = np.pad(quats, ((0,0), (0,1)), constant_values=1./quats.shape[0])
        elif quats.shape[1] == 5:
            quats[:,4] /= quats[:,4].sum()
        else:
            raise ValueError('Need 4- or 5-dimensional quaternions')
        self.quat.num_rot = quats.shape[0]

        self.quat.quats = <double*> malloc(self.num_rot * 5 * sizeof(double))
        quats = quats.ravel()
        for i in np.arange(5*self.num_rot):
            self.quat.quats[i] = quats[i]

    def divide(self, rank, num_proc, num_modes=1, num_nonrot_modes=0):
        '''Divide rotations among MPI processes.

        Args:
            rank (int): MPI rank.
            num_proc (int): Total number of MPI processes.
            num_modes (int): Number of modes. Default 1.
            num_nonrot_modes (int): Non-rotating modes. Default 0.

        Returns:
            int: Number of rotations for this process.
        '''
        tot_num_rot = num_modes * self.num_rot + num_nonrot_modes
        self.num_rot_p = tot_num_rot / num_proc
        if rank < (tot_num_rot % num_proc):
            self.quat.num_rot_p += 1
        if num_proc > 1:
            sys.stderr.write('%d: %s: num_rot_p = %d/%d\n' % (rank, gethostname(), self.num_rot_p, tot_num_rot))
            sys.stderr.flush()
        return self.quat.num_rot_p

    def reduce_icosahedral(self, return_sym=False):
        '''Reduce to icosahedral point group subset.

        Args:
            return_sym (bool): Return symmetry elements. Default False.

        Returns:
            int or :py:class:`numpy.ndarray`: Number of rotations or symmetry quaternions.
        '''
        if self.quat.quats == NULL:
            raise AttributeError('Generate quaternion first before reducing')
        if self.reduced:
            print('Already reduced')
            return self.num_rot

        cdef int r, t

        quats = self.quats
        vert_dist = distance.cdist(quats[:60,:4], quats[60:,:4])
        sel = np.where(vert_dist.argmin(axis=0) == 8)[0]
        self.num_rot = 1 + sel.shape[0]

        for t in range(5):
            self.quat.quats[t] = self.quat.quats[8*5 + t]
        for r in range(1, self.num_rot):
            for t in range(5):
                self.quat.quats[r*5 + t] = self.quat.quats[(sel[r-1]+60)*5 + t]

        self.reduced = True
        if return_sym:
            return quats[:60,:4]
        else:
            return self.num_rot

    def reduce_octahedral(self, return_sym=False):
        '''Reduce to octahedral point group subset.

        Args:
            return_sym (bool): Return symmetry elements. Default False.

        Returns:
            int or :py:class:`numpy.ndarray`: Number of rotations or symmetry quaternions.
        '''
        if self.quat.quats == NULL:
            raise AttributeError('Generate quaternion first before reducing')
        if self.reduced:
            print('Already reduced')
            return self.num_rot

        cube_quat = np.zeros((24, 4))
        cube_quat[:4] = np.identity(4)

        signs3 = np.stack(np.meshgrid([-1., 1.], [-1., 1.], [-1., 1.], indexing='ij'), -1).reshape(-1, 3)
        cube_quat[4:12] = 0.5
        cube_quat[4:12,1:] *= signs3

        val = np.sqrt(0.5)
        cube_quat[12:18, 0] = val
        cube_quat[12:15, 1:] = np.identity(3) * val
        cube_quat[15:18, 1:] = np.identity(3) * val * -1

        for i, pos in enumerate(itertools.permutations(range(3), 2)):
            cube_quat[18 + i, np.array(pos) + 1] = [val, -val] if pos[0] < pos[1] else [val, val]

        cdef int r, t
        quats = self.quats

        cube_dist = distance.cdist(cube_quat, quats[:,:4])
        sel = np.where(cube_dist.argmin(axis=0) == 0)[0]
        self.num_rot = sel.shape[0]
        for r in range(self.num_rot):
            for t in range(5):
                self.quat.quats[r*5 + t] = self.quat.quats[sel[r]*5 + t]

        self.reduced = True
        if return_sym:
            return cube_quat
        else:
            return self.num_rot

    def voronoi_subset(self, int coarse_num_div):
        '''Select subset using Voronoi cells.

        Args:
            coarse_num_div (int): Coarse division level.

        Returns:
            :py:class:`numpy.ndarray`: Nearest coarse division indices.
        '''
        if self.num_rot == 0:
            print('Generate fine quaternions first')
            return
        cdef c_quat.quaternion* coarse
        coarse = <c_quat.quaternion*> calloc(1, sizeof(c_quat.quaternion))
        coarse.num_div = coarse_num_div
        c_quat.quat_gen(coarse_num_div, coarse)

        nearest = np.empty(self.num_rot, dtype='i4')
        cdef int[:] nearest_view = nearest
        c_quat.voronoi_subset(coarse, self.quat, &nearest_view[0])

        free(coarse.quats)
        free(coarse)

        return nearest

    def free(self):
        '''Free allocated quaternion memory.'''
        if self.quat == NULL:
            return
        if self.quat.quats != NULL: free(self.quat.quats)
        free(self.quat)
        self.quat = NULL

    @property
    def num_div(self):
        '''Number of rotational subdivisions.'''
        return self.quat.num_div
    @property
    def num_rot(self):
        '''Total number of rotations.'''
        return self.quat.num_rot
    @num_rot.setter
    def num_rot(self, val):
        '''Set total number of rotations.'''
        self.quat.num_rot = val
    @property
    def num_rot_p(self):
        '''Number of rotations for current MPI process.'''
        return self.quat.num_rot_p
    @num_rot_p.setter
    def num_rot_p(self, val):
        '''Set rotations for current MPI process.'''
        self.quat.num_rot_p = val
    @property
    def reduced(self):
        '''Whether quaternions have been reduced to point group.'''
        return bool(self.quat.reduced)
    @reduced.setter
    def reduced(self, val):
        '''Set reduced flag.'''
        self.quat.reduced = int(val)
    @property
    def icosahedral_flag(self):
        '''Whether using icosahedral symmetry.'''
        return bool(self.quat.icosahedral_flag)
    @property
    def octahedral_flag(self):
        '''Whether using octahedral symmetry.'''
        return bool(self.quat.octahedral_flag)
    @property
    def quats(self):
        '''Quaternion array, shape (num_rot, 5).'''
        return np.asarray(<double[:self.num_rot*5]> self.quat.quats).reshape(-1,5) if self.quat.quats != NULL else None
    @quats.setter
    def quats(self, arr):
        '''Set quaternion array.'''
        if len(arr.shape) != 2 or arr.shape[1] != 5:
            raise ValueError('quats must be  2D array of shape (N, 5)')
        if arr.dtype != 'f8':
            raise TypeError('quats must be double precision floats (float64)')

        self.quat.quats = <double*> malloc(arr.size * sizeof(double))
        for i in range(arr.size):
            self.quat.quats[i] = arr.ravel()[i]
        self.quat.num_rot = arr.shape[0]
