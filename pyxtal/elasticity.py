# ======================================================================
# matscipy - Python materials science tools
# https://github.com/libAtoms/matscipy
#
# Copyright (2014) James Kermode, King's College London
#                  Lars Pastewka, Karlsruhe Institute of Technology
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# ======================================================================

from __future__ import print_function

import itertools
import warnings

import numpy as np
from numpy.linalg import inv, norm
try:
    import scipy.stats as scipy_stats
except:
    scipy_stats = None

import ase.units as units
from ase.atoms import Atoms
from ase.io import read
from time import time
###

# The indices of the full stiffness matrix of (orthorhombic) interest
Voigt_notation = [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)]

def full_3x3_to_Voigt_6_index(i, j):
    if i == j:
        return i
    return 6-i-j

###

def Voigt_6_to_full_3x3_strain(strain_vector):
    """
    Form a 3x3 strain matrix from a 6 component vector in Voigt notation
    """
    e1, e2, e3, e4, e5, e6 = np.transpose(strain_vector)
    return np.transpose([[1.0+e1, 0.5*e6, 0.5*e5],
                         [0.5*e6, 1.0+e2, 0.5*e4],
                         [0.5*e5, 0.5*e4, 1.0+e3]])


def Voigt_6_to_full_3x3_stress(stress_vector):
    """
    Form a 3x3 stress matrix from a 6 component vector in Voigt notation
    """
    s1, s2, s3, s4, s5, s6 = np.transpose(stress_vector)
    return np.transpose([[s1, s6, s5],
                         [s6, s2, s4],
                         [s5, s4, s3]])


def full_3x3_to_Voigt_6_strain(strain_matrix):
    """
    Form a 6 component strain vector in Voigt notation from a 3x3 matrix
    """
    strain_matrix = np.asarray(strain_matrix)
    return np.transpose([strain_matrix[...,0,0] - 1.0,
                         strain_matrix[...,1,1] - 1.0,
                         strain_matrix[...,2,2] - 1.0,
                         strain_matrix[...,1,2]+strain_matrix[...,2,1],
                         strain_matrix[...,0,2]+strain_matrix[...,2,0],
                         strain_matrix[...,0,1]+strain_matrix[...,1,0]])


def full_3x3_to_Voigt_6_stress(stress_matrix):
    """
    Form a 6 component stress vector in Voigt notation from a 3x3 matrix
    """
    stress_matrix = np.asarray(stress_matrix)
    return np.transpose([stress_matrix[...,0,0],
                         stress_matrix[...,1,1],
                         stress_matrix[...,2,2],
                         (stress_matrix[...,1,2]+stress_matrix[...,1,2])/2,
                         (stress_matrix[...,0,2]+stress_matrix[...,0,2])/2,
                         (stress_matrix[...,0,1]+stress_matrix[...,0,1])/2])


def Voigt_6x6_to_full_3x3x3x3(C):
    """
    Convert from the Voigt representation of the stiffness matrix to the full
    3x3x3x3 representation.

    Parameters
    ----------
    C : array_like
        6x6 stiffness matrix (Voigt notation).
    
    Returns
    -------
    C : array_like
        3x3x3x3 stiffness matrix.
    """
    
    C = np.asarray(C)
    C_out = np.zeros((3,3,3,3), dtype=float)
    for i, j, k, l in itertools.product(range(3), range(3), range(3), range(3)):
        Voigt_i = full_3x3_to_Voigt_6_index(i, j)
        Voigt_j = full_3x3_to_Voigt_6_index(k, l)
        C_out[i, j, k, l] = C[Voigt_i, Voigt_j]
    return C_out


def full_3x3x3x3_to_Voigt_6x6(C):
    """
    Convert from the full 3x3x3x3 representation of the stiffness matrix
    to the representation in Voigt notation. Checks symmetry in that process.
    """

    tol = 1e-3

    C = np.asarray(C)
    Voigt = np.zeros((6,6))
    for i in range(6):
        for j in range(6):
            k, l = Voigt_notation[i]
            m, n = Voigt_notation[j]
            Voigt[i,j] = C[k,l,m,n]

            #print '---'
            #print k,l,m,n, C[k,l,m,n]
            #print m,n,k,l, C[m,n,k,l]
            #print l,k,m,n, C[l,k,m,n]
            #print k,l,n,m, C[k,l,n,m]
            #print m,n,l,k, C[m,n,l,k]
            #print n,m,k,l, C[n,m,k,l]
            #print l,k,n,m, C[l,k,n,m]
            #print n,m,l,k, C[n,m,l,k]
            #print '---'

            # Check symmetries
            assert abs(Voigt[i,j]-C[m,n,k,l]) < tol, \
                'Voigt[{},{}] = {}, C[{},{},{},{}] = {}' \
                .format(i, j, Voigt[i,j], m, n, k, l, C[m,n,k,l])
            assert abs(Voigt[i,j]-C[l,k,m,n]) < tol, \
                'Voigt[{},{}] = {}, C[{},{},{},{}] = {}' \
                .format(i, j, Voigt[i,j], k, l, m, n, C[l,k,m,n])
            assert abs(Voigt[i,j]-C[k,l,n,m]) < tol, \
                'Voigt[{},{}] = {}, C[{},{},{},{}] = {}' \
                .format(i, j, Voigt[i,j], k, l, n, m, C[k,l,n,m])
            assert abs(Voigt[i,j]-C[m,n,l,k]) < tol, \
                'Voigt[{},{}] = {}, C[{},{},{},{}] = {}' \
                .format(i, j, Voigt[i,j], m, n, l, k, C[m,n,l,k])
            assert abs(Voigt[i,j]-C[n,m,k,l]) < tol, \
                'Voigt[{},{}] = {}, C[{},{},{},{}] = {}' \
                .format(i, j, Voigt[i,j], n, m, k, l, C[n,m,k,l])
            assert abs(Voigt[i,j]-C[l,k,n,m]) < tol, \
                'Voigt[{},{}] = {}, C[{},{},{},{}] = {}' \
                .format(i, j, Voigt[i,j], l, k, n, m, C[l,k,n,m])
            assert abs(Voigt[i,j]-C[n,m,l,k]) < tol, \
                'Voigt[{},{}] = {}, C[{},{},{},{}] = {}' \
                .format(i, j, Voigt[i,j], n, m, l, k, C[n,m,l,k])

    return Voigt


def Voigt_6x6_to_cubic(C):
    """
    Convert the Voigt 6x6 representation into the cubic elastic constants
    C11, C12 and C44.
    """

    tol = 1e-6

    C_check = np.zeros_like(C)
    C_check[np.diag_indices_from(C_check)] = C[np.diag_indices_from(C)]
    C_check[0:3,0:3] = C[0:3,0:3]
    if np.any(np.abs(C-C_check) > tol):
        raise ValueError('"C" does not have cubic symmetry.')

    C11s = np.array([C[0,0], C[1,1], C[2,2]])
    C12s = np.array([C[1,2], C[0,2], C[0,1]])
    C44s = np.array([C[3,3], C[4,4], C[5,5]])

    C11 = np.mean(C11s)
    C12 = np.mean(C12s)
    C44 = np.mean(C44s)

    if np.any(np.abs(C11-C11s) > tol) or np.any(np.abs(C12-C12s) > tol) or \
            np.any(np.abs(C44-C44s) > tol):
        raise ValueError('"C" does not have cubic symmetry.')

    return np.array([C11, C12, C44])


def cubic_to_Voigt_6x6(C11, C12, C44):
    return np.array([[C11,C12,C12,  0,  0,  0],
                     [C12,C11,C12,  0,  0,  0],
                     [C12,C12,C11,  0,  0,  0],
                     [  0,  0,  0,C44,  0,  0],
                     [  0,  0,  0,  0,C44,  0],
                     [  0,  0,  0,  0,  0,C44]])

def _invariants(s, syy=None, szz=None, syz=None, sxz=None, sxy=None,
               full_3x3_to_Voigt_6=full_3x3_to_Voigt_6_stress):
    """
    Receives a list of stress tensors and returns the three tensor invariants.
    """
    if syy is None:
        s = np.asarray(s)
        if s.shape == (6,):
            s = s.reshape(1,-1)
        elif s.shape == (3,3):
            s = full_3x3_to_Voigt_6(s)
        elif s.shape[-1] == 3 and s.shape[-2] == 3:
            s = full_3x3_to_Voigt_6(s)
    else:
        s = np.transpose([np.transpose(s),
                          np.transpose(syy),
                          np.transpose(szz),
                          np.transpose(syz),
                          np.transpose(sxz),
                          np.transpose(sxy)])
    I1 = s[...,0]+s[...,1]+s[...,2]
    I2 = s[...,0]*s[...,1]+s[...,1]*s[...,2]+s[...,2]*s[...,0]-s[...,3]**2- \
        s[...,4]**2-s[...,5]**2
    I3 = s[...,0]*s[...,1]*s[...,2]+2*s[...,3]*s[...,4]*s[...,5]- \
        s[...,3]**2*s[...,2]-s[...,4]**2*s[...,0]-s[...,5]**2*s[...,1]

    return I1, I2, I3

def invariants(s, syy=None, szz=None, syz=None, sxz=None, sxy=None,
               full_3x3_to_Voigt_6=full_3x3_to_Voigt_6_stress):
    I1, I2, I3 = _invariants(s, syy=syy, szz=szz, syz=syz, sxz=sxz, sxy=sxy,
                             full_3x3_to_Voigt_6=full_3x3_to_Voigt_6)

    J2 = I1**2/3-I2
    J3 = 2*I1**3/27-I1*I2/3+I3

    # Return hydrostatic pressure, octahedral shear stress and J3
    return -I1/3, np.sqrt(2*J2/3), J3

###

def rotate_cubic_elastic_constants(C11, C12, C44, A, tol=1e-6):
    """
    Return rotated elastic moduli for a cubic crystal given the elastic 
    constant in standard C11, C12, C44 notation.

    Parameters
    ----------
    C11, C12, C44 : float
        Cubic elastic moduli.
    A : array_like
        3x3 rotation matrix.

    Returns
    -------
    C : array
        6x6 matrix of rotated elastic constants (Voigt notation).
    """

    A = np.asarray(A)

    # Is this a rotation matrix?
    if np.sometrue(np.abs(np.dot(np.array(A), np.transpose(np.array(A))) - 
                          np.eye(3, dtype=float)) > tol):
        raise RuntimeError('Matrix *A* does not describe a rotation.')

    # Invariant elastic constants
    la = C12
    mu = C44
    al = C11 - la - 2*mu

    # Construct rotated C in Voigt notation
    C = [ ]
    for i, j in Voigt_notation:
        for k, l in Voigt_notation:
            h = 0.0
            if i == j and k == l:
                h += la
            if i == k and j == l:
                h += mu
            if i == l and j == k:
                h += mu
            h += al*np.sum(A[i,:]*A[j,:]*A[k,:]*A[l,:])
            C += [ h ]

    C = np.array(C)
    C.shape = (6, 6)
    return C

###

def rotate_elastic_constants(C, A, tol=1e-6):
    """
    Return rotated elastic moduli for a general crystal given the elastic 
    constant in Voigt notation.

    Parameters
    ----------
    C : array_like
        6x6 matrix of elastic constants (Voigt notation).
    A : array_like
        3x3 rotation matrix.

    Returns
    -------
    C : array
        6x6 matrix of rotated elastic constants (Voigt notation).
    """

    A = np.asarray(A)

    # Is this a rotation matrix?
    if np.sometrue(np.abs(np.dot(np.array(A), np.transpose(np.array(A))) - 
                          np.eye(3, dtype=float)) > tol):
        raise RuntimeError('Matrix *A* does not describe a rotation.')

    # Rotate
    return full_3x3x3x3_to_Voigt_6x6(np.einsum('ia,jb,kc,ld,abcd->ijkl',
                                               A, A, A, A,
                                               Voigt_6x6_to_full_3x3x3x3(C)))

###

class CubicElasticModuli:
    tol = 1e-6

    def __init__(self, C11, C12, C44):
        """
        Initialize a cubic system with elastic constants C11, C12, C44
        """

        warnings.warn('CubicElasticModuli is deprecated. Use '
                      'rotate_elastic_constants function instead.')

        # la, mu, al are the three invariant elastic constants
        self.la = C12
        self.mu = C44
        self.al = C11 - self.la - 2*self.mu

        A = np.eye(3, dtype=float)

        # Compute initial stiffness matrix
        self.rotate(A)


    def rotate(self, A):
        """
        Compute the rotated stiffness matrix
        """

        A = np.asarray(A)

        # Is this a rotation matrix?
        if np.sometrue(np.abs(np.dot(np.array(A), np.transpose(np.array(A))) - 
                              np.eye(3, dtype=float)) > self.tol):
            raise RuntimeError('Matrix *A* does not describe a rotation.')

        C = [ ]
        for i, j in Voigt_notation:
            for k, l in Voigt_notation:
                h = 0.0
                if i == j and k == l:
                    h += self.la
                if i == k and j == l:
                    h += self.mu
                if i == l and j == k:
                    h += self.mu
                h += self.al*np.sum(A[i,:]*A[j,:]*A[k,:]*A[l,:])
                C += [ h ]

        self.C = np.asarray(C)
        self.C.shape = (6, 6)
        return self.C


    def _rotate_explicit(self, A):
        """
        Compute the rotated stiffness matrix by applying the rotation to the
        full stiffness matrix. This function is for debugging purposes only.
        """

        A = np.asarray(A)

        # Is this a rotation matrix?
        if np.sometrue(np.abs(np.dot(np.array(A), np.transpose(np.array(A))) - 
                              np.eye(3, dtype=float) ) > self.tol):
            raise RuntimeError('Matrix *A* does not describe a rotation.')

        C = np.zeros((3, 3, 3, 3), dtype=float)

        # Construct unrotated stiffness matrix
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for m in range(3):
                        h = 0.0
                        if i == j and k == m:
                            h += self.la
                        if i == k and j == m:
                            h += self.mu
                        if i == m and j == k:
                            h += self.mu
                        if i == j and j == k and k == m:
                            h += self.al
                        C[i, j, k, m] = h

        # Rotate
        C = np.einsum('ia,jb,kc,ld,abcd->ijkl', A, A, A, A, C)

        self.C = full_3x3x3x3_to_Voigt_6x6(C)
        return self.C


    def stiffness(self):
        """
        Return the elastic constants
        """

        return self.C


    def compliance(self):
        """
        Return the compliance coefficients
        """

        return inv(self.C)

###

def measure_triclinic_elastic_constants(a, delta=0.001, optimizer=None, 
                                        logfile=None, **kwargs):
    """
    Brute-force measurement of elastic constants for a triclinic (general)
    unit cell.

    Parameters
    ----------
    a : ase.Atoms
        Atomic configuration.
    optimizer : ase.optimizer.*
        Optimizer to use for atomic position. Does not optimize atomic
        position if set to None.
    delta : float
        Strain increment for analytical derivatives of stresses.

    Returns
    -------
    C : array_like
        6x6 matrix of the elastic constants in Voigt notation.
    """

    if optimizer is not None:
        optimizer(a, logfile=logfile).run(**kwargs)

    r0 = a.positions.copy()

    cell = a.cell.copy()
    volume = a.get_volume()

    C = np.zeros((3,3,3,3), dtype=float)
    for i in range(3):
        for j in range(3):
            a.set_cell(cell, scale_atoms=True)
            a.set_positions(r0)
        
            D = np.eye(3)
            D[i, j] += 0.5*delta
            D[j, i] += 0.5*delta
            a.set_cell(np.dot(D, cell.T).T, scale_atoms=True)
            if optimizer is not None:
                optimizer(a, logfile=logfile).run(**kwargs)
            sp = Voigt_6_to_full_3x3_stress(a.get_stress()*a.get_volume())

            D = np.eye(3)
            D[i, j] -= 0.5*delta
            D[j, i] -= 0.5*delta
            a.set_cell(np.dot(D, cell.T).T, scale_atoms=True)
            if optimizer is not None:
                optimizer(a, logfile=logfile).run(**kwargs)
            sm = Voigt_6_to_full_3x3_stress(a.get_stress()*a.get_volume())

            C[:,:,i,j] = (sp-sm)/(2*delta*volume)

    a.set_cell(cell, scale_atoms=True)
    a.set_positions(r0)

    return full_3x3x3x3_to_Voigt_6x6(C)



Cij_symmetry = {
   'cubic':           np.array([[1, 7, 7, 0, 0, 0],
                                [7, 1, 7, 0, 0, 0],
                                [7, 7, 1, 0, 0, 0],
                                [0, 0, 0, 4, 0, 0],
                                [0, 0, 0, 0, 4, 0],
                                [0, 0, 0, 0, 0, 4]]),

   'trigonal_high':   np.array([[1, 7, 8, 9, 10, 0],
                                [7, 1, 8, 0,-9, 0],
                                [8, 8, 3, 0, 0, 0],
                                [9, -9, 0, 4, 0, 0],
                                [10, 0, 0, 0, 4, 0],
                                [0, 0, 0, 0, 0, 6]]),

   'trigonal_low':    np.array([[1,  7,  8,  9,  10,  0 ],
                                [7,  1,  8, -9, -10,  0 ],
                                [8,  8,  3,  0,   0,  0 ],
                                [9, -9,  0,  4,   0, -10],
                                [10,-10, 0,  0,   4,  9 ],
                                [0,  0,  0, -10 , 9,  6 ]]),

   'tetragonal_high': np.array([[1, 7, 8, 0, 0, 0],
                                [7, 1, 8, 0, 0, 0],
                                [8, 8, 3, 0, 0, 0],
                                [0, 0, 0, 4, 0, 0],
                                [0, 0, 0, 0, 4, 0],
                                [0, 0, 0, 0, 0, 6]]),

   'tetragonal_low':  np.array([[1, 7, 8, 0, 0, 11],
                                [7, 1, 8, 0, 0, -11],
                                [8, 8, 3, 0, 0, 0],
                                [0, 0, 0, 4, 0, 0],
                                [0, 0, 0, 0, 4, 0],
                                [11, -11, 0, 0, 0, 6]]),

   'orthorhombic':    np.array([[ 1,  7,  8,  0,  0,  0],
                                [ 7,  2, 12,  0,  0,  0],
                                [ 8, 12,  3,  0,  0,  0],
                                [ 0,  0,  0,  4,  0,  0],
                                [ 0,  0,  0,  0,  5,  0],
                                [ 0,  0,  0,  0,  0,  6]]),

   'monoclinic':      np.array([[ 1,  7,  8,  0,  10,  0],
                                [ 7,  2, 12,  0, 14,  0],
                                [ 8, 12,  3,  0, 17,  0],
                                [ 0,  0,  0,  4,  0,  20],
                                [10, 14, 17,  0,  5,  0],
                                [ 0,  0,  0, 20,  0,  6]]),
    
    'triclinic':       np.array([[ 1,  7,  8,  9,  10, 11],
                                 [ 7,  2, 12,  13, 14, 15],
                                 [ 8, 12,  3,  16, 17, 18],
                                 [ 9, 13, 16,  4,  19, 20],
                                 [10, 14, 17, 19,  5,  21],
                                 [11, 15, 18, 20,  21, 6 ]]),
   }

def _dec(pattern):
    # Decrease C_ij indices patterns by one.
    # Used to make strain_patterns definitions below more readable.
    return [(i-1, j-1) for (i,j) in pattern]

strain_patterns = {

   'cubic': [
      # strain pattern e1+e4, yields C11, C21, C31 and C44, then C12 is average of C21 and C31
      [ np.array([1,0,0,1,0,0]), _dec([(1,1), (2,1), (3,1), (4,4)])]
   ],

   'trigonal_high': [
      # strain pattern e3 yield C13, C23 and C33
      [ np.array([0,0,1,0,0,0]), _dec([(1,3), (2,3), (3,3)])],

      # strain pattern e1+e4 yields C11 C21 C31 and C44
      [ np.array([1,0,0,1,0,0]), _dec([(1,1), (2,1), (3,1), (4,4)])],

      # strain pattern e1 yields C11 C21 C31 C41 C51
      [ np.array([1,0,0,0,0,0]), _dec([(1,1), (2,1), (3,1), (4,1), (5,1)])],

      # strain pattern e3+e4
      [ np.array([0,0,1,1,0,0]), _dec([(3,3), (4,4)])]

   ],

   'trigonal_low': [
     # strain pattern e1, yields C11, C21, C31, C41, C51
     [ np.array([1,0,0,0,0,0]), _dec([(1,1), (2,1), (3,1), (4,1), (5,1)])],

     # strain pattern e3 + e4, yields C33, C44
     [ np.array([0,0,1,1,0,0]), _dec([(3,3), (4,4)]) ],

     [ np.array([0,0,0,0,0,1]), _dec([(6,6)]) ]
   ],

   'tetragonal': [
     # strain pattern e1+e4
     [ np.array([1,0,0,1,0,0]), _dec([(1,1), (2,1), (3,1), (6,1), (4,4)]) ],

     # strain pattern e3+e6
     [ np.array([0,0,1,0,0,1]), _dec([(3,3), (6,6)]) ]
   ],

   'orthorhombic': [
      # strain pattern e1+e4
      [ np.array([1,0,0,1,0,0]), _dec([(1,1), (2,1), (3,1), (4,4)]) ],

      # strain pattern e2+e5
      [ np.array([0,1,0,0,1,0]), _dec([(1,2), (2,2), (3,2), (5,5)]) ],

      # strain pattern e3+e6
      [ np.array([0,0,1,0,0,1]), _dec([(1,3), (2,3), (3,3), (6,6)]) ]
   ],

   'monoclinic': [
      # strain pattern e1+e4
      [ np.array([1,0,0,1,0,0]), _dec([(1,1), (2,1), (3,1), (4,4), (5,1), (6,4)]) ],

      # strain pattern e3+e6
      [ np.array([0,0,1,0,0,1]), _dec([(1,3), (2,3), (3,3), (5,3), (4,6), (6,6)]) ],

      # strain pattern e2
      [ np.array([0,1,0,0,0,0]), _dec([(1,2), (2,2), (3,2), (5,2)]) ],

      # strain pattern e5
      [ np.array([0,0,0,0,1,0]), _dec([(1,5), (2,5), (3,5), (5,5)]) ]
   ],

   'triclinic': [
      [ np.array([1,0,0,0,0,0]), _dec([(1,1), (2,1), (3,1), (4,1), (5,1), (6,1)])],
      [ np.array([0,1,0,0,0,0]), _dec([(1,2), (2,2), (3,2), (4,2), (5,2), (6,2)])],
      [ np.array([0,0,1,0,0,0]), _dec([(1,3), (2,3), (3,3), (4,3), (5,3), (6,3)])],
      [ np.array([0,0,0,1,0,0]), _dec([(1,4), (2,4), (3,4), (4,4), (5,4), (6,4)])],
      [ np.array([0,0,0,0,1,0]), _dec([(1,5), (2,5), (3,5), (4,5), (5,5), (6,5)])],
      [ np.array([0,0,0,0,0,1]), _dec([(1,6), (2,6), (3,6), (4,6), (5,6), (6,6)])],
   ]

   }

Cij_symmetry['hexagonal'] = Cij_symmetry['trigonal_high']
Cij_symmetry[None] = Cij_symmetry['triclinic']

strain_patterns['hexagonal'] = strain_patterns['trigonal_high']
strain_patterns['tetragonal_high'] = strain_patterns['tetragonal_low'] = strain_patterns['tetragonal']
strain_patterns[None] = strain_patterns['triclinic']


def generate_strained_configs(at0, symmetry='triclinic', N_steps=5, delta=1e-2):
    """
    Generate a sequence of strained configurations

    Parameters
    ----------
    a : ase.Atoms
        Bulk crystal configuration - both unit cell and atomic positions
        should be relaxed before calling this routine.
    symmetry : string
        Symmetry to use to determine which strain patterns are necessary.
        Default is 'triclininc', i.e. no symmetry.
    N_steps : int
        Number of atomic configurations to generate for each strain pattern.
        Default is 5. Absolute strain values range from -delta*N_steps/2
        to +delta*N_steps/2.
    delta : float
        Strain increment for analytical derivatives of stresses. Default 1e-2

    Returns
    -------
    Generator which yields a sequence of ase.Atoms instances corresponding
    to the minima set of strained conifugurations required to evaluate the
    full 6x6 C_ij matrix under the assumed symmetry.
    """

    if not symmetry in strain_patterns:
        raise ValueError('Unknown symmetry %s. Valid options are %s' % (symmetry, strain_patterns.keys()))

    for pindex, (pattern, fit_pairs) in enumerate(strain_patterns[symmetry]):
        for step in range(N_steps):
            strain = np.where(pattern == 1, delta*(step+1-(N_steps+1)/2.0), 0.0)
            at = at0.copy()
            if at0.get_calculator() is not None:
                at.set_calculator(at0.get_calculator())
            T = Voigt_6_to_full_3x3_strain(strain)
            at.set_cell(np.dot(T, at.cell.T).T, scale_atoms=False)
            at.positions[:] = np.dot(T, at.positions.T).T
            at.info['strain'] = T
            yield at


# Elastic constant calculation.

# Code adapted from elastics.py script, available from
# http://github.com/djw/elastic-constants
#
# Copyright (c) 2008, Dan Wilson
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of the copyright holder nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY DAN WILSON ''AS IS'' AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL DAN WILSON BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

def fit_elastic_constants(a, symmetry='triclinic', N_steps=5, delta=1e-2, 
                          optimizer=None, verbose=True, GPa=True, tag='test',
                          graphics=False, logfile=None, **kwargs):
    """
    Compute elastic constants by linear regression of stress vs. strain

    Parameters
    ----------
    a : ase.Atoms or list of ase.Atoms
        Either a single atomic configuration or a list of strained
        configurations. If a single configuration is given, it is
        passed :func:`generate_strained_configs` along with `symmetry`, `N_steps`,
        and `delta` to generate the set of strained configurations.
    symmetry : string
        Symmetry to use to determine which strain patterns are necessary.
        Default is 'triclininc', i.e. no symmetry.
    N_steps : int
        Number of atomic configurations to generate for each strain pattern.
        Default is 5. Absolute strain values range from -delta*N_steps/2
        to +delta*N_steps/2.    
    delta : float
        Strain increment for analytical derivatives of stresses.
        Default is 1e-2.
    optimizer : ase.optimizer.*
        Optimizer to use for atomic positions (internal degrees of freedom)
        for each applied strain. Initial config `a` is not optimised, and should
        already have relaxed cell and atomic positions. Does not optimize atomic
        positions if set to None.
    verbose : bool
        If True, print additional infomation about the quality of fit
        and summarise results of C_ij and estimated errors. Default True.
    GPa : bool
        If True, convert to GPa
    graphics : bool
        If True, use :mod:`matplotlib.pyplot` to plot the stress vs. strain
        curve for each C_ij component fitted. Default True. 
    logfile : bool
        Log file to write optimizer output to. Default None (i.e. suppress
        output).
    **kwargs : dict
        Additional arguments to pass to `optimizer.run()` method e.g. `fmax`.

    Returns
    -------
    C : array_like
        6x6 matrix of the elastic constants in Voigt notation.
    C_err : array_like
        If scipy.stats module is available then error estimates for each C_ij
        component are obtained from the accuracy of the linear regression.
        Otherwise an array of np.zeros((6,6)) is returned.
    
    Notes
    -----

    Code originally adapted from elastics.py script, available from
    http://github.com/djw/elastic-constants
    """

    if isinstance(a, Atoms):
        if a.constraints:
            raise ValueError('Atoms passed to fit_elastic_constants() '
                             'has constraints attached')
        # we've been passed a single Atoms object: use it to generate
        # set of strained configurations according to symmetry
        strained_configs = generate_strained_configs(a, symmetry, N_steps, delta)
    else:
        # assume we've been passed a list of strained configs
        strained_configs = a

    if graphics:
        import matplotlib.pyplot as plt

    def do_fit(index1, index2, stress, strain, patt):
        if verbose:
            print('Fitting C_%d%d' % (index1+1, index2+1))
            print('Strain %r' % strain[:,index2])
            print('Stress %r GPa' % (stress[:,index1]/units.GPa))

        if scipy_stats is not None:
            cijFitted, intercept,r,tt,stderr = scipy_stats.linregress(strain[:,index2],stress[:,index1])
        else:
            cijFitted, intercept = np.polyfit(strain[:,index2],stress[:,index1], 1)
            r, tt, stderr = 0.0, None, 0.0

        if verbose:
            # print info about the fit
            print('Cij (gradient) / GPa    :    ',cijFitted/units.GPa)
            if scipy_stats is not None:
                print('Error in Cij / GPa      :    ', stderr/units.GPa)
                if abs(r) > 0.9:
                    print('Correlation coefficient :    ',r)
                else:
                    print('Correlation coefficient :    ',r, '     <----- WARNING')
            else:
                print('scipy.stats not available, no error estimation performed')

        if graphics:
            # position this plot in a 6x6 grid
            sp = plt.subplot(6,6,6*index1+(index2+1))
            sp.set_axis_on()

            # change the labels on the axes
            xlabels = sp.get_xticklabels()
            plt.setp(xlabels,'rotation',90,fontsize=7)
            ylabels = sp.get_yticklabels()
            plt.setp(ylabels,fontsize=7)

            # colour the plot depending on the strain pattern
            colourDict = {0: '#BAD0EF', 1:'#FFCECE', 2:'#BDF4CB', 3:'#EEF093',4:'#FFA4FF',5:'#75ECFD'}
            sp.set_axis_bgcolor(colourDict[patt])

            # plot the data
            plt.plot([strain[0,index2],strain[-1,index2]],
                     [(cijFitted*strain[0,index2]+intercept)/units.GPa,
                      (cijFitted*strain[-1,index2]+intercept)/units.GPa])
            plt.plot(strain[:,index2],stress[:,index1]/units.GPa,'ro')
            plt.xticks(strain[:,index2])

        return cijFitted, stderr

    if not symmetry in strain_patterns:
        raise ValueError('Unknown symmetry %s. Valid options are %s' % (symmetry, strain_patterns.keys()))

    # There are 21 independent elastic constants
    Cijs = {}
    Cij_err = {}

    # Construct mapping from (i,j) to index into Cijs in range 1..21
    # (upper triangle only to start with)
    Cij_map = {}
    Cij_map_sym = {}
    for i in range(6):
        for j in range(i,6):
            Cij_map[(i,j)] = Cij_symmetry[None][i,j]
            Cij_map_sym[(i,j)] = Cij_symmetry[symmetry][i,j]

    # Reverse mapping, index 1..21 -> tuple (i,j) with i, j in range 0..5
    Cij_rev_map = dict(zip(Cij_map.values(), Cij_map.keys()))

    # Add the lower triangle to Cij_map, e.g. C21 = C12
    for (i1,i2) in Cij_map.copy().keys():
        Cij_map[(i2,i1)] = Cij_map[(i1,i2)]
        Cij_map_sym[(i2,i1)] = Cij_map_sym[(i1,i2)]


    N_pattern = len(strain_patterns[symmetry])
    configs = iter(strained_configs)

    strain = np.zeros((N_pattern, N_steps, 6))
    stress = np.zeros((N_pattern, N_steps, 6))

    if graphics:
        fig = plt.figure(num=1, figsize=(9.5,8),facecolor='white')
        fig.clear()
        fig.subplots_adjust(left=0.07,right=0.97,top=0.97,bottom=0.07,wspace=0.5,hspace=0.5)

        for index1 in range(6):
            for index2 in range(6):
                # position this plot in a 6x6 grid
                sp = plt.subplot(6,6,6*(index1)+index2+1)
                sp.set_axis_off()
                plt.text(0.4,0.4, "n/a")

    # Fill in strain and stress arrays from config Atoms list
    with open('Detail-'+tag+'.txt', 'w') as f:
        for pattern_index, (pattern, fit_pairs) in enumerate(strain_patterns[symmetry]):
            for step in range(N_steps):
                at = next(configs)
                t0 = time()
                E0 = at.get_potential_energy()
                if optimizer is not None:
                    optimizer(at, logfile=logfile).run(**kwargs)
                else:
                    # update position
                    pos = read('geo_end.gen').get_positions()
                    at.set_positions(pos)
                E1 = at.get_potential_energy()
                fmax = np.abs(at.get_forces()).max()
                t1 = time()-t0
                strs = "\n{:8s} {:2d}/{:2d}\n".format(tag, pattern_index, step)
                strs += "Eng:  {:.4f} -> {:.4f}, ".format(E0, E1)
                strs += "dE: {:.4f}\n".format(E1-E0)
                strs += "fmax: {:.5f}\n".format(fmax)
                strs += "time: {:.1f}\n".format(t1)
                strain_info = full_3x3_to_Voigt_6_strain(at.info['strain'])
                stress_info = at.get_stress()
                strain[pattern_index, step, :] = strain_info
                stress[pattern_index, step, :] = stress_info
                #print("Cell\n", at.get_cell())
                strs += "Strain {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f}\n".format(*strain_info)
                strs += "Stress (GPa) {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f}\n".format(*(stress_info/units.GPa))
                f.write(strs)               
                print(strs)

    # Do the linear regression
    for pattern_index, (pattern, fit_pairs) in enumerate(strain_patterns[symmetry]):
        for (index1, index2) in fit_pairs:
            fitted, err = do_fit(index1, index2,
                                 stress[pattern_index,:,:],
                                 strain[pattern_index,:,:],
                                 pattern_index)

            index = abs(Cij_map_sym[(index1, index2)])

            if not index in Cijs:
                if verbose:
                    print('Setting C%d%d (%d) to %f +/- %f' %
                          (index1+1, index2+1, index, fitted, err))
                Cijs[index] = [fitted]
                Cij_err[index] = [err]
            else:
                if verbose:
                    print('Updating C%d%d (%d) with value %f +/- %f' %
                          (index1+1, index2+1, index, fitted, err))
                Cijs[index].append(fitted)
                Cij_err[index].append(err)
            if verbose:
                print('\n')


    C = np.zeros((6,6))
    C_err = np.zeros((6,6))
    C_labels = np.zeros((6,6),dtype='S4')
    C_labels[:] = '    '

    # Convert lists to mean
    for k in Cijs:
        Cijs[k] = np.mean(Cijs[k])

    # Combine statistical errors
    for k, v in Cij_err.items():
        Cij_err[k] = np.sqrt(np.sum(np.array(v)**2))/np.sqrt(len(v))

    if symmetry.startswith('trigonal'):
        # Special case for trigonal lattice: C66 = (C11 - C12)/2
        Cijs[Cij_map[(5,5)]] = 0.5*(Cijs[Cij_map[(0,0)]]-Cijs[Cij_map[(0,1)]])
        Cij_err[Cij_map[(5,5)]] = np.sqrt(Cij_err[Cij_map[(0,0)]]**2 + Cij_err[Cij_map[(0,1)]]**2)

    # Generate the 6x6 matrix of elastic constants
    # - negative values signify a symmetry relation
    for i in range(6):
        for j in range(6):
            index = Cij_symmetry[symmetry][i,j]
            if index > 0:
                C[i,j] = Cijs[index]
                C_err[i,j] = Cij_err[index]
                ii, jj = Cij_rev_map[index]
                C_labels[i,j] = ' C%d%d' % (ii+1,jj+1)
                C_err[i,j] = Cij_err[index]
            elif index < 0:
                C[i,j] = -Cijs[-index]
                C_err[i,j] = Cij_err[-index]
                ii, jj = Cij_rev_map[-index]
                C_labels[i,j] = '-C%d%d' % (ii+1, jj+1)
    if GPa:
        C /= units.GPa
        C_err /= units.GPa

    if verbose:
        print(np.array2string(C_labels).replace("'",""))
        print('\n = \n')
        print(np.array2string(C/units.GPa,
                              suppress_small=True,
                              precision=2))
        # Summarise the independent components of C_ij matrix
        printed = {}
        for i in range(6):
            for j in range(6):
                index = Cij_symmetry[symmetry][i,j]
                if index <= 0 or index in printed: continue
                print('C_%d%d = %-4.2f +/- %-4.2f GPa' %
                      (i+1, j+1, C[i,j]/units.GPa, C_err[i,j]/units.GPa))
                printed[index] = 1

    return C, C_err


def youngs_modulus(C, l):
    """
    Calculate approximate Youngs modulus E_l from 6x6 elastic constants matrix C_ij

    This is the modulus for loading in the l direction. For the exact answer, taking
    into account elastic anisotropuy, rotate the C_ij matrix to the correct frame,
    compute the compliance matrix::

       C = ...  # 6x6 C_ij matrix in crystal frame
       A = ...  # rotation matrix
       Cr = rotate_elastic_constants(C, A)
       S = np.inv(Cr)
       E_x = 1/S[0, 0]  # Young's modulus for a pull in x direction
       E_y = 1/S[1, 1]  # Young's modulus for a pull in y direction
       E_z = 1/S[0, 0]  # Young's modulus for a pull in z direction

    Notes
    -----

    Formula is from W. Brantley, Calculated elastic constants for stress problems associated 
    with semiconductor devices. J. Appl. Phys., 44, 534 (1973).
    """  

    S = inv(C)        # Compliance matrix
    lhat = l/norm(l)  # Normalise directions

    # Youngs modulus in direction l, ratio of stress sigma_l 
    # to strain response epsilon_l
    E = 1.0/(S[0,0] - 2.0*(S[0,0]-S[0,1]-0.5*S[3,3])*(lhat[0]*lhat[0]*lhat[1]*lhat[1] +
         lhat[1]*lhat[1]*lhat[2]*lhat[2] +
         lhat[0]*lhat[0]*lhat[2]*lhat[2]))
    return E


def poisson_ratio(C, l, m):
    """
    Calculate approximate Poisson ratio \nu_{lm} from 6x6 elastic constant matrix C_{ij}

    This is the response in `m` direction to pulling in `l` direction. Result is dimensionless.

    Notes
    -----

    Formula is from W. Brantley, Calculated elastic constants for stress problems associated 
    with semiconductor devices. J. Appl. Phys., 44, 534 (1973).
    """
    
    S = inv(C)        # Compliance matrix
    lhat = l/norm(l)  # Normalise directions
    mhat = m/norm(m)

    # Poisson ratio v_lm: response in m direction to strain in 
    # l direction, v_lm = - epsilon_m/epsilon_l
    v = -((S[0,1] + (S[0,0]-S[0,1]-0.5*S[3,3])*(lhat[0]*lhat[0]*mhat[0]*mhat[0] +
         lhat[1]*lhat[1]*mhat[1]*mhat[1] +
         lhat[2]*lhat[2]*mhat[2]*mhat[2])) / 
         (S[0,0] - 2.0*(S[0,0]-S[0,1]-0.5*S[3,3])*(lhat[0]*lhat[0]*lhat[1]*lhat[1] + 
         lhat[1]*lhat[1]*lhat[2]*lhat[2] + 
         lhat[0]*lhat[0]*lhat[2]*lhat[2])))
    return v


def elastic_moduli(C, l=np.array([1, 0, 0]), R=None, tol=1e-6):
    """
    Calculate elastic moduli from 6x6 elastic constant matrix C_{ij}.
    
    The elastic moduli calculated are: Young's muduli, Poisson's ratios,
    shear moduli, bulk mudulus and bulk mudulus tensor.
    
    If a direction l is specified, the system is rotated to have it as its
    x direction (see Notes for details). If R is specified the system is
    rotated according to it.
    
    Parameters
    ----------
    C : array_like
        6x6 stiffness matrix (Voigt notation).
    l : array_like, optional
        3D direction vector for pull (the default is the x direction
        of the original system)
    R : array_like, optional
        3x3 rotation matrix.
    
    Returns
    -------
    E : array_like
        Young's modulus for a stress in each of the three directions
        of the rotated system.
    nu : array_like
         3x3 matrix with Poisson's ratios.
    Gm : array_like
         3x3 matrix with shear moduli.
    B : float
        Bulk modulus.
    K : array_like
        3x3 matrix with bulk modulus tensor.
    
    Other Parameters
    ----------------
    tol : float, optional
        tolerance for checking validity of rotation and comparison
        of vectors.
    
    Notes
    ---
    It works by rotating the elastic constant tensor to the desired
    direction, so it should be valid for arbitrary crystal structures.
    If only l is specified there is an infinite number of possible
    rotations. The chosen one is a rotation along the axis orthogonal
    to the plane defined by the vectors (1, 0, 0) and l.
    
    Bulk modulus tensor as defined in
    O. Rand and V. Rovenski, "Analytical Methods in Anisotropic
    Elasticity", Birkh\"auser (2005), pp. 71.
    
    """
    
    if R is not None:
        R = np.asarray(R)

        # Is this a rotation matrix?
        if np.sometrue(np.abs(np.dot(np.array(R), np.transpose(np.array(R))) - 
                              np.eye(3, dtype=float)) > tol):
            raise RuntimeError('Matrix *R* does not describe a rotation.')
    else:
        u_a = np.array([1, 0, 0])
        u_b = l/norm(l)  # Normalise directions
        R = np.eye(3)

        if not np.allclose(l, u_a, rtol=tol, atol=tol):
            u_v = np.cross(u_a, u_b)
            u_v_mat = np.array([[ 0,      -u_v[2], u_v[1]],
                                [ u_v[2], 0,      -u_v[0]],
                                [-u_v[1], u_v[0], 0]])

            R = R + u_v_mat + \
                np.dot(u_v_mat, u_v_mat) * (1 - np.dot(u_a, u_b)) / \
                np.linalg.norm(u_v)**2

    Cr = rotate_elastic_constants(C, R)
    S = np.linalg.inv(Cr)

    # Young's modulus for a stress in \alpha direction; \alpha = x, y, z
    E = np.zeros(3)
    E[0] = 1/S[0, 0]
    E[1] = 1/S[1, 1]
    E[2] = 1/S[2, 2]

    # Poisson's ratio ($\alpha$, $\beta$); $\alpha$, $\beta$ = x, y, z
    nu = np.array([
            [-1,               -S[1, 0]/S[0, 0], -S[2, 0]/S[0, 0]],
            [-S[0, 1]/S[1, 1], -1,               -S[2, 1]/S[1, 1]],
            [-S[0, 2]/S[2, 2], -S[1, 2]/S[2, 2], -1]
            ])

    # Shear modulus ($\alpha$, $\beta$); $\alpha$, $\beta$ = x, y, z
    G = np.zeros(3)
    G[0] = 1/S[3, 3]  # Shear modulus yz
    G[1] = 1/S[4, 4]  # Shear modulus zx
    G[2] = 1/S[5, 5]  # Shear modulus xy
    Gm = np.array([
            [E[0]/4,  G[2],    G[1]],
            [G[2],    E[1]/4,  G[0]],
            [G[1],    G[0],    E[2]/4]
            ])

    # Bulk modulus
    B = 1/np.sum(S[0:3, 0:3])
    
    # Bulk modulus tensor
    Crt = Voigt_6x6_to_full_3x3x3x3(Cr)
    K = np.einsum('ijkk', Crt)

    return E, nu, Gm, B, K

def elastic_properties(C):
    """
    A quick summary of elastic properties from the 6*6 matrix
    """
    Kv = C[:3,:3].mean()
    Gv = (C[0,0]+C[1,1]+C[2,2] - (C[0,1]+C[1,2]+C[2,0]) + 3*(C[3,3]+C[4,4]+C[5,5]))/15
    Ev = 1/((1/(3*Gv))+(1/(9*Kv)))
    vv  = 0.5*(1-((3*Gv)/(3*Kv+Gv))); 

    S = np.linalg.inv(C)
    Kr = 1/((S[0,0]+S[1,1]+S[2,2])+2*(S[0,1]+S[1,2]+S[2,0])) 
    Gr = 15/(4*(S[0,0]+S[1,1]+S[2,2])-4*(S[0,1]+S[1,2]+S[2,0])+3*(S[3,3]+S[4,4]+S[5,5])) 
    Er = 1/((1/(3*Gr))+(1/(9*Kr))) 
    vr = 0.5*(1-((3*Gr)/(3*Kr+Gr))) 

    Kh = (Kv+Kr)/2    
    Gh = (Gv+Gr)/2    
    Eh = (Ev+Er)/2    
    vh = (vv+vr)/2   
    return Kv, Gv, Ev, vv, Kr, Gr, Er, vr, Kh, Gh, Eh, vh


