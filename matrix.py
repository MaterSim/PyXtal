import numpy as np
from numpy import matrix
from numpy.linalg import eig
from numpy.linalg import det
import math
import random
from pymatgen.core.operations import SymmOp
from math import pi
rad = pi/180.
deg = 180./pi

#returns a 3x3 rotation matrix with random angle and orientation
'''def aa2matrix(axis='random',angle='random', radians=True):
    #Return a rotation matrix from an axis and angle
    #If axis and/or angle is not provided, a random value is chosen
    #TODO: Stop relying on SymmOp method - inaccurate
    #If axis is provided, it is assumed to be unit length
    if axis == 'random':
        phi = random.random() * 2. * math.pi
        N = random.random()
        theta = math.acos(2.*N - 1.)
        x = math.cos(phi) * math.sin(theta)
        y = math.sin(phi) * math.sin(theta)
        z = math.cos(theta)
    else:
        x = axis[0]
        y = axis[1]
        z = axis[2]
    if angle == 'random':
    	angle = random.random() * 2 * math.pi
    else:
        if radians is False:
            angle *= rad
    print(angle)
    return SymmOp.from_origin_axis_angle([0,0,0], [x,y,z], angle, angle_in_radians=True).rotation_matrix'''

def aa2matrix(axis, angle, radians=True):
    '''
    Given an axis and an angle, return a 3x3 rotation matrix
    Based on:
    https://en.wikipedia.org/wiki/Rotation_matrix#Axis_and_angle
    '''
    #Convert to radians if necessary
    if radians is not True:
        angle *= rad
    #Ensure axis is a unit vector
    axis = axis / np.linalg.norm(axis)
    #Define quantities which are reused
    x = axis[0]
    y = axis[1]
    z = axis[2]
    c = math.cos(angle)
    s = math.sin(angle)
    C = 1 - c
    #Define the rotation matrix
    Q = np.zeros([3,3])
    Q[0][0] = x*x*C + c
    Q[0][1] = x*y*C - z*s
    Q[0][2] = x*z*C + y*s
    Q[1][0] = y*x*C + z*s
    Q[1][1] = y*y*C + c
    Q[1][2] = y*z*C - x*s
    Q[2][0] = z*x*C - y*s
    Q[2][1] = z*y*C + x*s
    Q[2][2] = z*z*C + c
    return Q

    #def matrix2aa(m, radians=True):
    '''
    Return the axis and angle from a rotation matrix.
    '''
    '''#TODO: improve functionality with reflection and inversion matrices
    if type(m) == SymmOp:
        m = m.rotation_matrix
    if np.linalg.det(m) < 0:
        m *= -1.
    value = m[0][0]+m[1][1]+m[2][2]
    #Avoid problems caused by numerical errors
    if value > 1 or value < -1: value = np.round(value)
    angle = math.acos(( value - 1)/2)
    if radians is False: angle *= deg
    x = (m[2][1] - m[1][2])/np.sqrt((m[2][1] - m[1][2])**2+(m[0][2] - m[2][0])**2+(m[1][0] - m[0][1])**2)
    y = (m[0][2] - m[2][0])/np.sqrt((m[2][1] - m[1][2])**2+(m[0][2] - m[2][0])**2+(m[1][0] - m[0][1])**2)
    z = (m[1][0] - m[0][1])/np.sqrt((m[2][1] - m[1][2])**2+(m[0][2] - m[2][0])**2+(m[1][0] - m[0][1])**2)
    return [x,y,z], angle'''

def matrix2aa(m, radians=True):
    '''
    Return the axis and angle from a rotation matrix.
    '''
    #TODO: improve functionality with reflection and inversion matrices
    if type(m) == SymmOp:
        m = m.rotation_matrix
    #eigenvalues, eivenvectors = np.linalg.eig(m)


#Functions from 
def rotation_matrix(angle, direction, point=None):
    """Return matrix to rotate about axis defined by point and direction.

    >>> R = rotation_matrix(math.pi/2, [0, 0, 1], [1, 0, 0])
    >>> numpy.allclose(numpy.dot(R, [0, 0, 0, 1]), [1, -1, 0, 1])
    True
    >>> angle = (random.random() - 0.5) * (2*math.pi)
    >>> direc = numpy.random.random(3) - 0.5
    >>> point = numpy.random.random(3) - 0.5
    >>> R0 = rotation_matrix(angle, direc, point)
    >>> R1 = rotation_matrix(angle-2*math.pi, direc, point)
    >>> is_same_transform(R0, R1)
    True
    >>> R0 = rotation_matrix(angle, direc, point)
    >>> R1 = rotation_matrix(-angle, -direc, point)
    >>> is_same_transform(R0, R1)
    True
    >>> I = numpy.identity(4, numpy.float64)
    >>> numpy.allclose(I, rotation_matrix(math.pi*2, direc))
    True
    >>> numpy.allclose(2, numpy.trace(rotation_matrix(math.pi/2,
    ...                                               direc, point)))
    True

    """
    sina = math.sin(angle)
    cosa = math.cos(angle)
    direction = unit_vector(direction[:3])
    # rotation matrix around unit vector
    R = numpy.diag([cosa, cosa, cosa])
    R += numpy.outer(direction, direction) * (1.0 - cosa)
    direction *= sina
    R += numpy.array([[ 0.0,         -direction[2],  direction[1]],
                      [ direction[2], 0.0,          -direction[0]],
                      [-direction[1], direction[0],  0.0]])
    M = numpy.identity(4)
    M[:3, :3] = R
    if point is not None:
        # rotation not around origin
        point = numpy.array(point[:3], dtype=numpy.float64, copy=False)
        M[:3, 3] = point - numpy.dot(R, point)
    return M

def rotation_from_matrix(matrix):
    """Return rotation angle and axis from rotation matrix.

    >>> angle = (random.random() - 0.5) * (2*math.pi)
    >>> direc = numpy.random.random(3) - 0.5
    >>> point = numpy.random.random(3) - 0.5
    >>> R0 = rotation_matrix(angle, direc, point)
    >>> angle, direc, point = rotation_from_matrix(R0)
    >>> R1 = rotation_matrix(angle, direc, point)
    >>> is_same_transform(R0, R1)
    True

    """
    R = numpy.array(matrix, dtype=numpy.float64, copy=False)
    R33 = R[:3, :3]
    # direction: unit eigenvector of R33 corresponding to eigenvalue of 1
    w, W = numpy.linalg.eig(R33.T)
    i = numpy.where(abs(numpy.real(w) - 1.0) < 1e-8)[0]
    if not len(i):
        raise ValueError('no unit eigenvector corresponding to eigenvalue 1')
    direction = numpy.real(W[:, i[-1]]).squeeze()
    # point: unit eigenvector of R33 corresponding to eigenvalue of 1
    w, Q = numpy.linalg.eig(R)
    i = numpy.where(abs(numpy.real(w) - 1.0) < 1e-8)[0]
    if not len(i):
        raise ValueError('no unit eigenvector corresponding to eigenvalue 1')
    point = numpy.real(Q[:, i[-1]]).squeeze()
    point /= point[3]
    # rotation angle depending on direction
    cosa = (numpy.trace(R33) - 1.0) / 2.0
    if abs(direction[2]) > 1e-8:
        sina = (R[1, 0] + (cosa-1.0)*direction[0]*direction[1]) / direction[2]
    elif abs(direction[1]) > 1e-8:
        sina = (R[0, 2] + (cosa-1.0)*direction[0]*direction[2]) / direction[1]
    else:
        sina = (R[2, 1] + (cosa-1.0)*direction[1]*direction[2]) / direction[0]
    angle = math.atan2(sina, cosa)
    return angle, direction, point
    

#Test Functionality
if __name__ == "__main__":
#----------------------------------------------------
    identity = np.matrix([[1,0,0],[0,1,0],[0,0,1]])
    inverse = np.matrix([[-1,0,0],[0,-1,0],[0,0,-1]])
    x2 = np.matrix([[1,0,0],[0,-1,0],[0,0,-1]])
    y2 = np.matrix([[-1,0,0],[0,1,0],[0,0,-1]])
    mx = np.matrix([[-1,0,0],[0,1,0],[0,0,1]])
    my = np.matrix([[1,0,0],[0,-1,0],[0,0,1]])
    #Note: 2bar is just inverse
    x4 = [[],[],[]]
    x_4 = [[],[],[]]
    x3 = np.matrix([[],[],[]])
    x_3 = np.matrix([[],[],[]])

    '''for x in [identity, inverse, x2, y2, mx, my, x3, x_3]:
        print('------------------')
        print('Matrix:')
        print(x)
        print('Eigenvalues:')
        print(eig(x)[0])
        print('Eigenvectors:')
        print(eig(x)[1])'''

    print( aa2matrix([1,0,0],180, radians=False) )