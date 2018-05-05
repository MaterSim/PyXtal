import numpy as np
import math
import random
import pymatgen

#returns a 3x3 rotation matrix with random angle and orientation
def aa2matrix(axis='random',angle='random', radians=True):
    '''
    Return a rotation matrix from an axis and angle
    If axis and/or angle is not provided, a random value is chosen
    '''
    #If axis is provided, it is assumed to be unit length
    if axis == 'random':
        phi = random.random() * 2. * math.pi
        N = random.random()
        theta = math.acos(2.*N - 1.)
        x = math.cos(phi) * math.sin(theta)
        y = math.sin(phi) * math.sin(theta)
        z = math.cos(theta)
    if angle == 'random':
    	angle = random.random() * 2 * math.pi
    else:
        if radians is False:
            angle *= 180./math.pi
    return pymatgen.core.operations.SymmOp.from_origin_axis_angle([0,0,0], [x,y,z], angle, angle_in_radians=True).rotation_matrix

def matrix2aa(m):
    '''
    Return the axis and angle from a rotation matrix.
    '''
    #TODO: improve functionality with reflection and inversion matrices
    angle = math.acos(( m[0][0] + m[1][1] + m[2][2] - 1)/2)
    x = (m[2][1] - m[1][2])/np.sqrt((m[2][1] - m[1][2])**2+(m[0][2] - m[2][0])**2+(m[1][0] - m[0][1])**2)
    y = (m[0][2] - m[2][0])/np.sqrt((m[2][1] - m[1][2])**2+(m[0][2] - m[2][0])**2+(m[1][0] - m[0][1])**2)
    z = (m[1][0] - m[0][1])/np.sqrt((m[2][1] - m[1][2])**2+(m[0][2] - m[2][0])**2+(m[1][0] - m[0][1])**2)
    return [x,y,z], angle

#Test Functionality
#----------------------------------------------------
'''from structure import *
for x in range(100):
    mat1 = random_matrix()
    paras1 = matrix2para(mat1)
    mat2 = para2matrix(paras1, format='upper')
    paras2 = matrix2para(mat2)
    if not np.allclose(paras1, paras2):
        print(paras1)
        print(paras2)'''
