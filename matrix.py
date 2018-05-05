import numpy as np
import math
import random
from pymatgen.core.operations import SymmOp

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
    else:
        x = axis[0]
        y = axis[1]
        z = axis[2]
    if angle == 'random':
    	angle = random.random() * 2 * math.pi
    else:
        if radians is False:
            angle *= 180./math.pi
    return SymmOp.from_origin_axis_angle([0,0,0], [x,y,z], angle, angle_in_radians=True).rotation_matrix

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
if __name__ == "__main__":
#----------------------------------------------------
    m1 = aa2matrix()
    aa1 = matrix2aa(m1)
    print(aa1)
    m2 = aa2matrix(axis=aa1[0], angle=aa1[1])
    aa2 = matrix2aa(m2)
    print(aa2)
    m3 = aa2matrix(axis=aa2[0], angle=aa2[1])
    aa3 = matrix2aa(m3)
    print(aa3)
