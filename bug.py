from structure import *
import numpy as np

#It is stange that check_wyckoff returns false for the following set.
#should be 8j position for spg 97
coor = np.array(
[[0.85540127, 0.35540127, 0.25],
 [0.14459873, 0.64459873, 0.25],
 [0.64459873, 0.85540127, 0.25],
 [0.35540127, 0.14459873, 0.25],
 [0.35540127, 0.85540127, 0.75],
 [0.64459873, 0.14459873, 0.75],
 [0.14459873, 0.35540127, 0.75],
 [0.85540127, 0.64459873, 0.75]]
)
print(check_wyckoff_position(coor, 97))

#This set is to check the numerical tolerance of check_wyckoff
#coor = np.array(
#    [[-2.77555756e-17, -2.77555756e-17,  1.29634884e+00],
#     [-5.00000000e-01, -5.00000000e-01,  7.96348839e-01]])
#print(check_wyckoff_position(coor, 79))
#coor = np.around(coor*1000)/1000
#print(check_wyckoff_position(coor, 79))

