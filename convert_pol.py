#"This script calculates the polarization along cartesian coordinates. 
# It should be used when the unitcell is non-orthogonal. Workflow:
# 1. Reads the lattice vectors in geomerty file that is in aims format
# 2. Convert to reciprocal space
# 3. Normalize the rec. lattice vectors
# 4. Find the linear coefficients of the rec. lattice vectors that gives 
#    the three cartesian directions (3 equations, 3 unknowns).
# 5. polarization elements of the non-orthogonal unitcell are collected   
#    from aims output
# 6. px,py and pz are obtained
#"
from ase.build import bulk
from math import sqrt
from ase.io import read, write
import numpy as np
from numpy import dot, array, append
from sklearn import preprocessing

def split_line(lines):
    """Split input line"""
    line_array = array(lines.strip().split(" "))
    line_vals = line_array[line_array != ""]
    return line_vals


geo = read("geometry.in", 0, "aims")
filename = "aims.out"
#data = open(filename)
#lines = data.readlines()
# c=3 #change
#for line in lines:
#    if (
#        line.rfind(
#            "- Directive    1 in direction of rec. latt. vec.  1 yields the full polarization      :"
#        )
#        != -1
#    ):
#        p_x = np.float64(split_line(line)[-2:-1])
#    if (
#        line.rfind(
#            "- Directive    2 in direction of rec. latt. vec.  2 yields the full polarization      :"
#        )
#        != -1
#    ):
#        p_y = np.float64(split_line(line)[-2:-1])
#    if (
#        line.rfind(
#            "- Directive    3 in direction of rec. latt. vec.  3 yields the full polarization      :"
#        )
#        != -1
#    ):
#        p_z = np.float64(split_line(line)[-2:-1])
# Find the reciprocal lattice vectors
rec_lat = geo.get_reciprocal_cell()
print("reciprocal lattice vectors are: \n", rec_lat)


# l2-normalize the samples (rows).
X_normalized = preprocessing.normalize(rec_lat, norm='l2')
print("normalized reciprocal lattice vectors are: \n",X_normalized)
# Solving 3 equations, 3 unknowns:
x = np.linalg.solve(X_normalized, np.array([1, 0, 0]))
y = np.linalg.solve(X_normalized, np.array([0, 1, 0]))
z = np.linalg.solve(X_normalized, np.array([0, 0, 1]))
print(np.allclose(np.dot(X_normalized, x), np.array([1, 0, 0])))
# print('1?', np.dot(normalized_v,x))
# print('1?', np.dot(normalized_v,y))
# print('1?', np.dot(normalized_v,z))
exit()
R = np.array([x, y, z])
p = np.array([p_x, p_y, p_z])
print("Polarization of non-orthogonal is : \n", p)
print("R is: \n", R)
final = np.dot(p, R)

print("P_x is", final[0])
print("P_y is", final[1])
print("P_z is", final[2])
# now collect p1 p2 p3
