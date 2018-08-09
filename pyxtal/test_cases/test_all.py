"""
Test script for pyXtal version 0.1dev. Tests core functions for all modules.
"""
import sys

#Check if module and classes work correctly
def passed():
    global failed_module
    global failed
    if failed_module is False and failed is False:
        return True
    else:
        return False

#Reset flags for module and class
def reset():
    global failed_module
    global failed
    failed_module = False
    failed = False

#Set flags for package, module, class if error occurs
def fail():
    global failed_package
    global failed_module
    global failed
    failed_package = True
    failed_module = True
    failed = True
    print("    Error")

#Print whether module passed or failed
def check():
    if passed():
        pass#print("Success!")
    else:
        print("Failed module.")

#Call at end of script, or if module fails
def end():
    print("===")
    if failed_package is False:
        print("All modules passed!")
    else:
        print("One or more modules failed. Try reinstalling the package.")
    sys.exit(0)

#If importing fails
def import_fail():
    print("---Error: could not import---")
    fail()
    end()

print("====== Testing functionality for pyXtal version 0.1dev ======")

failed_package = False #Record if errors occur at any level

print("Importing sys...")
try:
    import sys
    print("Success!")
except:
    print("Error: could not import sys. Try reinstalling Python.")
    sys.exit(0)

print("Importing numpy...")
try:
    import numpy as np
    print("Success!")
except:
    print("Error: could not import numpy. Try reinstalling the package.")
    sys.exit(0)

print("Importing pymatgen...")
try:
    import pymatgen
    print("Success!")
except:
    print("Error: could not import pymatgen. Try reinstalling the package.")
    sys.exit(0)

try:
    from pymatgen.core.operations import SymmOp
except:
    print("Error: could not import SymmOp object from pymatgen. Try reinstalling the package.")
    sys.exit(0)

print("Importing pandas...")
try:
    import pandas
    print("Success!")
except:
    print("Error: could not import pandas. Try reinstalling the package.")
    sys.exit(0)

print("Importing spglib...")
try:
    import spglib
    print("Success!")
except:
    print("Error: could not import spglib. Try reinstalling the package.")
    sys.exit(0)

print("Importing ase...")
try:
    import ase
    print("Success!")
except:
    print("Error: could not import ase. Try reinstalling the package.")
    sys.exit(0)

print("Importing pyxtal...")
try:
    import pyxtal
    print("Success!")
except:
    print("Error: could not import pyxtal. Try reinstalling the package.")
    sys.exit(0)

print("=== Testing modules ===")

#=====database.element=====
print("pyxtal.database.element")
reset()
try:
    import pyxtal.database.element
except:
    import_fail()

print("  class Element")
try:
    from pyxtal.database.element import Element
except:
    fail()
    print("  Error: Could not import class Element")
if passed():
    for i in range(1, 95):
        if passed():
            try:
                ele = Element(i)
            except:
                fail()
                print("    Error: Could not access Element # "+str(i))
            try:
                y = ele.sf
                y = ele.z
                y = ele.short_name
                y = ele.long_name
                y = ele.valence
                y = ele.valence_electrons
                y = ele.covalent_radius
                y = ele.vdw_radius
                y = ele.get_all(0)
            except:
                fail()
                print("    Error accessing attribute for element # "+str(i))
            try:
                ele.all_z()
                ele.all_short_names()
                ele.all_long_names()
                ele.all_valences()
                ele.all_valence_electrons()
                ele.all_covalent_radii()
                ele.all_vdw_radii()
            except:
                fail()
                print("    Error accessing class methods")

check()

#=====database.hall=====
print("pyxtal.database.hall")
reset()
try:
    import pyxtal.database.hall
except:
    import_fail()

print("  function hall_from_hm")
try:
    from pyxtal.database.hall import hall_from_hm
except:
    import_fail()

if passed():
    for i in range(1, 230):
        if passed():
            try:
                hall_from_hm(i)
            except:
                fail()
                print("Could not access hm # "+str(i))

check()

#=====database.layergroup=====
print("pyxtal.database.layergroup")
reset()
try:
    import pyxtal.database.layergroup
except:
    import_fail()

print("  class Layergroup")
try:
    from pyxtal.database.layergroup import Layergroup
except:
    import_fail()

if passed():
    for i in range(1, 81):
        if passed():
            try:
                lgp = Layergroup(i)
            except:
                fail()
                print("    Error: Could not access layer group # "+str(i))
            try:
                lgp.input
                lgp.lg
                lgp.symbol
                lgp.sgnumber
                lgp.permutation
            except:
                fail()
                print("    Error accessing accessing attribute for layer group # "+str(i))

check()

#=====operations=====
print("pyxtal.operations")
reset()
try:
    import pyxtal.operations
except:
    import_fail()

print("  function random_vector")
try:
    from pyxtal.operations import random_vector
except:
    import_fail()

if passed():
    try:
        for i in range(10):
            random_vector()
    except:
        fail()

check()

print("  function angle")
try:
    from pyxtal.operations import angle
except:
    import_fail()

if passed():
    try:
        for i in range(10):
            v1 = random_vector()
            v2 = random_vector()
            angle(v1, v2)
    except:
        fail()

check()

print("  function random_shear_matrix")
try:
    from pyxtal.operations import random_shear_matrix
except:
    import_fail()

if passed():
    try:
        for i in range(10):
            random_shear_matrix()
    except:
        fail()

check()

print("  function is_orthogonal")
try:
    from pyxtal.operations import is_orthogonal
except:
    import_fail()

if passed():
    try:
        a = is_orthogonal([[1,0,0],[0,1,0],[0,0,1]])
        b = is_orthogonal([[0,0,1],[1,0,0],[1,0,0]])
        if a is True and b is False:
            pass
        else:
            fail()
    except:
        fail()

check()

print("  function aa2matrix")
try:
    from pyxtal.operations import aa2matrix
except:
    import_fail()

if passed():
    try:
        for i in range(10):
            aa2matrix(1, 1, random=True)
    except:
        fail()

check()

print("  function matrix2aa")
try:
    from pyxtal.operations import matrix2aa
except:
    import_fail()

if passed():
    try:
        for i in range(10):
            m = aa2matrix(1, 1, random=True)
            aa = matrix2aa(m)
    except:
        fail()

check()

print("  function rotate_vector")
try:
    from pyxtal.operations import rotate_vector
except:
    import_fail()

if passed():
    try:
        for i in range(10):
            v1 = random_vector()
            v2 = random_vector()
            rotate_vector(v1, v2)
    except:
        fail()

check()

print("  function are_equal")
try:
    from pyxtal.operations import are_equal
except:
    import_fail()

if passed():
    try:
        op1 = SymmOp.from_xyz_string('x,y,z')
        op2 = SymmOp.from_xyz_string('x,y,z+1')
        a = are_equal(op1, op2, allow_pbc=True)
        b = are_equal(op1, op2, allow_pbc=False)
        if a is True and b is False:
            pass
        else:
            fail()
    except:
        fail()

check()


print("  class OperationAnalyzer")
try:
    from pyxtal.operations import OperationAnalyzer
except:
    import_fail()

if passed():
    try:
        for i in range(10):
            m = aa2matrix(1,1,random=True)
            t = random_vector()
            op1 = SymmOp.from_rotation_and_translation(m, t)
            OperationAnalyzer(op1)
    except:
        fail()

check()

print("  class orientation")
try:
    from pyxtal.operations import orientation
except:
    import_fail()

if passed():
    try:
        for i in range(10):
            v1 = random_vector()
            c1 = random_vector()
            o = orientation.from_constraint(v1, c1)
    except:
        fail()

check()



end()
