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

print("Importing pyxtal...")
try:
    import pyxtal
    print("Success!")
except:
    print("Error: could not import pyxtal. Try reinstalling the package.")
    sys.exit(0)

print("=== Testing modules ===")
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
                print("Error accessing accessing attribute for layer group # "+str(i))

check()

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
        print("Error generating random vector")

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
        print("Error calculating angle")


end()
