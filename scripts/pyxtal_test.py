#!/usr/bin/env  python
# encoding: utf-8

# Test script for pyXtal v-0.1.4. Tests core functions for all modules.


import sys
import numpy as np
import warnings

from time import time
from copy import deepcopy

from spglib import get_symmetry_dataset
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pyxtal.symmetry import (
    Group,
    get_wyckoffs,
    get_layer,
    get_rod,
    get_point,
)
from pyxtal import pyxtal
from pyxtal.operations import distance, filtered_coords


_summary_text_ = ""


def fprint(text):
    """Custom print function for output to file
    """
    global _summary_text_
    print(text)
    if _summary_text_ != "":
        _summary_text_ += "\n"
    _summary_text_ += text


sys.settrace(None)

outstructs = []
outstrings = []


class Logger(object):
    def __init__(self):
        self.terminal = sys.stdout
        self.log = open("test_summary.txt", "w")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)

    def flush(self):
        # this flush method is needed for python 3 compatibility.
        # this handles the flush command by doing nothing.
        # you might want to specify some extra behavior here.
        pass


sys.stdout = Logger()


# Check if module and classes work correctly
def passed():
    global failed_module
    global failed
    if failed_module is False and failed is False:
        return True
    else:
        return False


# Reset flags for module and class
def reset():
    global failed_module
    global failed
    failed_module = False
    failed = False


# Set flags for package, module, class if error occurs
def fail(*argv):
    if argv != ():
        e = argv[0]
    else:
        e = "Unknown error"
    global failed_package
    global failed_module
    global failed
    failed_package = True
    failed_module = True
    failed = True
    try:
        fprint("~~~ Error:")
        import pdb, traceback

        extype, value, tb = sys.exc_info()
        traceback.print_exc()
    except:
        fprint("~~~ Error: ", e)


# Print whether module passed or failed
def check():
    if passed():
        pass  # fprint("Success!")
    else:
        fprint("~~~ Failed module ~~~")


# Call at end of script, or if module fails
def end(condition=1):
    fprint("===")
    if failed_package is False:
        fprint("All modules passed!")
        if condition == 1:
            sys.exit(0)
        elif condition == 2:
            pass
    else:
        fprint("One or more modules failed. Try reinstalling the package.")
        sys.exit(0)


def compare_wyckoffs(num1, num2, dim=3):
    """Given 2 groups, return whether the second point
    group has equal or greater symmetry than the first group."""

    if num1 == "???":
        fprint("Error: invalid value for num1 passed to compare_wyckoffs")
        return
    if num2 == "???":
        return False
    # Get general positions for both groups
    if dim == 3:

        g1 = get_wyckoffs(num1)[0]
        g2 = get_wyckoffs(num2)[0]
    elif dim == 2:

        g1 = get_layer(num1)[0]
        g2 = get_layer(num2)[0]
    elif dim == 1:

        g1 = get_rod(num1)[0]
        g2 = get_rod(num2)[0]
    elif dim == 0:

        g1 = get_point(num1)[0]
        g2 = get_point(num2)[0]
    # If group 2 has higher symmetry
    if len(g2) > len(g1):
        return True
    # Compare point group operations
    for i, op2 in enumerate(g2):
        op1 = g1[i]
        m1 = op1.rotation_matrix
        m2 = op2.rotation_matrix
        if not np.allclose(m1, m2):
            return False
    return True


def check_struct_group(crystal, group, dim=3, tol=1e-2):
    # Supress pymatgen/numpy complex casting warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        """Given a pymatgen structure, group number, and dimension, return
        whether or not the structure matches the group number."""
        if isinstance(crystal, pyxtal):
            pmg_struc = crystal.to_pymatgen()
            if dim > 0:
                lattice = pmg_struc.lattice.matrix
            else:
                lattice = crystal.lattice.matrix
            if dim != 0:
                old_coords = deepcopy(pmg_struc.frac_coords)
                old_species = deepcopy(pmg_struc.atomic_numbers)
            elif dim == 0:
                old_coords = deepcopy(pmg_struc.cart_coords)
                old_species = deepcopy(pmg_struc.species)
        else:
            lattice = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
            old_coords = np.array(crystal)
            old_species = ["C"] * len(old_coords)

        PBC = [1, 1, 1]

        # Obtain the generators for the group
        if dim == 3:

            generators = get_wyckoffs(group)[0]

        elif dim == 2:

            generators = get_layer(group)[0]
            PBC = [1, 1, 0]
        elif dim == 1:

            generators = get_rod(group)[0]
            PBC = [0, 0, 1]
        elif dim == 0:

            generators = Group(group, dim=0)[0]
            PBC = [0, 0, 0]

        # TODO: Add check for lattice symmetry

        # Apply SymmOps to generate new points
        # old_coords = filtered_coords(struct.frac_coords,PBC=PBC)

        new_coords = []
        new_species = []
        for i, point in enumerate(old_coords):
            for j, op in enumerate(generators):
                if j != 0:
                    new_coords.append(op.operate(point))
                    new_species.append(old_species[i])
        # new_coords = filtered_coords(new_coords,PBC=PBC)

        # Check that all points in new list are still in old

        failed = False
        i_list = list(range(len(new_coords)))
        for i, point1 in enumerate(new_coords):
            found = False
            for j, point2 in enumerate(old_coords):
                if new_species[i] == old_species[j]:
                    difference = filtered_coords(point2 - point1, PBC=PBC)
                    if distance(difference, lattice, PBC=PBC) <= tol:
                        found = True
                        break
            if found is False:
                failed = True
                break

        if failed is False:
            return True
        else:
            return False


def test_atomic():
    global outstructs
    global outstrings
    fprint("=== Testing generation of atomic 3D crystals. This may take some time. ===")
    my_crystal1 = pyxtal()
    my_crystal1.from_random(3, 99, ['Ba','Ti','O'], [1,1,3], 1.0, sites=[["1b"], ["1b"], ["2c", "1b"]])
    print(my_crystal1)
    #my_crystal1.to_file("1.cif")
    
    my_crystal2 = pyxtal()
    my_crystal2.from_random(3, 225, ['C'], [12], 1.0, sites=[["4a", "8c"]])
    print(my_crystal2)
    #my_crystal2.to_file("2.cif")
    
    my_crystal3 = pyxtal()
    my_crystal3.from_random(3, 225, ['C','Si'], [12, 4], 1.0, sites=[["4a", "8c"], None])
    print(my_crystal3)
    #my_crystal3.to_file("3.cif")


    slow = []
    failed = []
    fprint("  Spacegroup #  |Generated (SPG)|Generated (PMG)|  Time Elapsed")
    skip = []
    # skip = (
    #     [124, 139, 166, 167, 196, 202, 203, 204, 207, 209, 210, 216, 217,
    #     219, 220, 221, 223, 225, 226, 227, 228, 229, 230] #slow to generate
    # )
    for sg in range(1, 231):
        if sg not in skip:
            multiplicity = len(get_wyckoffs(sg)[0]) 
            start = time()
            rand_crystal = pyxtal()
            rand_crystal.from_random(3, sg, ["C"], [multiplicity], 1.0)
            end = time()
            timespent = np.around((end - start), decimals=2)
            t = str(timespent)
            if len(t) == 3:
                t += "0"
            t += " s"
            if timespent >= 1.0:
                t += " ~"
            if timespent >= 3.0:
                t += "~"
            if timespent >= 10.0:
                t += "~"
            if timespent >= 60.0:
                t += "~"
                slow.append(sg)
            if rand_crystal.valid:
                check = False
                #ans1 = get_symmetry_dataset(rand_crystal.spg_struct, symprec=1e-1)
                ans1 = get_symmetry_dataset(rand_crystal.to_ase(), symprec=1e-1)
                if ans1 is None:
                    ans1 = "???"
                else:
                    ans1 = ans1["number"]
                sga = SpacegroupAnalyzer(rand_crystal.to_pymatgen())
                ans2 = "???"
                if sga is not None:
                    try:
                        ans2 = sga.get_space_group_number()
                    except:
                        ans2 = "???"
                if ans2 is None:
                    ans2 = "???"

                # Compare expected and detected groups
                if ans1 == "???" and ans2 == "???":
                    check = True
                elif ans1 == "???":
                    if int(ans2) > sg:
                        pass
                elif ans2 == "???":
                    if int(ans1) > sg:
                        pass
                else:
                    if ans1 < sg and ans2 < sg:
                        if compare_wyckoffs(sg, ans1) or compare_wyckoffs(sg, ans2):
                            pass
                        else:
                            check = True

                # output cif files for incorrect space groups
                if check is True:
                    if check_struct_group(rand_crystal, sg, dim=3):
                        pass
                    else:
                        t += " xxxxx"
                        outstructs.append(rand_crystal.struct)
                        outstrings.append(str("3D_Atomic_" + str(sg) + ".vasp"))
                fprint("\t{}\t|\t{}\t|\t{}\t|\t{}".format(sg, ans1, ans2, t))
            else:
                fprint(
                    "~~~~ Error: Could not generate space group {} after {}".format(sg, t)
                )
                failed.append(sg)
    if slow != []:
        fprint("~~~~ The following space groups took more than 60 seconds to generate:")
        for i in slow:
            fprint("     " + str(i))
    if failed != []:
        fprint("~~~~ The following space groups failed to generate:")
        for i in failed:
            fprint("     " + str(i))

def test_molecular():
    global outstructs
    global outstrings
    fprint(
        "=== Testing generation of molecular 3D crystals. This may take some time. ==="
    )

    slow = []
    failed = []
    fprint("  Spacegroup #  |Generated (SPG)|Generated (PMG)|  Time Elapsed")
    skip = [
        225,
        226,
        227,
        228,
    ]  # [24, 183, 202, 203, 209, 210, 216, 219, 225, 226, 227, 228, 229, 230] #slow
    for sg in range(1, 231):
        if sg not in skip:
            multiplicity = len(get_wyckoffs(sg)[0])
            start = time()
            rand_crystal = pyxtal(molecular=True)
            rand_crystal.from_random(3, sg, ["H2O"], [multiplicity], 2.5)
            end = time()
            timespent = np.around((end - start), decimals=2)
            t = str(timespent)
            if len(t) == 3:
                t += "0"
            t += " s"
            if timespent >= 1.0:
                t += " ~"
            if timespent >= 3.0:
                t += "~"
            if timespent >= 10.0:
                t += "~"
            if timespent >= 60.0:
                t += "~"
                slow.append(sg)
            if rand_crystal.valid:
                check = False
                ans1 = get_symmetry_dataset(rand_crystal.to_ase(), symprec=1e-1)
                if ans1 is None:
                    ans1 = "???"
                else:
                    ans1 = ans1["number"]
                sga = SpacegroupAnalyzer(rand_crystal.to_pymatgen())
                ans2 = "???"
                if sga is not None:
                    try:
                        ans2 = sga.get_space_group_number()
                    except:
                        ans2 = "???"
                if ans2 is None:
                    ans2 = "???"

                # Compare expected and detected groups
                if ans1 == "???" and ans2 == "???":
                    check = True
                elif ans1 == "???":
                    if int(ans2) > sg:
                        pass
                elif ans2 == "???":
                    if int(ans1) > sg:
                        pass
                else:
                    if ans1 < sg and ans2 < sg:
                        if compare_wyckoffs(sg, ans1) or compare_wyckoffs(sg, ans2):
                            pass
                        else:
                            check = True

                # output cif files for incorrect space groups
                if check is True:
                    if check_struct_group(rand_crystal, sg, dim=3):
                        pass
                    else:
                        t += " xxxxx"
                        # rand_crystal.to_file("poscar", "1.vasp")
                        # import sys
                        # sys.exit()
                        outstructs.append(rand_crystal.to_pymatgen())
                        outstrings.append(str("3D_Molecular_" + str(sg) + ".vasp"))
                fprint("\t{}\t|\t{}\t|\t{}\t|\t{}".format(sg, ans1, ans2, t))
            else:
                fprint(
                    "~~~~ Error: Could not generate space group {} after {}".format(sg, t)
                )
                failed.append(sg)
    if slow != []:
        fprint("~~~~ The following space groups took more than 60 seconds to generate:")
        for i in slow:
            fprint("     " + str(i))
    if failed != []:
        fprint("~~~~ The following space groups failed to generate:")
        for i in failed:
            fprint("     " + str(i))


def test_atomic_2D():
    global outstructs
    global outstrings
    fprint("=== Testing generation of atomic 2D crystals. This may take some time. ===")

    slow = []
    failed = []
    fprint("  Layer group # |     Symbol    |  Time Elapsed")
    skip = []
    for sg in range(1, 81):
        if sg not in skip:
            g = Group(sg, dim=2)
            multiplicity = len(g[0])  # multiplicity of the general position
            start = time()
            rand_crystal = pyxtal()
            rand_crystal.from_random(2, sg, ["C"], [multiplicity], 4.0)
            end = time()
            timespent = np.around((end - start), decimals=2)
            t = str(timespent)
            if len(t) == 3:
                t += "0"
            t += " s"
            if timespent >= 1.0:
                t += " ~"
            if timespent >= 3.0:
                t += "~"
            if timespent >= 10.0:
                t += "~"
            if timespent >= 60.0:
                t += "~"
                slow.append(sg)
            if rand_crystal.valid:
                if check_struct_group(rand_crystal, sg, dim=2):
                    pass
                else:
                    t += " xxxxx"
                    outstructs.append(rand_crystal.to_pymatgen())
                    outstrings.append(str("atomic_2D_" + str(sg) + ".vasp"))
                symbol = g.symbol
                fprint("\t{}\t|\t{}\t|\t{}".format(sg, symbol, t))
            else:
                fprint(
                    "~~~~ Error: Could not generate layer group {} after {}".format(sg, t)
                )
                failed.append(sg)
    if slow != []:
        fprint("~~~~ The following layer groups took more than 60 seconds to generate:")
        for i in slow:
            fprint("     " + str(i))
    if failed != []:
        fprint("~~~~ The following layer groups failed to generate:")
        for i in failed:
            fprint("     " + str(i))


def test_molecular_2D():
    global outstructs
    global outstrings
    fprint(
        "=== Testing generation of molecular 2D crystals. This may take some time. ==="
    )

    slow = []
    failed = []
    fprint("  Layer group # |     Symbol    |  Time Elapsed")
    skip = []
    for sg in range(1, 81):
        if sg not in skip:
            g = Group(sg, dim=2)
            multiplicity = len(g[0])  # multiplicity of the general position
            start = time()
            rand_crystal = pyxtal(molecular=True)
            rand_crystal.from_random(2, sg, ["H2O"], [multiplicity], 4.0)
            end = time()
            timespent = np.around((end - start), decimals=2)
            t = str(timespent)
            if len(t) == 3:
                t += "0"
            t += " s"
            if timespent >= 1.0:
                t += " ~"
            if timespent >= 3.0:
                t += "~"
            if timespent >= 10.0:
                t += "~"
            if timespent >= 60.0:
                t += "~"
                slow.append(sg)
            if rand_crystal.valid:
                if check_struct_group(rand_crystal, sg, dim=2):
                    pass
                else:
                    t += " xxxxx"
                    outstructs.append(rand_crystal.to_pymatgen())
                    outstrings.append(str("molecular_2D_" + str(sg) + ".vasp"))
                symbol = g.symbol
                fprint("\t{}\t|\t{}\t|\t{}".format(sg, symbol, t))
            else:
                fprint(
                    "~~~~ Error: Could not generate layer group {} after {}".format(sg, t)
                )
                failed.append(sg)
    if slow != []:
        fprint("~~~~ The following layer groups took more than 60 seconds to generate:")
        for i in slow:
            fprint("     " + str(i))
    if failed != []:
        fprint("~~~~ The following layer groups failed to generate:")
        for i in failed:
            fprint("     " + str(i))


def test_atomic_1D():
    global outstructs
    global outstrings
    fprint("=== Testing generation of atomic 1D crystals. This may take some time. ===")

    slow = []
    failed = []
    fprint("    Rod group   | Gen sg. (SPG) | Gen. sg (PMG) |Time Elapsed")
    skip = []  # slow to generate
    for num in range(1, 76):
        if num not in skip:
            multiplicity = len(get_rod(num)[0])  # multiplicity of the general position
            start = time()
            rand_crystal = pyxtal()
            rand_crystal.from_random(1, num, ["H"], [multiplicity], 4.0)
            end = time()
            timespent = np.around((end - start), decimals=2)
            t = str(timespent)
            if len(t) == 3:
                t += "0"
            t += " s"
            if timespent >= 1.0:
                t += " ~"
            if timespent >= 3.0:
                t += "~"
            if timespent >= 10.0:
                t += "~"
            if timespent >= 60.0:
                t += "~"
                slow.append(num)
            if rand_crystal.valid:
                try:
                    ans1 = get_symmetry_dataset(rand_crystal.to_ase(), symprec=1e-1)
                except:
                    ans1 = "???"
                if ans1 is None or ans1 == "???":
                    ans1 = "???"
                else:
                    ans1 = ans1["number"]
                sga = SpacegroupAnalyzer(rand_crystal.to_pymatgen())
                try:
                    ans2 = sga.get_space_group_number()
                except:
                    ans2 = "???"
                if ans2 is None:
                    ans2 = "???"

                check = True

                # output cif files for incorrect space groups
                if check is True:
                    if check_struct_group(rand_crystal, num, dim=1):
                        pass
                    else:
                        t += " xxxxx"
                        outstructs.append(rand_crystal.to_pymatgen)
                        outstrings.append(str("1D_Atomic_" + str(num) + ".vasp"))
                fprint("\t{}\t|\t{}\t|\t{}\t|\t{}".format(num, ans1, ans2, t))
            else:
                fprint(
                    "~~~~ Error: Could not generate layer group {} after {}".format(
                        num, t
                    )
                )
                failed.append(num)
    if slow != []:
        fprint("~~~~ The following layer groups took more than 60 seconds to generate:")
        for i in slow:
            fprint("     " + str(i))
    if failed != []:
        fprint("~~~~ The following layer groups failed to generate:")
        for i in failed:
            fprint("     " + str(i))


def test_molecular_1D():
    global outstructs
    global outstrings
    fprint(
        "=== Testing generation of molecular 1D crystals. This may take some time. ==="
    )

    slow = []
    failed = []
    fprint("    Rod group   | Gen sg. (SPG) | Gen. sg (PMG) |Time Elapsed")
    skip = []  # slow to generate
    for num in range(1, 76):
        if num not in skip:
            multiplicity = len(get_rod(num)[0])  # multiplicity of the general position
            start = time()
            rand_crystal = pyxtal(molecular=True)
            rand_crystal.from_random(1, num, ["H2O"], [multiplicity], 4.0)
            end = time()
            timespent = np.around((end - start), decimals=2)
            t = str(timespent)
            if len(t) == 3:
                t += "0"
            t += " s"
            if timespent >= 1.0:
                t += " ~"
            if timespent >= 3.0:
                t += "~"
            if timespent >= 10.0:
                t += "~"
            if timespent >= 60.0:
                t += "~"
                slow.append(num)
            if rand_crystal.valid:
                try:
                    ans1 = get_symmetry_dataset(rand_crystal.to_ase(), symprec=1e-1)
                except:
                    ans1 = "???"
                if ans1 is None or ans1 == "???":
                    ans1 = "???"
                else:
                    ans1 = ans1["number"]
                sga = SpacegroupAnalyzer(rand_crystal.to_pymatgen())
                try:
                    ans2 = sga.get_space_group_number()
                except:
                    ans2 = "???"
                if ans2 is None:
                    ans2 = "???"

                check = True

                # output cif files for incorrect space groups
                if check is True:
                    if check_struct_group(rand_crystal, num, dim=1):
                        pass
                    else:
                        t += " xxxxx"
                        outstructs.append(rand_crystal.to_pymatgen())
                        outstrings.append(str("1D_Molecular_" + str(num) + ".vasp"))
                fprint("\t{}\t|\t{}\t|\t{}\t|\t{}".format(num, ans1, ans2, t))
            else:
                fprint(
                    "~~~~ Error: Could not generate layer group {} after {}".format(
                        num, t
                    )
                )
                failed.append(num)
    if slow != []:
        fprint("~~~~ The following layer groups took more than 60 seconds to generate:")
        for i in slow:
            fprint("     " + str(i))
    if failed != []:
        fprint("~~~~ The following layer groups failed to generate:")
        for i in failed:
            fprint("     " + str(i))


def test_cluster():
    global outstructs
    global outstrings
    fprint("=== Testing generation of point group clusters. This may take some time. ===")

    slow = []
    failed = []
    fprint("  Point group # |     Symbol    |  Time Elapsed")
    skip = [56]  # [32,55,56]#[28,29,30,31,32,55,56]
    for sg in range(1, 57):
        if sg not in skip:
            multiplicity = len(
                Group(sg, dim=0)[0]
            )  # multiplicity of the general position
            start = time()
            rand_crystal = pyxtal()
            rand_crystal.from_random(0, sg, ["C"], [multiplicity], 1.0)
            end = time()
            timespent = np.around((end - start), decimals=2)
            t = str(timespent)
            if len(t) == 3:
                t += "0"
            t += " s"
            if timespent >= 1.0:
                t += " ~"
            if timespent >= 3.0:
                t += "~"
            if timespent >= 10.0:
                t += "~"
            if timespent >= 60.0:
                t += "~"
                slow.append(sg)
            if rand_crystal.valid:
                if check_struct_group(rand_crystal, sg, dim=0):
                    pass
                else:
                    t += " xxxxx"
                    outstructs.append(rand_crystal.to_pymatgen())
                    outstrings.append(str("Cluster_" + str(sg) + ".vasp"))
                pgsymbol = Group(sg, dim=0).symbol
                fprint("\t{}\t|\t{}\t|\t{}".format(sg, pgsymbol, t))
            else:
                fprint(
                    "~~~~ Error: Could not generate space group {} after {}".format(sg, t)
                )
                failed.append(sg)
    if slow != []:
        fprint("~~~~ The following space groups took more than 60 seconds to generate:")
        for i in slow:
            fprint("     " + str(i))
    if failed != []:
        fprint("~~~~ The following space groups failed to generate:")
        for i in failed:
            fprint("     " + str(i))


def test_modules():
    fprint("====== Testing functionality for pyXtal version 0.1dev ======")

    global failed_package
    failed_package = False  # Record if errors occur at any level

    reset()

    fprint("Importing sys...")
    try:
        import sys

        fprint("Success!")
    except Exception as e:
        fail(e)
        sys.exit(0)

    fprint("Importing numpy...")
    try:
        import numpy as np

        fprint("Success!")
    except Exception as e:
        fail(e)
        sys.exit(0)

    fprint("Importing pymatgen...")
    try:
        import pymatgen

        fprint("Success!")
    except Exception as e:
        fail(e)
        sys.exit(0)

    try:
        from pymatgen.core.operations import SymmOp
    except Exception as e:
        fail(e)
        sys.exit(0)

    fprint("Importing pandas...")
    try:
        import pandas

        fprint("Success!")
    except Exception as e:
        fail(e)
        sys.exit(0)

    fprint("Importing spglib...")
    try:
        import spglib

        fprint("Success!")
    except Exception as e:
        fail(e)
        sys.exit(0)

    fprint("Importing ase...")
    try:
        import ase

        fprint("Success!")
    except:
        fprint("Error: could not import openbabel. Try reinstalling the package.")

    fprint("=== Testing modules ===")

    # =====database.element=====
    fprint("pyxtal.database.element")
    reset()
    try:
        import pyxtal.database.element
    except Exception as e:
        fail(e)

    fprint("  class Element")
    try:
        from pyxtal.database.element import Element
    except Exception as e:
        fail(e)
    if passed():
        for i in range(1, 95):
            if passed():
                try:
                    ele = Element(i)
                except:
                    fail("Could not access Element # " + str(i))
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
                    fail("Could not access attribute for element # " + str(i))
                try:
                    ele.all_z()
                    ele.all_short_names()
                    ele.all_long_names()
                    ele.all_valences()
                    ele.all_valence_electrons()
                    ele.all_covalent_radii()
                    ele.all_vdw_radii()
                except:
                    fail("Could not access class methods")

    check()

    # =====database.hall=====
    fprint("pyxtal.database.hall")
    reset()
    try:
        import pyxtal.database.hall
    except Exception as e:
        fail(e)

    fprint("  hall_from_hm")
    try:
        from pyxtal.database.hall import hall_from_hm
    except Exception as e:
        fail(e)

    if passed():
        for i in range(1, 230):
            if passed():
                try:
                    hall_from_hm(i)
                except:
                    fail("Could not access hm # " + str(i))

    check()

    # =====database.collection=====
    fprint("pyxtal.database.collection")
    reset()
    try:
        import pyxtal.database.collection
    except Exception as e:
        fail(e)

    fprint("  Collection")
    try:
        from pyxtal.database.collection import Collection
    except Exception as e:
        fail(e)

    if passed():
        for i in range(1, 230):
            if passed():
                try:
                    molecule_collection = Collection("molecules")
                except:
                    fail("Could not access hm # " + str(i))

    check()

    # =====operations=====
    fprint("pyxtal.operations")
    reset()
    try:
        import pyxtal.operations
    except Exception as e:
        fail(e)
    from pyxtal.lattice import random_shear_matrix, random_vector

    fprint("  angle")
    try:
        from pyxtal.operations import angle
    except Exception as e:
        fail(e)

    if passed():
        try:
            for i in range(10):
                v1 = random_vector()
                v2 = random_vector()
                angle(v1, v2)
        except Exception as e:
            fail(e)

    check()

    fprint("  is_orthogonal")
    try:
        from pyxtal.operations import is_orthogonal
    except Exception as e:
        fail(e)

    if passed():
        try:
            a = is_orthogonal([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
            b = is_orthogonal([[0, 0, 1], [1, 0, 0], [1, 0, 0]])
            if a is True and b is False:
                pass
            else:
                fail()
        except Exception as e:
            fail(e)

    check()

    fprint("  rotate_vector")
    try:
        from pyxtal.operations import rotate_vector
    except Exception as e:
        fail(e)

    if passed():
        try:
            for i in range(10):
                v1 = random_vector()
                v2 = random_vector()
                rotate_vector(v1, v2)
        except Exception as e:
            fail(e)

    check()

    fprint("  are_equal")
    try:
        from pyxtal.operations import are_equal
    except Exception as e:
        fail(e)

    if passed():
        try:
            op1 = SymmOp.from_xyz_string("x,y,z")
            op2 = SymmOp.from_xyz_string("x,y,z+1")
            a = are_equal(op1, op2, PBC=[0, 0, 1])
            b = are_equal(op1, op2, PBC=[1, 0, 0])
            if a is True and b is False:
                pass
            else:
                fail()
        except Exception as e:
            fail(e)

    check()

    fprint("  class OperationAnalyzer")
    try:
        from pyxtal.operations import OperationAnalyzer
    except Exception as e:
        fail(e)

    if passed():
        try:
            m = np.eye(3)
            t = random_vector()
            op1 = SymmOp.from_rotation_and_translation(m, t)
            OperationAnalyzer(op1)
        except Exception as e:
            fail(e)

    check()

    # =====symmetry=====
    fprint("pyxtal.symmetry")
    reset()
    try:
        import pyxtal.symmetry
    except Exception as e:
        fail(e)

    fprint("  get_wyckoffs (may take a moment)")
    try:
        from pyxtal.symmetry import get_wyckoffs
    except Exception as e:
        fail(e)

    if passed():
        try:
            for i in [1, 2, 229, 230]:
                get_wyckoffs(i)
                get_wyckoffs(i, organized=True)
        except:
            fail(" Could not access Wyckoff positions for space group # " + str(i))

    check()

    fprint("  get_wyckoff_symmetry (may take a moment)")
    try:
        from pyxtal.symmetry import get_wyckoff_symmetry
    except Exception as e:
        fail(e)

    if passed():
        try:
            for i in [1, 2, 229, 230]:
                get_wyckoff_symmetry(i)
                get_wyckoff_symmetry(i, molecular=True)
        except:
            fail("Could not access Wyckoff symmetry for space group # " + str(i))

    check()

    fprint("  get_wyckoffs_generators (may take a moment)")
    try:
        from pyxtal.symmetry import get_wyckoff_generators
    except Exception as e:
        fail(e)

    if passed():
        try:
            for i in [1, 2, 229, 230]:
                get_wyckoff_generators(i)
        except:
            fail("Could not access Wyckoff generators for space group # " + str(i))

    check()

    fprint("  letter_from_index")
    try:
        from pyxtal.symmetry import letter_from_index
    except Exception as e:
        fail(e)

    if passed():
        try:
            if letter_from_index(0, get_wyckoffs(47)) == "A":
                pass
            else:
                fail()
        except Exception as e:
            fail(e)

    check()

    fprint("  index_from_letter")
    try:
        from pyxtal.symmetry import index_from_letter
    except Exception as e:
        fail(e)

    if passed():
        try:
            if index_from_letter("A", get_wyckoffs(47)) == 0:
                pass
            else:
                fail()
        except Exception as e:
            fail(e)

    check()

    fprint("  jk_from_i")
    try:
        from pyxtal.symmetry import jk_from_i
    except Exception as e:
        fail(e)

    if passed():
        try:
            w = get_wyckoffs(2, organized=True)
            j, k = jk_from_i(1, w)
            if j == 1 and k == 0:
                pass
            else:
                fprint(j, k)
                fail()
        except Exception as e:
            fail(e)

    check()

    fprint("  i_from_jk")
    try:
        from pyxtal.symmetry import i_from_jk
    except Exception as e:
        fail(e)

    if passed():
        try:
            w = get_wyckoffs(2, organized=True)
            j, k = jk_from_i(1, w)
            i = i_from_jk(j, k, w)
            if i == 1:
                pass
            else:
                fprint(j, k)
                fail()
        except Exception as e:
            fail(e)

    check()

    fprint("  ss_string_from_ops")
    try:
        from pyxtal.symmetry import ss_string_from_ops
    except Exception as e:
        fail(e)

    if passed():
        try:
            strings = ["1", "4 . .", "2 3 ."]
            for i, sg in enumerate([1, 75, 195]):
                ops = get_wyckoffs(sg)[0]
                ss_string_from_ops(ops, sg, dim=3)
        except Exception as e:
            fail(e)

    check()

    fprint("  Wyckoff_position")
    try:
        from pyxtal.symmetry import Wyckoff_position
    except Exception as e:
        fail(e)

    if passed():
        try:
            wp = Wyckoff_position.from_group_and_index(20, 1)
        except Exception as e:
            fail(e)

    check()

    fprint("  Group")
    try:
        from pyxtal.symmetry import Group
    except Exception as e:
        fail(e)

    if passed():
        try:
            g3 = Group(230)
            g2 = Group(80, dim=2)
            g1 = Group(75, dim=1)
        except Exception as e:
            fail(e)

    check()

    # =====molecule=====
    fprint("pyxtal.molecule")
    reset()
    try:
        from pyxtal.molecule import pyxtal_molecule
    except Exception as e:
        fail(e)

    if passed():
        try:
            h2o = pyxtal_molecule("H2O").mol
            ch4 = pyxtal_molecule("CH4").mol
        except Exception as e:
            fail(e)
    check()

    fprint("  reoriented_molecule")
    try:
        from pyxtal.molecule import reoriented_molecule
    except Exception as e:
        fail(e)

    if passed():
        try:
            reoriented_molecule(h2o)
            reoriented_molecule(ch4)
        except Exception as e:
            fail(e)

    check()

    fprint("  orientation_in_wyckoff_position")
    try:
        from pyxtal.molecule import orientation_in_wyckoff_position
    except Exception as e:
        fail(e)

    if passed():
        try:
            w = get_wyckoffs(20)
            ws = get_wyckoff_symmetry(20, molecular=True)
            wp = Wyckoff_position.from_group_and_index(20, 1)
            orientation_in_wyckoff_position(h2o, wp)
            orientation_in_wyckoff_position(ch4, wp)
        except Exception as e:
            fail(e)

    check()

    end(condition=2)


if __name__ == "__main__":
    from pyxtal import print_logo
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument(
        "-m",
        "--module",
        dest="module",
        metavar="module",
        default="all",
        type=str,
        help="modules options: 'all', 'atomic', 'molecular', \
                               'atomic_2D', 'molecular_2D', 'atomic_1D',\
                               'molecular_1D', 'cluster' ",
    )
    options = parser.parse_args()

    print_logo()

    modules_lib = {
        "atomic": "test_atomic()",
        "molecular": "test_molecular()",
        "atomic_2D": "test_atomic_2D()",
        "molecular_2D": "test_molecular_2D()",
        "atomic_1D": "test_atomic_1D()",
        "molecular_1D": "test_molecular_1D()",
        "cluster": "test_cluster()",
    }
    if options.module == "all":
        modules = modules_lib
    else:
        if options.module in modules_lib.keys():
            modules = [options.module]
        else:
            fprint("please choose the modules from the followings:")
            for module in modules_lib.keys():
                fprint(module)

    masterstart = time()

    test_modules()

    for module in modules:
        eval(modules_lib[module])

    masterend = time()
    mastertime = masterend - masterstart

    fprint("TEST COMPLETE")
    fprint("\nTotal time elapsed: {:.2f}s".format(mastertime))
