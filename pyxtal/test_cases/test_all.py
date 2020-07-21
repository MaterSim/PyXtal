"""
Test script for pyXtal version 0.1dev. Tests core functions for all modules.
"""

# Custom print function for output to file
_summary_text_ = ""

oldprint = print


def newprint(text):
    global _summary_text_
    oldprint(text)
    if _summary_text_ != "":
        _summary_text_ += "\n"
    _summary_text_ += text


print = newprint

import sys

sys.settrace(None)

outstructs = []
outstrings = []


class Logger(object):
    def __init__(self):
        self.terminal = sys.stdout
        self.log = open("Summary.txt", "w")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)

    def flush(self):
        # this flush method is needed for python 3 compatibility.
        # this handles the flush command by doing nothing.
        # you might want to specify some extra behavior here.
        pass


sys.stdout = Logger()


def compare_wyckoffs(num1, num2, dim=3):
    """Given 2 groups, return whether the second point
    group has equal or greater symmetry than the first group."""
    from numpy import allclose

    if num1 == "???":
        print("Error: invalid value for num1 passed to compare_wyckoffs")
        return
    if num2 == "???":
        return False
    # Get general positions for both groups
    if dim == 3:
        from pyxtal.symmetry import get_wyckoffs

        g1 = get_wyckoffs(num1)[0]
        g2 = get_wyckoffs(num2)[0]
    elif dim == 2:
        from pyxtal.symmetry import get_layer

        g1 = get_layer(num1)[0]
        g2 = get_layer(num2)[0]
    elif dim == 1:
        from pyxtal.symmetry import get_rod

        g1 = get_rod(num1)[0]
        g2 = get_rod(num2)[0]
    elif dim == 0:
        from pyxtal.symmetry import get_point

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
        if not allclose(m1, m2):
            return False
    return True


def check_struct_group(crystal, group, dim=3, tol=1e-2):
    # Supress pymatgen/numpy complex casting warnings
    from pyxtal.crystal import random_crystal
    from pyxtal.molecular_crystal import molecular_crystal
    from copy import deepcopy
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        """Given a pymatgen structure, group number, and dimension, return
        whether or not the structure matches the group number."""
        if isinstance(crystal, (random_crystal, molecular_crystal)):
            lattice = crystal.struct.lattice.matrix
            if dim != 0:
                old_coords = deepcopy(crystal.struct.frac_coords)
                old_species = deepcopy(crystal.struct.atomic_numbers)
            elif dim == 0:
                old_coords = deepcopy(crystal.cart_coords)
                old_species = deepcopy(crystal.species)
        else:
            lattice = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
            old_coords = np.array(crystal)
            old_species = ["C"] * len(old_coords)

        from pyxtal.symmetry import distance
        from pyxtal.symmetry import filtered_coords
        from copy import deepcopy

        PBC = [1, 1, 1]

        # Obtain the generators for the group
        if dim == 3:
            from pyxtal.symmetry import get_wyckoffs

            generators = get_wyckoffs(group)[0]

        elif dim == 2:
            from pyxtal.symmetry import get_layer

            generators = get_layer(group)[0]
            PBC = [1, 1, 0]
        elif dim == 1:
            from pyxtal.symmetry import get_rod

            generators = get_rod(group)[0]
            PBC = [0, 0, 1]
        elif dim == 0:
            from pyxtal.symmetry import Group

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
        print("~~~ Error:")
        import pdb, traceback

        extype, value, tb = sys.exc_info()
        traceback.print_exc()
    except:
        print("~~~ Error: ", e)


# Print whether module passed or failed
def check():
    if passed():
        pass  # print("Success!")
    else:
        print("~~~ Failed module ~~~")


# Call at end of script, or if module fails
def end(condition=1):
    print("===")
    if failed_package is False:
        print("All modules passed!")
        if condition == 1:
            sys.exit(0)
        elif condition == 2:
            pass
    else:
        print("One or more modules failed. Try reinstalling the package.")
        sys.exit(0)


def test_atomic():
    global outstructs
    global outstrings
    print("=== Testing generation of atomic 3D crystals. This may take some time. ===")
    from time import time
    from spglib import get_symmetry_dataset
    from pyxtal.symmetry import get_wyckoffs
    from pyxtal.crystal import random_crystal
    from pyxtal.crystal import cellsize
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

    slow = []
    failed = []
    print("  Spacegroup #  |Generated (SPG)|Generated (PMG)|  Time Elapsed")
    skip = (
        []
    )  # [124, 139, 166, 167, 196, 202, 203, 204, 207, 209, 210, 216, 217, 219, 220, 221, 223, 225, 226, 227, 228, 229, 230] #slow to generate
    for sg in range(1, 231):
        if sg not in skip:
            multiplicity = len(get_wyckoffs(sg)[0]) / cellsize(
                sg
            )  # multiplicity of the general position
            start = time()
            rand_crystal = random_crystal(sg, ["C"], [multiplicity], 1.0)
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
                ans1 = get_symmetry_dataset(rand_crystal.spg_struct, symprec=1e-1)
                if ans1 is None:
                    ans1 = "???"
                else:
                    ans1 = ans1["number"]
                sga = SpacegroupAnalyzer(rand_crystal.struct)
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
                print(
                    "\t"
                    + str(sg)
                    + "\t|\t"
                    + str(ans1)
                    + "\t|\t"
                    + str(ans2)
                    + "\t|\t"
                    + t
                )
            else:
                print(
                    "~~~~ Error: Could not generate space group "
                    + str(sg)
                    + " after "
                    + t
                )
                failed.append(sg)
    if slow != []:
        print("~~~~ The following space groups took more than 60 seconds to generate:")
        for i in slow:
            print("     " + str(i))
    if failed != []:
        print("~~~~ The following space groups failed to generate:")
        for i in failed:
            print("     " + str(i))


def test_molecular():
    global outstructs
    global outstrings
    print(
        "=== Testing generation of molecular 3D crystals. This may take some time. ==="
    )
    from time import time
    from spglib import get_symmetry_dataset
    from pyxtal.symmetry import get_wyckoffs
    from pyxtal.crystal import cellsize
    from pyxtal.molecular_crystal import molecular_crystal
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

    slow = []
    failed = []
    print("  Spacegroup #  |Generated (SPG)|Generated (PMG)|  Time Elapsed")
    skip = (
        []
    )  # [24, 183, 202, 203, 209, 210, 216, 219, 225, 226, 227, 228, 229, 230] #slow
    for sg in range(1, 231):
        if sg not in skip:
            multiplicity = len(get_wyckoffs(sg)[0]) / cellsize(
                sg
            )  # multiplicity of the general position
            start = time()
            rand_crystal = molecular_crystal(sg, ["H2O"], [multiplicity], 2.5)
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
                ans1 = get_symmetry_dataset(rand_crystal.spg_struct, symprec=1e-1)
                if ans1 is None:
                    ans1 = "???"
                else:
                    ans1 = ans1["number"]
                sga = SpacegroupAnalyzer(rand_crystal.struct)
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
                        outstrings.append(str("3D_Molecular_" + str(sg) + ".vasp"))
                print(
                    "\t"
                    + str(sg)
                    + "\t|\t"
                    + str(ans1)
                    + "\t|\t"
                    + str(ans2)
                    + "\t|\t"
                    + t
                )
            else:
                print(
                    "~~~~ Error: Could not generate space group "
                    + str(sg)
                    + " after "
                    + t
                )
                failed.append(sg)
    if slow != []:
        print("~~~~ The following space groups took more than 60 seconds to generate:")
        for i in slow:
            print("     " + str(i))
    if failed != []:
        print("~~~~ The following space groups failed to generate:")
        for i in failed:
            print("     " + str(i))


def test_atomic_2D():
    global outstructs
    global outstrings
    print("=== Testing generation of atomic 2D crystals. This may take some time. ===")
    from time import time
    from pyxtal.symmetry import Group
    from pyxtal.crystal import random_crystal_2D
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

    slow = []
    failed = []
    print("  Layer group # |     Symbol    |  Time Elapsed")
    skip = []
    for sg in range(1, 81):
        if sg not in skip:
            g = Group(sg, dim=2)
            multiplicity = len(g[0])  # multiplicity of the general position
            start = time()
            rand_crystal = random_crystal_2D(sg, ["C"], [multiplicity], 4.0)
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
                    outstructs.append(rand_crystal.struct)
                    outstrings.append(str("atomic_2D_" + str(sg) + ".vasp"))
                symbol = g.symbol
                print("\t" + str(sg) + "\t|\t" + symbol + "\t|\t" + t)
            else:
                print(
                    "~~~~ Error: Could not generate layer group "
                    + str(sg)
                    + " after "
                    + t
                )
                failed.append(sg)
    if slow != []:
        print("~~~~ The following layer groups took more than 60 seconds to generate:")
        for i in slow:
            print("     " + str(i))
    if failed != []:
        print("~~~~ The following layer groups failed to generate:")
        for i in failed:
            print("     " + str(i))


def test_molecular_2D():
    global outstructs
    global outstrings
    print(
        "=== Testing generation of molecular 2D crystals. This may take some time. ==="
    )
    from time import time
    from pyxtal.symmetry import Group
    from pyxtal.molecular_crystal import molecular_crystal_2D
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

    slow = []
    failed = []
    print("  Layer group # |     Symbol    |  Time Elapsed")
    skip = []
    for sg in range(1, 81):
        if sg not in skip:
            g = Group(sg, dim=2)
            multiplicity = len(g[0])  # multiplicity of the general position
            start = time()
            rand_crystal = molecular_crystal_2D(sg, ["H2O"], [multiplicity], 4.0)
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
                    outstructs.append(rand_crystal.struct)
                    outstrings.append(str("molecular_2D_" + str(sg) + ".vasp"))
                symbol = g.symbol
                print("\t" + str(sg) + "\t|\t" + symbol + "\t|\t" + t)
            else:
                print(
                    "~~~~ Error: Could not generate layer group "
                    + str(sg)
                    + " after "
                    + t
                )
                failed.append(sg)
    if slow != []:
        print("~~~~ The following layer groups took more than 60 seconds to generate:")
        for i in slow:
            print("     " + str(i))
    if failed != []:
        print("~~~~ The following layer groups failed to generate:")
        for i in failed:
            print("     " + str(i))


def test_atomic_1D():
    global outstructs
    global outstrings
    print("=== Testing generation of atomic 1D crystals. This may take some time. ===")
    from time import time
    from spglib import get_symmetry_dataset
    from pyxtal.symmetry import get_rod
    from pyxtal.crystal import random_crystal_1D
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

    slow = []
    failed = []
    print("    Rod group   | Gen sg. (SPG) | Gen. sg (PMG) |Time Elapsed")
    skip = []  # slow to generate
    for num in range(1, 76):
        if num not in skip:
            multiplicity = len(get_rod(num)[0])  # multiplicity of the general position
            start = time()
            rand_crystal = random_crystal_1D(num, ["H"], [multiplicity], 4.0)
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
                    ans1 = get_symmetry_dataset(rand_crystal.spg_struct, symprec=1e-1)
                except:
                    ans1 = "???"
                if ans1 is None or ans1 == "???":
                    ans1 = "???"
                else:
                    ans1 = ans1["number"]
                sga = SpacegroupAnalyzer(rand_crystal.struct)
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
                        outstructs.append(rand_crystal.struct)
                        outstrings.append(str("1D_Atomic_" + str(num) + ".vasp"))
                print(
                    "\t"
                    + str(num)
                    + "\t|\t"
                    + str(ans1)
                    + "\t|\t"
                    + str(ans2)
                    + "\t|\t"
                    + t
                )
            else:
                print(
                    "~~~~ Error: Could not generate layer group "
                    + str(num)
                    + " after "
                    + t
                )
                failed.append(num)
    if slow != []:
        print("~~~~ The following layer groups took more than 60 seconds to generate:")
        for i in slow:
            print("     " + str(i))
    if failed != []:
        print("~~~~ The following layer groups failed to generate:")
        for i in failed:
            print("     " + str(i))


def test_molecular_1D():
    global outstructs
    global outstrings
    print(
        "=== Testing generation of molecular 1D crystals. This may take some time. ==="
    )
    from time import time
    from spglib import get_symmetry_dataset
    from pyxtal.symmetry import get_rod
    from pyxtal.molecular_crystal import molecular_crystal_1D
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

    slow = []
    failed = []
    print("    Rod group   | Gen sg. (SPG) | Gen. sg (PMG) |Time Elapsed")
    skip = []  # slow to generate
    for num in range(1, 76):
        if num not in skip:
            multiplicity = len(get_rod(num)[0])  # multiplicity of the general position
            start = time()
            rand_crystal = molecular_crystal_1D(num, ["H2O"], [multiplicity], 4.0)
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
                    ans1 = get_symmetry_dataset(rand_crystal.spg_struct, symprec=1e-1)
                except:
                    ans1 = "???"
                if ans1 is None or ans1 == "???":
                    ans1 = "???"
                else:
                    ans1 = ans1["number"]
                sga = SpacegroupAnalyzer(rand_crystal.struct)
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
                        outstructs.append(rand_crystal.struct)
                        outstrings.append(str("1D_Molecular_" + str(num) + ".vasp"))
                print(
                    "\t"
                    + str(num)
                    + "\t|\t"
                    + str(ans1)
                    + "\t|\t"
                    + str(ans2)
                    + "\t|\t"
                    + t
                )
            else:
                print(
                    "~~~~ Error: Could not generate layer group "
                    + str(num)
                    + " after "
                    + t
                )
                failed.append(num)
    if slow != []:
        print("~~~~ The following layer groups took more than 60 seconds to generate:")
        for i in slow:
            print("     " + str(i))
    if failed != []:
        print("~~~~ The following layer groups failed to generate:")
        for i in failed:
            print("     " + str(i))


def test_cluster():
    global outstructs
    global outstrings
    print(
        "=== Testing generation of point group clusters. This may take some time. ==="
    )
    from time import time
    from spglib import get_symmetry_dataset
    from pyxtal.symmetry import Group
    from pyxtal.crystal import random_cluster
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

    slow = []
    failed = []
    print("  Point group # |     Symbol    |  Time Elapsed")
    skip = []  # [32,55,56]#[28,29,30,31,32,55,56]
    for sg in range(1, 57):
        if sg not in skip:
            multiplicity = len(
                Group(sg, dim=0)[0]
            )  # multiplicity of the general position
            start = time()
            rand_crystal = random_cluster(sg, ["C"], [multiplicity], 1.0)
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
                    outstructs.append(rand_crystal.struct)
                    outstrings.append(str("Cluster_" + str(sg) + ".vasp"))
                pgsymbol = Group(sg, dim=0).symbol
                print("\t" + str(sg) + "\t|\t" + pgsymbol + "\t|\t" + t)
            else:
                print(
                    "~~~~ Error: Could not generate space group "
                    + str(sg)
                    + " after "
                    + t
                )
                failed.append(sg)
    if slow != []:
        print("~~~~ The following space groups took more than 60 seconds to generate:")
        for i in slow:
            print("     " + str(i))
    if failed != []:
        print("~~~~ The following space groups failed to generate:")
        for i in failed:
            print("     " + str(i))


def test_modules():
    print("====== Testing functionality for pyXtal version 0.1dev ======")

    global failed_package
    failed_package = False  # Record if errors occur at any level

    reset()

    print("Importing sys...")
    try:
        import sys

        print("Success!")
    except Exception as e:
        fail(e)
        sys.exit(0)

    print("Importing numpy...")
    try:
        import numpy as np

        print("Success!")
    except Exception as e:
        fail(e)
        sys.exit(0)

    I = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])

    print("Importing pymatgen...")
    try:
        import pymatgen

        print("Success!")
    except Exception as e:
        fail(e)
        sys.exit(0)

    try:
        from pymatgen.core.operations import SymmOp
    except Exception as e:
        fail(e)
        sys.exit(0)

    print("Importing pandas...")
    try:
        import pandas

        print("Success!")
    except Exception as e:
        fail(e)
        sys.exit(0)

    print("Importing spglib...")
    try:
        import spglib

        print("Success!")
    except Exception as e:
        fail(e)
        sys.exit(0)

    print("Importing openbabel...")
    try:
        import ase

        print("Success!")
    except:
        print("Error: could not import openbabel. Try reinstalling the package.")

    print("Importing pyxtal...")
    try:
        import pyxtal

        print("Success!")
    except Exception as e:
        fail(e)
        sys.exit(0)

    print("=== Testing modules ===")

    # =====database.element=====
    print("pyxtal.database.element")
    reset()
    try:
        import pyxtal.database.element
    except Exception as e:
        fail(e)

    print("  class Element")
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
    print("pyxtal.database.hall")
    reset()
    try:
        import pyxtal.database.hall
    except Exception as e:
        fail(e)

    print("  hall_from_hm")
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
    print("pyxtal.database.collection")
    reset()
    try:
        import pyxtal.database.collection
    except Exception as e:
        fail(e)

    print("  Collection")
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
    print("pyxtal.operations")
    reset()
    try:
        import pyxtal.operations
    except Exception as e:
        fail(e)

    print("  random_vector")
    try:
        from pyxtal.operations import random_vector
    except Exception as e:
        fail(e)

    if passed():
        try:
            for i in range(10):
                random_vector()
        except Exception as e:
            fail(e)

    check()

    print("  angle")
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

    print("  random_shear_matrix")
    try:
        from pyxtal.operations import random_shear_matrix
    except Exception as e:
        fail(e)

    if passed():
        try:
            for i in range(10):
                random_shear_matrix()
        except Exception as e:
            fail(e)

    check()

    print("  is_orthogonal")
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

    print("  aa2matrix")
    try:
        from pyxtal.operations import aa2matrix
    except Exception as e:
        fail(e)

    if passed():
        try:
            for i in range(10):
                aa2matrix(1, 1, random=True)
        except Exception as e:
            fail(e)

    check()

    print("  matrix2aa")
    try:
        from pyxtal.operations import matrix2aa
    except Exception as e:
        fail(e)

    if passed():
        try:
            for i in range(10):
                m = aa2matrix(1, 1, random=True)
                aa = matrix2aa(m)
        except Exception as e:
            fail(e)

    check()

    print("  rotate_vector")
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

    print("  are_equal")
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

    print("  class OperationAnalyzer")
    try:
        from pyxtal.operations import OperationAnalyzer
    except Exception as e:
        fail(e)

    if passed():
        try:
            for i in range(10):
                m = aa2matrix(1, 1, random=True)
                t = random_vector()
                op1 = SymmOp.from_rotation_and_translation(m, t)
                OperationAnalyzer(op1)
        except Exception as e:
            fail(e)

    check()

    print("  class Orientation")
    try:
        from pyxtal.operations import Orientation
    except Exception as e:
        fail(e)

    if passed():
        try:
            for i in range(10):
                v1 = random_vector()
                c1 = random_vector()
                o = Orientation.from_constraint(v1, c1)
        except Exception as e:
            fail(e)

    check()

    # =====symmetry=====
    print("pyxtal.symmetry")
    reset()
    try:
        import pyxtal.symmetry
    except Exception as e:
        fail(e)

    print("  get_wyckoffs (may take a moment)")
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

    print("  get_wyckoff_symmetry (may take a moment)")
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

    print("  get_wyckoffs_generators (may take a moment)")
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

    print("  letter_from_index")
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

    print("  index_from_letter")
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

    print("  jk_from_i")
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
                print(j, k)
                fail()
        except Exception as e:
            fail(e)

    check()

    print("  i_from_jk")
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
                print(j, k)
                fail()
        except Exception as e:
            fail(e)

    check()

    print("  ss_string_from_ops")
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

    print("  Wyckoff_position")
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

    print("  Group")
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

    # =====crystal=====
    print("pyxtal.crystal")
    reset()
    try:
        import pyxtal.crystal
    except Exception as e:
        fail(e)

    print("  random_crystal")
    try:
        from pyxtal.crystal import random_crystal
    except Exception as e:
        fail(e)

    if passed():
        try:
            c = random_crystal(1, ["H"], [1], 10.0)
            if c.valid is True:
                pass
            else:
                fail()
        except Exception as e:
            fail(e)

    check()

    print("  random_crystal_2D")
    try:
        from pyxtal.crystal import random_crystal_2D
    except Exception as e:
        fail(e)

    if passed():
        try:
            c = random_crystal_2D(1, ["H"], [1], 10.0)
            if c.valid is True:
                pass
            else:
                fail()
        except Exception as e:
            fail(e)

    check()

    # =====molecule=====
    print("pyxtal.molecule")
    reset()
    try:
        import pyxtal.molecule
    except Exception as e:
        fail(e)

    check()

    print("  Collections")
    try:
        from pyxtal.molecule import mol_from_collection
    except Exception as e:
        fail(e)

    if passed():
        try:
            h2o = mol_from_collection("H2O")
            ch4 = mol_from_collection("CH4")
        except Exception as e:
            fail(e)

    print("  get_inertia_tensor")
    try:
        from pyxtal.molecule import get_inertia_tensor
    except Exception as e:
        fail(e)

    if passed():
        try:
            get_inertia_tensor(h2o)
            get_inertia_tensor(ch4)
        except Exception as e:
            fail(e)

    check()

    print("  get_moment_of_inertia")
    try:
        from pyxtal.molecule import get_moment_of_inertia
    except Exception as e:
        fail(e)

    if passed():
        try:
            v = random_vector()
            get_moment_of_inertia(h2o, v)
            get_moment_of_inertia(ch4, v)
        except Exception as e:
            fail(e)

    check()

    print("  reoriented_molecule")
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

    print("  orientation_in_wyckoff_position")
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

    # =====molecular_crystal=====
    print("pyxtal.molecular_crystal")
    reset()
    try:
        import pyxtal.crystal
    except Exception as e:
        fail(e)

    print("  molecular_crystal")
    try:
        from pyxtal.molecular_crystal import molecular_crystal
    except Exception as e:
        fail(e)

    if passed():
        try:
            c = molecular_crystal(1, ["H2O"], [1], 10.0)
            if c.valid is True:
                pass
            else:
                fail()
        except Exception as e:
            fail(e)

    check()

    print("  molecular_crystal_2D")
    try:
        from pyxtal.molecular_crystal import molecular_crystal_2D
    except Exception as e:
        fail(e)

    if passed():
        try:
            c = molecular_crystal_2D(1, ["H2O"], [1], 10.0)
            if c.valid is True:
                pass
            else:
                fail()
        except Exception as e:
            fail(e)

    check()

    end(condition=2)


from optparse import OptionParser

if __name__ == "__main__":
    import sys
    from time import time

    parser = OptionParser()
    parser.add_option(
        "-m",
        "--module",
        dest="module",
        metavar="module",
        default="all",
        type=str,
        help="modules options: 'all', 'atomic', 'molecular', 'atomic_2D', 'molecular_2D', 'atomic_1D', 'molecular_1D' ",
    )
    parser.add_option(
        "-s",
        "--summary",
        action="store_true",
        dest="output",
        default=False,
        help="output summary.txt file",
    )
    (options, args) = parser.parse_args()

    output = options.output

    try:
        import numpy as np
    except Exception as e:
        fail(e)
        sys.exit(0)
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
            print("please choose the modules from the followings:")
            for module in modules_lib.keys():
                print(module)

    masterstart = time()

    test_modules()

    for module in modules:
        eval(modules_lib[module])

    masterend = time()
    mastertime = np.around((masterend - masterstart), decimals=2)

    print("TEST COMPLETE")
    print("Total time elapsed: " + str(mastertime) + " s")

    if outstructs != []:
        output = True
    if output is True:
        # from pymatgen.io.cif import CifWriter
        from os import mkdir
        from os.path import isdir

        outdir0 = "test_out_"
        i = 1
        while True:
            outdir = outdir0 + str(i)
            if not isdir(outdir):
                mkdir(outdir)
                break
            i += 1
            if i > 100:
                break
        if outstructs != []:
            print("Some generated space groups did not match the expected group.")
            print(
                "POSCAR files for these groups will be output to the directory "
                + outdir
                + ":"
            )
        for struct, string in zip(outstructs, outstrings):
            fpath = outdir + "/" + string
            try:
                struct = struct.get_sorted_structure()
            except:
                pass
            struct.to(filename=fpath, fmt="poscar")
            # CifWriter(struct, symprec=0.1).write_file(filename = fpath)
            print("  " + string)
        # Output summary text file
        txtfile = open(outdir + "/summary.txt", "w")
        txtfile.write(_summary_text_)
