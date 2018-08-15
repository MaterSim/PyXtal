"""
Test script for pyXtal version 0.1dev. Tests core functions for all modules.
"""
import sys
sys.settrace(None)

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
def fail(e):
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

#Print whether module passed or failed
def check():
    if passed():
        pass#print("Success!")
    else:
        print("~~~ Failed module ~~~")

#Call at end of script, or if module fails
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
    print("=== Testing generation of atomic 3D crystals. This may take some time. ===")
    from time import time
    from spglib import get_symmetry_dataset
    from pyxtal.crystal import get_wyckoffs
    from pyxtal.crystal import random_crystal
    slow = []
    print("Spacegroup # | Spacegroup Generated | Time Elapsed")
    skip = [202, 216, 225, 226, 227, 228, 229, 230] #slow to generate
    for sg in range(1, 231):
        if sg not in skip:
            multiplicity = len(get_wyckoffs(sg)[0]) #multiplicity of the general position
            start = time()
            rand_crystal = random_crystal(sg, ['C'], [multiplicity], 1.0)
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
                ans = get_symmetry_dataset(rand_crystal.spg_struct, symprec=1e-1)
                if ans is not None:
                    print("\t"+str(sg)+"\t|\t"+str(ans['number'])+"\t|\t"+t)
                else:
                    print("\t"+str(sg)+"\t|\t"+"???"+"\t|\t"+t)
            else:
                print("~~~~ Error: Could not generate space group "+str(sg)+" after "+t)
    if slow != []:
        print("~~~~ The following space groups took more than 60 seconds to generate:")
        for i in slow:
            print("     "+str(i))

def test_molecular():
    print("=== Testing generation of molecular 3D crystals. This may take some time. ===")
    from time import time
    from spglib import get_symmetry_dataset
    from pyxtal.crystal import get_wyckoffs
    from pyxtal.molecular_crystal import molecular_crystal
    slow = []
    print("Spacegroup # | Spacegroup Generated | Time Elapsed")
    skip = [183, 202, 203, 209, 210, 216, 219, 225, 226, 227, 228, 229, 230] #slow
    for sg in range(1, 231):
        if sg not in skip:
            multiplicity = len(get_wyckoffs(sg)[0]) #multiplicity of the general position
            start = time()
            rand_crystal = molecular_crystal(sg, ['H2O'], [multiplicity], 2.5)
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
                ans = get_symmetry_dataset(rand_crystal.spg_struct, symprec=1e-1)
                if ans is not None:
                    print("\t"+str(sg)+"\t|\t"+str(ans['number'])+"\t|\t"+t)
                else:
                    print("\t"+str(sg)+"\t|\t"+"???"+"\t|\t"+t)
            else:
                print("~~~~ Error: Could not generate space group "+str(sg)+" after "+t)
    if slow != []:
        print("~~~~ The following space groups took more than 60 seconds to generate:")
        for i in slow:
            print("     "+str(i))

def test_atomic_2D():
    print("=== Testing generation of atomic 2D crystals. This may take some time. ===")
    from time import time
    from spglib import get_symmetry_dataset
    from pyxtal.crystal import get_wyckoffs
    from pyxtal.crystal import random_crystal_2D
    from pyxtal.database.layergroup import Layergroup
    slow = []
    print("Layergroup | Spacegroup Expected | Spacegroup Generated | Time Elapsed")
    skip = [13, 18, 22, 24, 25, 26, 30, 33, 39, 40, 42, 43, 45, 47, 48, 52, 53, 54, 56, 57, 60, 61, 62, 63, 64, 72, 75, 76, 78, 79, 80] #slow to generate
    for num in range(1, 81):
        if num not in skip:
            sg = Layergroup(num).sgnumber
            multiplicity = len(get_wyckoffs(sg)[0]) #multiplicity of the general position
            start = time()
            rand_crystal = random_crystal_2D(num, ['H'], [multiplicity], 3.0, 1.0)
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
                ans = get_symmetry_dataset(rand_crystal.spg_struct, symprec=1e-1)
                if ans is not None:
                    print("\t"+str(num)+"\t|\t"+str(sg)+"\t|\t"+str(ans['number'])+"\t|\t"+t)
                else:
                    print("\t"+str(num)+"\t|\t"+str(sg)+"\t|\t"+"???"+"\t|\t"+t)
            else:
                print("~~~~ Error: Could not generate layer group "+str(num)+" after "+t)
    if slow != []:
        print("~~~~ The following layer groups took more than 60 seconds to generate:")
        for i in slow:
            print("     "+str(i))

def test_molecular_2D():
    print("=== Testing generation of molecular 2D crystals. This may take some time. ===")
    from time import time
    from spglib import get_symmetry_dataset
    from pyxtal.crystal import get_wyckoffs
    from pyxtal.molecular_crystal import molecular_crystal_2D
    from pyxtal.database.layergroup import Layergroup
    slow = []
    print("Layergroup | Spacegroup Expected | Spacegroup Generated | Time Elapsed")
    skip = [12, 64, 65, 80] #slow to generate
    for num in range(1, 81):
        if num not in skip:
            sg = Layergroup(num).sgnumber
            multiplicity = len(get_wyckoffs(sg)[0]) #multiplicity of the general position
            start = time()
            rand_crystal = molecular_crystal_2D(num, ['H2O'], [multiplicity], 3.0, 1.0)
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
                ans = get_symmetry_dataset(rand_crystal.spg_struct, symprec=1e-1)
                if ans is not None:
                    print("\t"+str(num)+"\t|\t"+str(sg)+"\t|\t"+str(ans['number'])+"\t|\t"+t)
                else:
                    print("\t"+str(num)+"\t|\t"+str(sg)+"\t|\t"+"???"+"\t|\t"+t)
            else:
                print("~~~~ Error: Could not generate layer group "+str(num)+" after "+t)
    if slow != []:
        print("~~~~ The following layer groups took more than 60 seconds to generate:")
        for i in slow:
            print("     "+str(i))

def test_modules():
    print("====== Testing functionality for pyXtal version 0.1dev ======")

    global failed_package
    failed_package = False #Record if errors occur at any level

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

    I = np.array([[1,0,0],[0,1,0],[0,0,1]])

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

    print("Importing pybtex (needed for ase)...")
    try:
        import pybtex
        print("Success!")
    except:
        print("Error: could not import pybtex. Try reinstalling the package.")
        print("PyXtal will still run, but cannot import molecules from ase.")

    print("Importing ase...")
    try:
        import ase
        print("Success!")
    except:
        print("Error: could not import ase. Try reinstalling the package.")
        print("PyXtal will still run, but cannot import molecules from ase.")

    print("Importing pyxtal...")
    try:
        import pyxtal
        print("Success!")
    except Exception as e:
        fail(e)
        sys.exit(0)

    print("=== Testing modules ===")

    #=====database.element=====
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
                    fail("Could not access Element # "+str(i))
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
                    fail("Could not access attribute for element # "+str(i))
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

    #=====database.hall=====
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
                    fail("Could not access hm # "+str(i))

    check()

    #=====database.layergroup=====
    print("pyxtal.database.layergroup")
    reset()
    try:
        import pyxtal.database.layergroup
    except Exception as e:
        fail(e)

    print("  class Layergroup")
    try:
        from pyxtal.database.layergroup import Layergroup
    except Exception as e:
        fail(e)

    if passed():
        for i in range(1, 81):
            if passed():
                try:
                    lgp = Layergroup(i)
                except:
                    fail("Could not access layer group # "+str(i))
                try:
                    lgp.input
                    lgp.lg
                    lgp.symbol
                    lgp.sgnumber
                    lgp.permutation
                except:
                    fail("Could not access attribute for layer group # "+str(i))

    check()

    #=====operations=====
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
            a = is_orthogonal([[1,0,0],[0,1,0],[0,0,1]])
            b = is_orthogonal([[0,0,1],[1,0,0],[1,0,0]])
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
            op1 = SymmOp.from_xyz_string('x,y,z')
            op2 = SymmOp.from_xyz_string('x,y,z+1')
            a = are_equal(op1, op2, allow_pbc=True)
            b = are_equal(op1, op2, allow_pbc=False)
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
                m = aa2matrix(1,1,random=True)
                t = random_vector()
                op1 = SymmOp.from_rotation_and_translation(m, t)
                OperationAnalyzer(op1)
        except Exception as e:
            fail(e)

    check()

    print("  class orientation")
    try:
        from pyxtal.operations import orientation
    except Exception as e:
        fail(e)

    if passed():
        try:
            for i in range(10):
                v1 = random_vector()
                c1 = random_vector()
                o = orientation.from_constraint(v1, c1)
        except Exception as e:
            fail(e)

    check()

    #=====molecule=====
    print("pyxtal.molecule")
    reset()
    try:
        import pyxtal.molecule
    except Exception as e:
        fail(e)

    print("  get_ase_mol")
    try:
        from pyxtal.molecule import get_ase_mol
    except Exception as e:
        fail(e)

    if passed():
        try:
            h2 = get_ase_mol("H2")
            h2o = get_ase_mol("H2O")
            ch4 = get_ase_mol("CH4")
        except Exception as e:
            fail(e)

    check()

    print("  get_inertia_tensor")
    try:
        from pyxtal.molecule import get_inertia_tensor
    except Exception as e:
        fail(e)

    if passed():
        try:
            get_inertia_tensor(h2)
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
            get_moment_of_inertia(h2, v)
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
            reoriented_molecule(h2)
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
            orientation_in_wyckoff_position(h2, 20, 1)
            orientation_in_wyckoff_position(h2o, 20, 1)
            orientation_in_wyckoff_position(ch4, 20, 1)
        except Exception as e:
            fail(e)

    check()

    #=====crystal=====
    print("pyxtal.crystal")
    reset()
    try:
        import pyxtal.crystal
    except Exception as e:
        fail(e)

    print("  get_wyckoffs (may take a moment)")
    try:
        from pyxtal.crystal import get_wyckoffs
    except Exception as e:
        fail(e)

    if passed():
        try:
            for i in [1, 2, 229, 230]:
                get_wyckoffs(i)
                get_wyckoffs(i, organized=True)
        except:
            fail(" Could not access Wyckoff positions for space group # "+str(i))

    check()

    print("  get_wyckoff_symmetry (may take a moment)")
    try:
        from pyxtal.crystal import get_wyckoff_symmetry
    except Exception as e:
        fail(e)

    if passed():
        try:
            for i in [1, 2, 229, 230]:
                get_wyckoff_symmetry(i)
                get_wyckoff_symmetry(i, molecular=True)
        except:
            fail("Could not access Wyckoff symmetry for space group # "+str(i))

    check()

    print("  get_wyckoffs_generators (may take a moment)")
    try:
        from pyxtal.crystal import get_wyckoff_generators
    except Exception as e:
        fail(e)

    if passed():
        try:
            for i in [1, 2, 229, 230]:
                get_wyckoff_generators(i)
        except:
            fail("Could not access Wyckoff generators for space group # "+str(i))

    check()

    print("  letter_from_index")
    try:
        from pyxtal.crystal import letter_from_index
    except Exception as e:
        fail(e)

    if passed():
        try:
            if letter_from_index(0, 47) == "A":
                pass
            else:
                fail()
        except Exception as e:
            fail(e)

    check()

    print("  index_from_letter")
    try:
        from pyxtal.crystal import index_from_letter
    except Exception as e:
        fail(e)

    if passed():
        try:
            if index_from_letter("A", 47) == 0:
                pass
            else:
                fail()
        except Exception as e:
            fail(e)

    check()

    print("  jk_from_i")
    try:
        from pyxtal.crystal import jk_from_i
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
        from pyxtal.crystal import i_from_jk
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
        from pyxtal.crystal import ss_string_from_ops
    except Exception as e:
        fail(e)

    if passed():
        try:
            strings = ['1','4 . .','2 3 .']
            for i, sg in enumerate([1, 75, 195]):
                ops = get_wyckoffs(sg)[0]
                ss_string_from_ops(ops, sg)
        except Exception as e:
            fail(e)

    check()

    print("  random_crystal")
    try:
        from pyxtal.crystal import random_crystal
    except Exception as e:
        fail(e)

    if passed():
        try:
            c = random_crystal(1, ['H'], [1], 10.0)
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
            c = random_crystal_2D(1, ['H'], [1], 1.0, 10.0)
            if c.valid is True:
                pass
            else:
                fail()
        except Exception as e:
            fail(e)

    check()

    #=====molecular_crystal=====
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
            c = molecular_crystal(1, ['H2O'], [1], 10.0)
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
            c = molecular_crystal_2D(1, ['H2O'], [1], 1.0, 10.0)
            if c.valid is True:
                pass
            else:
                fail()
        except Exception as e:
            fail(e)

    check()

    end(condition=2)

if __name__ == "__main__":
    import sys
    from time import time
    try:
        import numpy as np
    except Exception as e:
        fail(e)
        sys.exit(0)

    masterstart = time()

    test_modules()

    test_atomic()

    test_molecular()

    test_atomic_2D()

    test_molecular_2D()

    masterend = time()
    mastertime = np.around((masterend-masterstart), decimals=2)

    print("TEST COMPLETE")
    print("Total time elapsed: "+str(mastertime)+" s")
