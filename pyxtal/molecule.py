"""
Module for handling molecules. Uses the pymatgen.core.structure.Molecule
class as a base. Has a function for reorienting molecules
(reoriented_molecule), and for calculating valid orientations within a Wyckoff
position based on symmetry (orientation_in_wyckoff_position).
"""
#Imports
#------------------------------
#External Libraries
from pymatgen.core.structure import Molecule
from pymatgen.symmetry.analyzer import PointGroupAnalyzer
from pymatgen.symmetry.analyzer import generate_full_symmops

#PyXtal imports
from pyxtal.operations import *
from pyxtal.database.collection import Collection

#Define functions
#------------------------------
def find_ellipsoid(mol):
    printx("Error: bounding ellipsoid calculator not yet implemented.", priority=0)
    pass

def mol_from_file(fname):
    """
    Reads a file into a pymatgen Molecule. Supported formats include xyz, gaussian,
    and pymatgen JSON. Openbabel is optional but adds additional format options.
    
    Args:
        fname: the file path string
    
    Returns:
        a pymatgen Molecule object
    """
    try:
        return Molecule.from_file(fname)
    except:
        printx("Error: could not import file "+str(fname)+" to Molecule.\n"
            +"Default supported formats are xyz, gaussian and pymatgen JSON molecules.\n"
            +"Installing openbabel allows for more extensions.", priority=1)
        return

def mol_from_string(string, fmt):
    """
    Reads a string into a pymatgen Molecule. Uses the pymatgen IMolecule method from_str.
    
    Args:
        string: a string containing the molecular data
        fmt: the conversion format to use
    
    Returns:
        a pymatgen Molecule object
    """
    try:
        return Molecule.from_str(string, fmt)
    except:
        printx("Error: could not convert string '"+str(fmt)+"' to Molecule.\n"
            +"Default supported formats are xyz, gaussian and pymatgen JSON molecules.\n"
            +"Installing openbabel allows for more extensions.", priority=1)
        return

def mol_from_collection(mname):
    """
    Get a molecule from pyxtal.database.collection

    Args:
        mname: the name of the molecule
    
    Returns:
        a pymatgen Molecule object
    """
    return Collection('molecules')[mname]

def get_inertia_tensor(mol):
    """
    Calculate the symmetric inertia tensor for a Molecule object. Used to find
    the principal axes of symmetry.

    Args:
        mol: a Molecule object

    Returns:
        a 3x3 numpy array representing a moment of inertia tensor
    """
    mo = mol.get_centered_molecule()
    # Initialize elements of the inertia tensor
    I11 = I22 = I33 = I12 = I13 = I23 = 0.0
    for i in range(len(mo)):
        x, y, z = mo.cart_coords[i]
        m = mo[i].specie.number
        I11 += m * (y ** 2 + z ** 2)
        I22 += m * (x ** 2 + z ** 2)
        I33 += m * (x ** 2 + y ** 2)
        I12 += -m * x * y
        I13 += -m * x * z
        I23 += -m * y * z
    return np.array([[I11, I12, I13],
                  [I12, I22, I23],
                  [I13, I23, I33]])

def get_moment_of_inertia(mol, axis, scale=1.0):
    """
    Calculate the moment of inertia of a molecule about an axis.

    Args:
        mol: a Molecule object
        axis: a 3d axis (list or array) to compute the moment about
        scale: changes the length scale of the molecule

    Returns:
        a scalar value for the moment of inertia about the axis
    """
    #convert axis to unit vector
    axis = axis / np.linalg.norm(axis)
    moment = 0
    for i, a in enumerate(mol):
        v = a.coords
        moment += (scale * np.linalg.norm(np.cross(axis, v)) ) ** 2
    return moment

def reoriented_molecule(mol, nested=False):
    """
    Reorient a molecule so that its principal axes are aligned with the
    identity matrix.

    Args:
        mol: a Molecule object
        nested: internal variable to keep track of how many times the function
            has been called recursively

    Returns:
        new_mol, P: new_mol is a reoriented copy of the original molecule. P is
            the 3x3 rotation matrix used to obtain it.
    """
    def reorient(mol):
        new_mol = mol.get_centered_molecule()
        A = get_inertia_tensor(new_mol)
        #Store the eigenvectors of the inertia tensor
        P = np.transpose(np.linalg.eigh(A)[1])
        if np.linalg.det(P) < 0:
            P[0] *= -1
        #reorient the molecule
        P = SymmOp.from_rotation_and_translation(P,[0,0,0])
        new_mol.apply_operation(P)
        #Our molecule should never be inverted during reorientation.
        if np.linalg.det(P.rotation_matrix) < 0:
            printx("Error: inverted reorientation applied.\n"
            +"(Within reoriented_molecule)", priority=0)
        return new_mol, P
    #If needed, recursively apply reorientation (due to numerical errors)
    iterations = 1
    max_iterations = 100
    new_mol, P = reorient(mol)
    while iterations < max_iterations:
        is_okay = True
        for i in range(3):
            for j in range(3):
                x = np.linalg.eigh(get_inertia_tensor(new_mol))[1][i][j]
                okay = True
                if i == j:
                    #Check that diagonal elements are 0 or 1
                    if (not np.isclose(x, 0)) and (not np.isclose(x, 1)):
                        okay = False
                else:
                    #Check that off-diagonal elements are 0
                    if (not np.isclose(x, 0)):
                        okay = False
        if okay is False:
            #If matrix is not diagonal with 1's and/or 0's, reorient
            new_mol, Q = reorient(new_mol)
            P = Q*P
            iterations += 1
        elif okay is True:
            break
    if iterations == max_iterations:
        printx("Error: Could not reorient molecule after "+str(max_iterations)+" attempts\n"
            +str(new_mol)+"\n"
            +str(get_inertia_tensor(new_mol)), priority=0)
        return False
    return new_mol, P

def get_symmetry(mol, already_oriented=False):
    """
    Return a molecule's point symmetry.
    Note: for linear molecules, infinitessimal rotations are treated as 6-fold
    rotations, which works for 3d and 2d point groups.

    Args:
        mol: a Molecule object
        already_oriented: whether or not the principle axes of mol are already
            reoriented. Can save time if True, but is not required.

    Returns:
        a list of SymmOp objects which leave the molecule unchanged when applied
    """
    #For single atoms, we cannot represent the point group using a list of operations
    if len(mol) == 1:
        return []
    pga = PointGroupAnalyzer(mol)
    #Handle linear molecules
    if '*' in pga.sch_symbol:
        if already_oriented == False:
            #Reorient the molecule
            oriented_mol, P = reoriented_molecule(mol)
            pga = PointGroupAnalyzer(oriented_mol)
        pg = pga.get_pointgroup()
        symm_m = []
        for op in pg:
            symm_m.append(op)
        #Add 12-fold  and reflections in place of ininitesimal rotation
        for axis in [[1,0,0],[0,1,0],[0,0,1]]:
            op = SymmOp.from_rotation_and_translation(aa2matrix(axis, pi/6), [0,0,0])
            if pga.is_valid_op(op):
                symm_m.append(op)
                #Any molecule with infinitesimal symmetry is linear;
                #Thus, it possess mirror symmetry for any axis perpendicular
                #To the rotational axis. pymatgen does not add this symmetry
                #for all linear molecules - for example, hydrogen
                if axis == [1,0,0]:
                    symm_m.append(SymmOp.from_xyz_string('x,-y,z'))
                    symm_m.append(SymmOp.from_xyz_string('x,y,-z'))
                    r = SymmOp.from_xyz_string('-x,y,-z')
                elif axis == [0,1,0]:
                    symm_m.append(SymmOp.from_xyz_string('-x,y,z'))
                    symm_m.append(SymmOp.from_xyz_string('x,y,-z'))
                    r = SymmOp.from_xyz_string('-x,-y,z')
                elif axis == [0,0,1]:
                    symm_m.append(SymmOp.from_xyz_string('-x,y,z'))
                    symm_m.append(SymmOp.from_xyz_string('x,-y,z'))
                    r = SymmOp.from_xyz_string('x,-y,-z')
                #Generate a full list of SymmOps for the molecule's pointgroup
                symm_m = generate_full_symmops(symm_m, 1e-3)
                break
        #Reorient the SymmOps into mol's original frame
        if not already_oriented:
            new = []
            for op in symm_m:
                new.append(P.inverse*op*P)
            return new
        elif already_oriented:
            return symm_m
    #Handle nonlinear molecules
    else:
        pg = pga.get_pointgroup()
        symm_m = []
        for op in pg:
            symm_m.append(op)
        return symm_m

def orientation_in_wyckoff_position(mol, wyckoff_position, randomize=True,
    exact_orientation=False, already_oriented=False, allow_inversion=False):
    """
    Tests if a molecule meets the symmetry requirements of a Wyckoff position,
    and returns the valid orientations.

    Args:
        mol: a Molecule object. Orientation is arbitrary
        wyckoff_position: a pyxtal.symmetry.Wyckoff_position object
        randomize: whether or not to apply a random rotation consistent with
            the symmetry requirements
        exact_orientation: whether to only check compatibility for the provided
            orientation of the molecule. Used within general case for checking.
            If True, this function only returns True or False
        already_oriented: whether or not to reorient the principle axes
            when calling get_symmetry. Setting to True can remove redundancy,
            but is not necessary
        allow_inversion: whether or not to allow chiral molecules to be
            inverted. Should only be True if the chemical and biological
            properties of the mirror image are known to be suitable for the
            desired application

    Returns:
        a list of operations.Orientation objects which can be applied to the
        molecule while allowing it to satisfy the symmetry requirements of the
        Wyckoff position. If no orientations are found, returns False.
    """
    #For single atoms, there are no constraints
    if len(mol) == 1:
        return [Orientation([[1,0,0],[0,1,0],[0,0,1]], degrees=2)]
    wyckoffs = wyckoff_position.ops
    w_symm = wyckoff_position.symmetry_m
    index = wyckoff_position.index

    #Obtain the Wyckoff symmetry
    symm_w = w_symm[0]
    pga = PointGroupAnalyzer(mol)

    #Check exact orientation
    if exact_orientation is True:
        mo = deepcopy(mol)
        valid = True
        for op in symm_w:
            if not pga.is_valid_op(op):
                valid = False
        if valid is True:
            return True
        elif valid is False:
            return False

    #Obtain molecular symmetry, exact_orientation==False
    symm_m = get_symmetry(mol, already_oriented=already_oriented)
    #Store OperationAnalyzer objects for each molecular SymmOp
    chiral = True
    opa_m = []
    for op_m in symm_m:
        opa = OperationAnalyzer(op_m)
        opa_m.append(opa)
        if opa.type == "rotoinversion":
            chiral = False
        elif opa.type == "inversion":
            chiral = False
    #If molecule is chiral and allow_inversion is False,
    #check if WP breaks symmetry
    if chiral is True:
        if allow_inversion is False:
            for op in wyckoffs:
                if np.linalg.det(op.rotation_matrix) < 0:
                    printx("Warning: cannot place chiral molecule in spagegroup", priority=2)
                    return False
    #Store OperationAnalyzer objects for each Wyckoff symmetry SymmOp
    opa_w = []
    for op_w in symm_w:
        opa_w.append(OperationAnalyzer(op_w))


    #Check for constraints from the Wyckoff symmetry...
    #If we find ANY two constraints (SymmOps with unique axes), the molecule's
    #point group MUST contain SymmOps which can be aligned to these particular
    #constraints. However, there may be multiple compatible orientations of the
    #molecule consistent with these constraints
    constraint1 = None
    constraint2 = None
    for i, op_w in enumerate(symm_w):
        if opa_w[i].axis is not None:
            constraint1 = opa_w[i]
            for j, op_w in enumerate(symm_w):
                if opa_w[j].axis is not None:
                    dot = np.dot(opa_w[i].axis, opa_w[j].axis)
                    if (not np.isclose(dot, 1, rtol=.01)) and (not np.isclose(dot, -1, rtol=.01)):
                        constraint2 = opa_w[j]
                        break
            break
    #Indirectly store the angle between the constraint axes
    if (constraint1 is not None
        and constraint2 is not None):
        dot_w = np.dot(constraint1.axis, constraint2.axis)
    #Generate 1st consistent molecular constraints
    constraints_m = []
    if constraint1 is not None:
        for i, opa1 in enumerate(opa_m):
            if opa1.is_conjugate(constraint1):
                constraints_m.append([opa1, []])
                #Generate 2nd constraint in opposite direction
                extra = deepcopy(opa1)
                extra.axis = [opa1.axis[0]*-1, opa1.axis[1]*-1, opa1.axis[2]*-1]
                constraints_m.append([extra, []])

    #Remove redundancy for the first constraints
    list_i = list(range(len(constraints_m)))
    list_j = list(range(len(constraints_m)))
    copy = deepcopy(constraints_m)
    for i , c1 in enumerate(copy):
        if i in list_i:
            for j , c2 in enumerate(copy):
                if i > j and j in list_j and j in list_i:
                    #Check if axes are colinear
                    if np.isclose(np.dot(c1[0].axis, c2[0].axis), 1, rtol=.01):
                        list_i.remove(j)
                        list_j.remove(j)
                    #Check if axes are symmetrically equivalent
                    else:
                        cond1 = False
                        cond2 = False
                        for opa in opa_m:
                            if opa.type == "rotation":
                                op = opa.op
                                if np.isclose(np.dot(op.operate(c1[0].axis), c2[0].axis), 1, rtol=.05):
                                    cond1 = True
                                    break
                        if cond1 is True: # or cond2 is True:
                            list_i.remove(j)
                            list_j.remove(j)             
    c_m = deepcopy(constraints_m)
    constraints_m = []
    for i in list_i:
        constraints_m.append(c_m[i])

    #Generate 2nd consistent molecular constraints
    valid = list(range(len(constraints_m)))
    if constraint2 is not None:
        for i, c in enumerate(constraints_m):
            opa1 = c[0]
            for j, opa2 in enumerate(opa_m):
                if opa2.is_conjugate(constraint2):
                    dot_m = np.dot(opa1.axis, opa2.axis)
                    #Ensure that the angles are equal
                    if abs(dot_m - dot_w) < .02 or abs(dot_m + dot_w) < .02:
                        constraints_m[i][1].append(opa2)
                        #Generate 2nd constraint in opposite direction
                        extra = deepcopy(opa2)
                        extra.axis = [opa2.axis[0]*-1, opa2.axis[1]*-1, opa2.axis[2]*-1]
                        constraints_m[i][1].append(extra)
            #If no consistent constraints are found, remove first constraint
            if constraints_m[i][1] == []:
                valid.remove(i)
    copy = deepcopy(constraints_m)
    constraints_m = []
    for i in valid:
        constraints_m.append(copy[i])

    #Generate orientations consistent with the possible constraints
    orientations = []
    #Loop over molecular constraint sets
    for c1 in constraints_m:
        v1 = c1[0].axis
        v2 = constraint1.axis
        T = rotate_vector(v1, v2)
        #Loop over second molecular constraints
        for opa in c1[1]:
            phi = angle(constraint1.axis, constraint2.axis)
            phi2 = angle(constraint1.axis, np.dot(T, opa.axis))
            if np.isclose(phi, phi2, rtol=.01):
                r = math.sin(phi)
                c = np.linalg.norm(np.dot(T, opa.axis) - constraint2.axis)
                theta = math.acos(1 - (c**2)/(2*(r**2)))
                R = aa2matrix(constraint1.axis, theta)
                T2 = np.dot(R, T)
                a = angle(np.dot(T2, opa.axis), constraint2.axis)
                if not np.isclose(a, 0, rtol=.01):
                    T2 = np.dot(np.linalg.inv(R), T)
                a = angle(np.dot(T2, opa.axis), constraint2.axis)
                if not np.isclose(a, 0, rtol=.01):
                    printx("Error: Generated incorrect rotation: "+str(theta)+"\n"
                    +"(Within orientation_in_wyckoff_position)", priority=0)
                o = Orientation(T2, degrees=0)
                orientations.append(o)
        #If there is only one constraint
        if c1[1] == []:
            o = Orientation(T, degrees=1, axis=constraint1.axis)
            orientations.append(o)
    #Ensure the identity orientation is checked if no constraints are found
    if constraints_m == []:
        o = Orientation(np.identity(3), degrees=2)
        orientations.append(o)
    
    #Remove redundancy from orientations
    list_i = list(range(len(orientations)))
    list_j = list(range(len(orientations)))
    for i , o1 in enumerate(orientations):
        if i in list_i:
            for j , o2 in enumerate(orientations):
                if i > j and j in list_j and j in list_i:
                    m1 = o1.get_matrix(angle=0)
                    m2 = o2.get_matrix(angle=0)
                    new_op = SymmOp.from_rotation_and_translation(np.dot(m2, np.linalg.inv(m1)), [0,0,0])
                    P = SymmOp.from_rotation_and_translation(np.linalg.inv(m1), [0,0,0])
                    old_op = P*new_op*P.inverse
                    if pga.is_valid_op(old_op):
                        list_i.remove(j)
                        list_j.remove(j)
    copy = deepcopy(orientations)
    orientations = []
    for i in list_i:
        orientations.append(copy[i])

    #Check each of the found orientations for consistency with the Wyckoff pos.
    #If consistent, put into an array of valid orientations
    allowed = []
    for o in orientations:
        if randomize is True:
            op = o.get_op()
        elif randomize is False:
            op = o.get_op(angle=0)
        mo = deepcopy(mol)
        mo.apply_operation(op)
        if orientation_in_wyckoff_position(mo, wyckoff_position, exact_orientation=True, randomize=False, allow_inversion=allow_inversion) is True:
            allowed.append(o)
    #Return the array of allowed orientations. If there are none, return False
    if allowed == []:
        return False
    else:
        return allowed

#Test Functionality
if __name__ == "__main__":
#---------------------------------------------------

    #Testing water
    mol = Collection('molecules')['H2O']
    print("Original molecule:")
    print(mol)
    print()
    #Apply random rotation to avoid lucky results
    R = aa2matrix(1,1,random=True)
    R_op = SymmOp.from_rotation_and_translation(R,[0,0,0])
    mol.apply_operation(R_op)
    print("Rotated molecule:")
    print(mol)
    print()

    #pga = PointGroupAnalyzer(mol)
    #mol = pga.symmetrize_molecule()['sym_mol']    
