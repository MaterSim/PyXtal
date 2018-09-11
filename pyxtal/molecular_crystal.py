"""
Module for generation of random molecular crystals which meet symmetry
constraints. A pymatgen- or spglib-type structure object is created, which can
be saved to a .cif file. Options (preceded by two dashes) are provided for
command-line usage of the module:  

    spacegroup (-s): the international spacegroup number (between 1 and 230)
        to be generated. In the case of a 2D crystal (using option '-d 2'),
        this will instead be the layer group number (between 1 and 80).
        Defaults to 36.  

    molecule (-e): the chemical formula of the molecule to use. For multiple
        molecule types, separate entries with commas. Ex: "C60", "H2O, CH4,
        NH3". Defaults to H2O  

    numMols (-n): the number of molecules in the PRIMITIVE unit cell (For
        P-type spacegroups, this is the same as the number of molecules in the
        conventional unit cell. For A, B, C, and I-centered spacegroups, this
        is half the number of the conventional cell. For F-centered unit cells,
        this is one fourth the number of the conventional cell.). For multiple
        molecule types, separate entries with commas. Ex: "8", "1, 4, 12".
        Defaults to 4  

    factor (-f): the relative volume factor used to generate the unit cell.
        Larger values result in larger cells, with molecules spaced further
        apart. If generation fails after max attempts, consider increasing this
        value. Defaults to 3.0  

    verbosity (-v): the amount of information which should be printed for each
        generated structure. For 0, only prints the requested and generated
        spacegroups. For 1, also prints the Wyckoff positions and time elapsed.
        For 2, also prints the contents of the generated pymatgen structure.
        Defaults to 0  

    attempts (-a): the number of structures to generate. Note: if any of the
        attempts fail, the number of generated structures will be less than this
        value. Structures will be output to separate cif files. Defaults to 1  

    outdir (-o): the file directory where cif files will be output to. Defaults
        to "out"  

    checkatoms (-c): whether or not to check inter-atomic distances at each step
        of generation. When True, produces more accurate results, but requires
        more computation time for larger molecules. When False, produces less
        accurate results and may require a larger volume factor, but does not
        require more computation time for large molecules. Generally, the flag
        should only be set to False for large, approximately spherical molecules
        like C60. Defaults to True  

    allowinversion (-i): whether or not to allow inversion of chiral molecules
        for spacegroups which contain inversional and/or rotoinversional
        symmetry. This should only be True if the chemical and biological
        properties of the mirror molecule are known and suitable for the desired
        application. Defaults to False  

    dimension (-d): 3 for 3D, or 2 for 2D. If 2D, generates a 2D crystal using
        a layer group number instead of a space group number.  

    thickness (-t): The thickness, in Angstroms, to use when generating a
        2D crystal. Note that this will not necessarily be one of the lattice
        vectors, but will represent the perpendicular distance along the non-
        periodic direction. Defaults to 2.0  
"""
from pyxtal.crystal import *
from pyxtal.molecule import *
from pyxtal.operations import *
from pyxtal.database.collection import Collection
from time import time

molecule_collection = Collection('molecules')
max1 = 30 #Attempts for generating lattices
max2 = 30 #Attempts for a given lattice
max3 = 30 #Attempts for a given Wyckoff position

tol_m = 1.0

def estimate_volume_molecular(numMols, boxes, factor=2.0):
    """
    Estimate the volume needed for a molecular crystal conventional unit cell.

    Args:
        numMols: A list with the number of each type of molecule
        boxes: A list of bounding boxes for each molecule. Obtained from get_box
        factor: a factor to multiply the final result by. Used to increase space
        between molecules

    Returns:
        the estimated volume (in cubic Angstroms) needed for the unit cell
    """
    volume = 0
    for numMol, box in zip(numMols, boxes):
        volume += numMol*(box[1]-box[0])*(box[3]-box[2])*(box[5]-box[4])
    return abs(factor*volume)

def get_sg_orientations(mol, sg, allow_inversion=False, PBC=[1,2,3]):
    """
    Calculate the valid orientations for each Molecule and Wyckoff position.
    Returns a list with 3 indices:

        index 1: the Wyckoff position's 1st index (based on multiplicity)  

        index 2: the WP's 2nd index (within the group of equal multiplicity)  

        index 3: the index of the valid orientation for the molecule/WP pair

    For example, self.valid_orientations[i][j] would be a list of valid
    orientations for self.molecules[i], in the Wyckoff position
    self.wyckoffs_organized[i][j]

    Args:
        mol: a pymatgen Molecule object.
        sg: the international spacegroup number
        allow_inversion: whether or not to allow inversion operations for chiral
            molecules
        PBC: a list of periodic axes (1,2,3)->(x,y,z)

    Returns:
        a list of operations orientation objects for each Wyckoff position. 1st
            and 2nd indices correspond to the Wyckoff position
    """
    wyckoffs = get_wyckoffs(sg, PBC=PBC, organized=True)
    w_symm_all = get_wyckoff_symmetry(sg, PBC=PBC, molecular=True)
    valid_orientations = []
    wp_index = -1
    for i, x in enumerate(wyckoffs):
        valid_orientations.append([])
        for j, wp in enumerate(x):
            wp_index += 1
            allowed = orientation_in_wyckoff_position(mol, wyckoffs, w_symm_all, wp_index, already_oriented=True, allow_inversion=allow_inversion)
            if allowed is not False:
                valid_orientations[-1].append(allowed)
            else:
                valid_orientations[-1].append([])
    return valid_orientations

def get_box(mol):
    """
    Given a molecule, find a minimum orthorhombic box containing it.
    Size is calculated using min and max x, y, and z values, plus the padding defined by the vdw radius
    For best results, call oriented_molecule first.
    
    Args:
        mol: a pymatgen Molecule object

    Returns:
        a list [x1,x2,y1,y2,z1,z2] where x1 is the relative displacement in
        the negative x direction, x2 is the displacement in the positive x
        direction, and so on
    """
    minx, miny, minz, maxx, maxy, maxz = 0.,0.,0.,0.,0.,0.
    #for p in mol:
    for p in mol:
        x, y, z = p.coords
        r = Element(p.species_string).vdw_radius
        if x-r < minx: minx = x-r
        if y-r < miny: miny = y-r
        if z-r < minz: minx = z-r
        if x+r > maxx: maxx = x+r
        if y+r > maxy: maxy = y+r
        if z+r > maxz: maxz = z+r
    return [minx,maxx,miny,maxy,minz,maxz]

def check_distance_molecular(coord1, coord2, indices1, index2, lattice, radii, factor = 1.0, PBC=[1,2,3]):
    """
    Check the distances between two set of molecules. The first set is generally
    larger than the second. Distances between coordinates within the first set
    are not checked, and distances between coordinates within the second set are
    not checked. Only distances between points from different sets are checked.

    Args:
        coord1: multiple lists of fractional coordinates e.g. [[[.1,.6,.4],
            [.3,.8,.2]],[[.4,.4,.4],[.3,.3,.3]]]
        coord2: a list of new fractional coordinates e.g. [[.7,.8,.9],
            [.4,.5,.6]]
        indices1: the corresponding molecular indices of coord1, e.g. [1, 3].
            Indices correspond to which value in radii to use
        index2: the molecular index for coord2. Corresponds to which value in
            radii to use
        lattice: matrix describing the unit cell vectors
        radii: a list of radii used to judge whether or not two molecules
            overlap
        d_factor: the tolerance is multiplied by this amount. Larger values
            mean molecules must be farther apart
        PBC: a list of periodic axes (1,2,3)->(x,y,z)

    Returns:
        a bool for whether or not the atoms are sufficiently far enough apart
    """
    #add PBC
    coord2s = []
    matrix = create_matrix(PBC=PBC)
    for coord in coord2:
        for m in matrix:
            coord2s.append(coord+m)
    coord2 = np.array(coord2s)

    coord2 = np.dot(coord2, lattice)
    if len(coord1)>0:
        for coord, index1 in zip(coord1, indices1):
            coord = np.dot(coord, lattice)
            d_min = np.min(cdist(coord, coord2))

            tol = (radii[index1]+radii[index2])

            if d_min < tol:
                return False
        return True
    else:
        return True

def check_wyckoff_position_molecular(points, orientations, wyckoffs, w_symm_all, exact_translation=False, PBC=[1,2,3]):
    """
    Given a list of points, returns the index of the Wyckoff position within
    the spacegroup.

    Args:
        points: a list of 3d fractional coordinates or SymmOps to check
        wyckoffs: a list of unorganized Wyckoff positions obtained from
            get_wyckoffs
        w_symm_all: a list of site symmetry obtained from get_wyckoff_symmetry
        exact_translation: whether we require two SymmOps to have exactly equal
            translational components. If false, translations related by +-1 are
            considered equal
        PBC: a list of periodic axes (1,2,3)->(x,y,z)

    Returns:
        a single index corresponding to the detected Wyckoff position. If no
        valid Wyckoff position is found, returns False
    """
    index = check_wyckoff_position(points, wyckoffs, w_symm_all, exact_translation=exact_translation, PBC=PBC)
    if index is not False:
        generators = wyckoffs[index]
        point = find_generating_point(points, generators, PBC=PBC)
        if point is not None:
            j, k = jk_from_i(index, orientations)
            if orientations[j][k] != []:
                return index
            else:
                return False
        else:
            print("(Inside check_wyckoff_position_molecular)")
            print("Error: Could not generate points from Wyckoff generators")
            print("wp_index: "+str(index))
            print("Coordinates:")
            for c in points:
                print(c)
            print("Generators:")
            for g in generators:
                print(g.as_xyz_string())
    return False

def merge_coordinate_molecular(coor, lattice, wyckoffs, w_symm_all, tol, orientations, PBC=[1,2,3]):
    """
    Given a list of fractional coordinates, merges them within a given
    tolerance, and checks if the merged coordinates satisfy a Wyckoff
    position. Used for merging general Wyckoff positions into special Wyckoff
    positions within the random_crystal (and its derivative) classes.

    Args:
        coor: a list of fractional coordinates
        lattice: a 3x3 matrix representing the unit cell
        wyckoffs: an unorganized list of Wyckoff positions to check
        w_symm_all: A list of Wyckoff site symmetry obtained from
            get_wyckoff_symmetry
        tol: the cutoff distance for merging coordinates
        orientations: a list of valid molecular orientations within the
            space group
        PBC: a list of periodic axes (1,2,3)->(x,y,z)

    Returns:
        coor, index: (coor) is the new list of fractional coordinates after
        merging, and index is a single index of the Wyckoff position within
        the spacegroup. If merging is unsuccesful, or no index is found,
        returns the original coordinates and False
    """
    while True:
        pairs, graph = find_short_dist(coor, lattice, tol, PBC=PBC)
        index = None
        valid = True
        if len(pairs)>0 and valid is True:
            if len(coor) > len(wyckoffs[-1]):
                merged = []
                groups = connected_components(graph)
                for group in groups:
                    merged.append(get_center(coor[group], lattice, PBC=PBC))
                merged = np.array(merged)
                index = check_wyckoff_position_molecular(merged, orientations, wyckoffs, w_symm_all, exact_translation=False, PBC=PBC)
                if index is False:
                    return coor, False
                elif index is None:
                    valid = False
                else:
                    #Check each possible merged Wyckoff position for orientaitons
                    coor = merged

            else:#no way to merge
                return coor, False
        else:
            if index is None:
                index = check_wyckoff_position_molecular(coor, orientations, wyckoffs, w_symm_all, exact_translation=False, PBC=PBC)
            return coor, index

def choose_wyckoff_molecular(wyckoffs, number, orientations):
    """
    Choose a Wyckoff position to fill based on the current number of molecules
    needed to be placed within a unit cell

    Rules:

        1) The new position's multiplicity is equal/less than (number).
        2) We prefer positions with large multiplicity.
        3) The site must admit valid orientations for the desired molecule.

    Args:
        wyckoffs: an unsorted list of Wyckoff positions
        number: the number of molecules still needed in the unit cell
        orientations: the valid orientations for a given molecule. Obtained
            from get_sg_orientations, which is called within molecular_crystal

    Returns:
        a single index for the Wyckoff position. If no position is found,
        returns False
    """
    if np.random.random()>0.5: #choose from high to low
        for j, wyckoff in enumerate(wyckoffs):
            if len(wyckoff[0]) <= number:
                good_wyckoff = []
                for k, w in enumerate(wyckoff):
                    if orientations[j][k] != []:
                        good_wyckoff.append([j,k])
                if len(good_wyckoff) > 0:
                    return choose(good_wyckoff)
        return False
    else:
        good_wyckoff = []
        for j, wyckoff in enumerate(wyckoffs):
            if len(wyckoff[0]) <= number:
                for k, w in enumerate(wyckoff):
                    if orientations[j][k] != []:
                        good_wyckoff.append([j,k])
        if len(good_wyckoff) > 0:
            return choose(good_wyckoff)
        else:
            return False

class mol_site():
    """
    Class for storing molecular Wyckoff positions and orientations within
    the molecular_crystal class. Each mol_site object represenents an
    entire Wyckoff position, not necessarily a single molecule.
    """
    def __init__(self, mol, position, wp, wp_generators, lattice, PBC=[1,2,3]):
        self.mol = mol
        """A Pymatgen molecule object"""
        self.position = position
        """Relative coordinates of the molecule's center within the unit cell"""
        self.wp = wp
        """The Wyckoff position for the site"""
        self.wp = wp_generators
        """The Wyckoff pgenerators for the site"""
        self.multiplicity = len(wp)
        """The multiplicity of the molecule's Wyckoff position"""
        self.PBC = PBC
        """The periodic boundary condition direction"""

class molecular_crystal():
    """
    Class for storing and generating molecular crystals based on symmetry
    constraints. Based on the crystal.random_crystal class for atomic crystals.
    Given a spacegroup, list of molecule objects, molecular stoichiometry, and
    a volume factor, generates a molecular crystal consistent with the given
    constraints. This crystal is stored as a pymatgen struct via self.struct
    
    Args:
        sg: The international spacegroup number
        molecules: a list of pymatgen.core.structure.Molecule objects for
            each type of molecule
        numMols: A list of the number of each type of molecule within the
            primitive cell (NOT the conventioal cell)
        volume_factor: A volume factor used to generate a larger or smaller
            unit cell. Increasing this gives extra space between molecules
        allow_inversion: Whether or not to allow chiral molecules to be
            inverted. If True, the final crystal may contain mirror images of
            the original molecule. Unless the chemical properties of the mirror
            image are known, it is highly recommended to keep this value False
        orientations: Once a crystal with the same spacegroup and molecular
            stoichiometry has been generated, you may pass its
            valid_orientations attribute here to avoid repeating the
            calculation, but this is not required
        check_atomic_distances: If True, checks the inter-atomic distances
            after each Wyckoff position is added. This requires slightly more
            time, but vastly improves accuracy. For approximately spherical
            molecules, or for large inter-molecular distances, this may be
            turned off
    """
    def __init__(self, sg, molecules, numMols, volume_factor, allow_inversion=False, orientations=None, check_atomic_distances=True):
        
        #Necessary input
        self.Msgs()
        """A list of warning messages to use during generation."""
        numMols = np.array(numMols) #must convert it to np.array
        self.factor = volume_factor
        """The supplied volume factor for the unit cell."""
        self.numMols0 = numMols
        self.sg = sg
        """The international spacegroup number of the crystal."""
        #Reorient the molecules along their principle axes
        oriented_molecules = []
        #Allow support for generating molecules from text via openbable
        for i, mol in enumerate(molecules):
            if type(mol) == str:
                #Read strings into molecules, try collection first,
                #If string not in collection, use SMILES format
                try:
                    mo = molecule_collection[mol]
                except:
                    mo = ob_mol_from_string(mol)
                    mo = pmg_from_ob(mo)
                molecules[i] = mo
        for mol in molecules:
            pga = PointGroupAnalyzer(mol)
            mo = pga.symmetrize_molecule()['sym_mol']
            oriented_molecules.append(mo)
        self.molecules = oriented_molecules
        """A list of pymatgen.core.structure.Molecule objects, symmetrized and
        oriented along their symmetry axes."""
        self.boxes = []
        """A list of bounding boxes for each molecule. Used for estimating
        volume of the unit cell."""
        self.radii = []
        """A list of approximated radii for each molecule type. Used for
        checking inter-molecular distances."""
        for mol in self.molecules:
            self.boxes.append(get_box(reoriented_molecule(mol)[0]))
            max_r = 0
            for site in mol:
                radius = math.sqrt( site.x**2 + site.y**2 + site.z**2 )
                if radius > max_r: max_r = radius
            self.radii.append(max_r+1.0)
        self.numMols = numMols * cellsize(self.sg)
        """The number of each type of molecule in the CONVENTIONAL cell"""
        self.volume = estimate_volume_molecular(self.numMols, self.boxes, self.factor)
        """The volume of the generated unit cell"""
        self.wyckoffs = get_wyckoffs(self.sg)
        """The Wyckoff positions for the crystal's spacegroup."""
        self.wyckoffs_organized = get_wyckoffs(self.sg, organized=True)
        """The Wyckoff positions for the crystal's spacegroup. Sorted by
        multiplicity."""
        self.w_symm = get_wyckoff_symmetry(self.sg, molecular=True)
        """A list of site symmetry operations for the Wyckoff positions, obtained
            from get_wyckoff_symmetry."""
        self.wyckoff_generators = get_wyckoff_generators(self.sg)
        self.check_atomic_distances = check_atomic_distances
        """Whether or not inter-atomic distances are checked at each step."""
        self.allow_inversion = allow_inversion
        """Whether or not to allow chiral molecules to be inverted."""
        #When generating multiple crystals of the same stoichiometry and sg,
        #allow the user to re-use the allowed orientations, to reduce time cost
        if orientations is None:
            self.get_orientations()
        else:
            self.valid_orientations = orientations
            """The valid orientations for each molecule and Wyckoff position.
            May be copied when generating a new molecular_crystal to save a
            small amount of time"""
        self.generate_crystal()


    def Msgs(self):
        self.Msg1 = 'Error: the stoichiometry is incompatible with the wyckoff sites choice'
        self.Msg2 = 'Error: failed in the cycle of generating structures'
        self.Msg3 = 'Warning: failed in the cycle of adding species'
        self.Msg4 = 'Warning: failed in the cycle of choosing wyckoff sites'
        self.Msg5 = 'Finishing: added the specie'
        self.Msg6 = 'Finishing: added the whole structure'

    def get_orientations(self):
        """
        Calculates the valid orientations for each Molecule and Wyckoff
        position. Returns a list with 4 indices:

        index 1: the molecular prototype's index within self.molecules

        index 2: the Wyckoff position's 1st index (based on multiplicity)

        index 3: the WP's 2nd index (within the group of equal multiplicity)

        index 4: the index of the valid orientation for the molecule/WP pair

        For example, self.valid_orientations[i][j][k] would be a list of valid
        orientations for self.molecules[i], in the Wyckoff position
        self.wyckoffs_organized[j][k]
        """
        self.valid_orientations = []
        for mol in self.molecules:
            self.valid_orientations.append([])
            wp_index = -1
            for i, x in enumerate(self.wyckoffs_organized):
                self.valid_orientations[-1].append([])
                for j, wp in enumerate(x):
                    wp_index += 1
                    allowed = orientation_in_wyckoff_position(mol, self.wyckoffs, self.w_symm, wp_index, already_oriented=True, allow_inversion=self.allow_inversion)
                    if allowed is not False:
                        self.valid_orientations[-1][-1].append(allowed)
                    else:
                        self.valid_orientations[-1][-1].append([])

    def check_compatible(self):
        """
        Checks if the number of molecules is compatible with the Wyckoff
        positions. Considers the number of degrees of freedom for each Wyckoff
        position, and makes sure at least one valid combination of WP's exists.
        """
        N_site = [len(x[0]) for x in self.wyckoffs_organized]
        has_freedom = False
        #remove WP's with no freedom once they are filled
        removed_wyckoffs = []
        for i, numMol in enumerate(self.numMols):
            #Check that the number of molecules is a multiple of the smallest Wyckoff position
            if numMol % N_site[-1] > 0:
                return False
            else:
                #Check if smallest WP has at least one degree of freedom
                op = self.wyckoffs_organized[-1][-1][0]
                if op.rotation_matrix.all() != 0.0:
                    if self.valid_orientations[i][-1][-1] != []:
                        has_freedom = True
                else:
                    #Subtract from the number of ions beginning with the smallest Wyckoff positions
                    remaining = numMol
                    for j, x in enumerate(self.wyckoffs_organized):
                        for k, wp in enumerate(x):
                            while remaining >= len(wp) and wp not in removed_wyckoffs:
                                if self.valid_orientations[i][j][k] != []:
                                    #Check if WP has at least one degree of freedom
                                    op = wp[0]
                                    remaining -= len(wp)
                                    if np.allclose(op.rotation_matrix, np.zeros([3,3])):
                                        if (len(self.valid_orientations[i][j][k]) > 1 or
                                            self.valid_orientations[i][j][k][0].degrees > 0):
                                            #NOTE: degrees of freedom may be inaccurate for linear molecules
                                            has_freedom = True
                                        else:
                                            removed_wyckoffs.append(wp)
                                    else:
                                        has_freedom = True
                                else:
                                    removed_wyckoffs.append(wp)
                    if remaining != 0:
                        return False
        if has_freedom:
            return True
        else:
            #Wyckoff Positions have no degrees of freedom
            return 0

        return True

    def generate_crystal(self, max1=max1, max2=max2, max3=max3):
        """
        The main code to generate a random molecular crystal. If successful,
        stores a pymatgen.core.structure object in self.struct and sets
        self.valid to True. If unsuccessful, sets self.valid to False and
        outputs an error message.

        Args:
            max1: the number of attempts for generating a lattice
            max2: the number of attempts for a given lattice
            max3: the number of attempts for a given Wyckoff position
        """
        #Check the minimum number of degrees of freedom within the Wyckoff positions
        degrees = self.check_compatible()
        if degrees is False:
            print(self.Msg1)
            self.struct = None
            self.valid = False
            return
        else:
            if degrees == 0:
                max1 = 10
                max2 = 10
                max3 = 10
            #Calculate a minimum vector length for generating a lattice
            #minvector = max(radius*2 for radius in self.radii)
            all_lengths = []
            for box in self.boxes:
                lengths = []
                for a, b in [[0,1],[2,3],[4,5]]:
                    lengths.append(abs(box[b]-box[a]))
                all_lengths.append(min(lengths))
            minvector = max(all_lengths)
            for cycle1 in range(max1):
                #1, Generate a lattice
                cell_para = generate_lattice(self.sg, self.volume, minvec=minvector)
                if cell_para is None:
                    break
                else:
                    cell_matrix = para2matrix(cell_para)
                    if abs(self.volume - np.linalg.det(cell_matrix)) > 1.0: 
                        print('Error, volume is not equal to the estimated value: ', self.volume, ' -> ', np.linalg.det(cell_matrix))
                        print('cell_para:  ', cell_para)
                        sys.exit(0)

                    molecular_coordinates_total = [] #to store the added molecular coordinates
                    molecular_sites_total = []      #to store the corresponding molecular specie
                    atomic_coordinates_total = [] #to store the added atomic coordinates
                    atomic_sites_total = []      #to store the corresponding atomic specie
                    wps_total = []      #to store corresponding Wyckoff position indices
                    points_total = []   #to store the generating x,y,z points
                    mol_generators_total = []
                    good_structure = False

                    for cycle2 in range(max2):
                        molecular_coordinates_tmp = deepcopy(molecular_coordinates_total)
                        molecular_sites_tmp = deepcopy(molecular_sites_total)
                        atomic_coordinates_tmp = deepcopy(atomic_coordinates_total)
                        atomic_sites_tmp = deepcopy(atomic_sites_total)
                        wps_tmp = deepcopy(wps_total)
                        points_tmp = deepcopy(points_total)
                        mol_generators_tmp = []
                        
                	    #Add molecules specie by specie
                        for numMol, mol in zip(self.numMols, self.molecules):
                            i = self.molecules.index(mol)
                            numMol_added = 0

                            #Now we start to add the specie to the wyckoff position
                            for cycle3 in range(max3):
                                #Choose a random Wyckoff position for given multiplicity: 2a, 2b, 2c
                                #NOTE: The molecular version return wyckoff indices, not ops
                                indices = choose_wyckoff_molecular(self.wyckoffs_organized, numMol-numMol_added, self.valid_orientations[i])
                                if indices is not False:
                                    j, k = indices
                                    if self.valid_orientations[i][j][k] == []:
                                        print("Error: Failed to catch empty set...")
                                        print(i,j,k)
                	    	    #Generate a list of coords from ops
                                    ops = self.wyckoffs_organized[j][k]
                                    point = np.random.random(3)
                                    coords = np.array([op.operate(point) for op in ops])
                                    #merge_coordinate if the atoms are close
                                    if self.check_atomic_distances is False:
                                        mtol = self.radii[i]*2
                                    elif self.check_atomic_distances is True:
                                        mtol = 3.0
                                    coords_toadd, good_merge = merge_coordinate_molecular(coords, cell_matrix, 
                                            self.wyckoffs, self.w_symm, mtol, self.valid_orientations[i])
                                    if good_merge is not False:
                                        wp_index = good_merge
                                        j, k = jk_from_i(wp_index, self.wyckoffs_organized)
                                        wyckoffs = self.wyckoffs_organized[j][k]

                                        #Check that coords_toadd are generated by point
                                        generators = get_wyckoffs(self.sg)[wp_index]
                                        point = find_generating_point(coords_toadd, generators)
                                        if point is None:
                                            print("Error: Could not generate merged coordinates from Wyckoff generators")
                                            print("Space group: "+str(self.sg))
                                            print("wp_index: "+str(wp_index))
                                            print("Coordinates:")
                                            for c in coords_toadd:
                                                print(c)
                                            print("Generators:")
                                            for g in generators:
                                                print(g.as_xyz_string())
                                            self.valid = False
                                            return
                                        coords_toadd = filtered_coords(coords_toadd) #scale the coordinates to [0,1], very important!

                                        #Check inter-molecular distances
                                        if self.check_atomic_distances is False:
                                            if check_distance_molecular(molecular_coordinates_tmp, coords_toadd, molecular_sites_tmp, i, cell_matrix, self.radii):
                                                molecular_coordinates_tmp.append(coords_toadd)
                                                molecular_sites_tmp.append(i)
                                                wps_tmp.append(wp_index)
                                                points_tmp.append(point)
                                                numMol_added += len(coords_toadd)
                                                if numMol_added == numMol:
                                                    molecular_coordinates_total = deepcopy(molecular_coordinates_tmp)
                                                    molecular_sites_total = deepcopy(molecular_sites_tmp)
                                                    wps_total = deepcopy(wps_tmp)
                                                    points_total = deepcopy(points_tmp)
                                                    break

                                        #Check inter-atomic distances
                                        elif self.check_atomic_distances is True:
                                            #Generate atomic coordinates from molecules
                                            mo = deepcopy(self.molecules[i])
                                            #get j, k from wp_index
                                            num = 0
                                            found = False
                                            j, k = jk_from_i(wp_index, self.wyckoffs_organized)
                                            op1 = choose(self.valid_orientations[i][j][k]).get_op()
                                            mo.apply_operation(op1)
                                            ms0 = mol_site(mo, point, self.wyckoffs[wp_index], self.wyckoff_generators[wp_index], cell_matrix)
                                            wp_atomic_sites = [] #The species for the Wyckoff position
                                            wp_atomic_coords = [] #The coords for the Wyckoff position
                                            flag1 = True
                                            for point_index, op2 in enumerate(self.wyckoff_generators[wp_index]):
                                                current_atomic_sites = []
                                                current_atomic_coords = []
                                                for site in mo:
                                                    #Place molecular coordinates in relative coordinates
                                                    #relative_coords = np.dot(np.linalg.inv(np.transpose(cell_matrix)), site.coords)
                                                    #relative_coords = np.dot(np.linalg.inv(cell_matrix), site.coords)
                                                    relative_coords = np.dot(site.coords, np.linalg.inv(cell_matrix))
                                                    center1 = op2.operate(point)
                                                    rot = SymmOp.from_rotation_and_translation(op2.rotation_matrix,[0,0,0])
                                                    relative_coords = rot.operate(relative_coords)
                                                    new_vector = center1 + relative_coords
                                                    new_vector -= np.floor(new_vector)
                                                    current_atomic_sites.append(site.specie.name)
                                                    current_atomic_coords.append(new_vector)
                                                wp_atomic_sites.append(current_atomic_sites)
                                                wp_atomic_coords.append(current_atomic_coords)
                                                #Check distances between molecules in current WP
                                                if point_index == 1:
                                                    for a_index, specie2 in enumerate(current_atomic_sites):
                                                        flag1 = check_distance([wp_atomic_coords[0]], [current_atomic_coords[a_index]], wp_atomic_sites[0], specie2, cell_matrix, d_factor=2.0)
                                                        if flag1 is False:
                                                            break
                                                    if flag1 is False:
                                                        break
                                            if flag1 is True:
                                                #Check distances between current and previous molecular atoms
                                                a = []
                                                for x in wp_atomic_coords:
                                                    a += x
                                                b = []
                                                for x in wp_atomic_sites:
                                                    b += x
                                                flag2 = True
                                                for a_index, specie2 in enumerate(b):
                                                    flag2 = check_distance([atomic_coordinates_total], [a[a_index]], atomic_sites_total, specie2, cell_matrix, d_factor=2.0)
                                                    if flag2 is False: break
                                                if flag2 is True:
                                                    mol_generators_tmp.append(ms0)
                                                    molecular_coordinates_tmp.append(coords_toadd)
                                                    molecular_sites_tmp.append(i)
                                                    atomic_coordinates_tmp += a
                                                    atomic_sites_tmp += b
                                                    wps_tmp.append(wp_index)
                                                    points_tmp.append(point)
                                                    numMol_added += len(coords_toadd)
                                                    if numMol_added == numMol:
                                                        mol_generators_total = deepcopy(mol_generators_tmp)
                                                        molecular_coordinates_total = deepcopy(molecular_coordinates_tmp)
                                                        atomic_sites_total = deepcopy(atomic_sites_tmp)
                                                        atomic_coordinates_total = deepcopy(atomic_coordinates_tmp)
                                                        molecular_sites_total = deepcopy(molecular_sites_tmp)
                                                        wps_total = deepcopy(wps_tmp)
                                                        points_total = deepcopy(points_tmp)
                                                        break

                            if numMol_added != numMol:
                                break  #need to repeat from the 1st species

                        if numMol_added == numMol:
                            #print(self.Msg6)
                            good_structure = True
                            break
                        else: #reset the coordinates and sites
                            molecular_coordinates_total = []
                            molecular_sites_total = []
                            wps_total = []
                    #placing molecules here
                    if good_structure:
                        final_lattice = cell_matrix 
                        final_coor = []
                        final_site = []
                        final_number = []
                        self.mol_generators = []
                        """A list of mol_site objects which can be used to regenerate the crystal."""

                        if self.check_atomic_distances is False:
                            for center0, i, wp_index in zip(points_total, molecular_sites_total, wps_total):
                                mo = deepcopy(self.molecules[i])
                                #get j, k from wp_index
                                num = 0
                                found = False
                                j, k = jk_from_i(wp_index, self.wyckoffs_organized)
                                op1 = choose(self.valid_orientations[i][j][k]).get_op()
                                mo.apply_operation(op1)
                                ms0 = mol_site(mo, center0, self.wyckoffs[wp_index], self.wyckoff_generators[wp_index], cell_matrix)
                                self.mol_generators.append(ms0)
                                for index, op2 in enumerate(self.wyckoff_generators[wp_index]):
                                    for site in mo:
                                        #Place molecular coordinates in relative coordinates
                                        #relative_coords = np.dot(np.linalg.inv(np.transpose(cell_matrix)), site.coords)
                                        #relative_coords = np.dot(np.linalg.inv(cell_matrix), site.coords)
                                        relative_coords = np.dot(site.coords, np.linalg.inv(cell_matrix))
                                        center1 = op2.operate(center0)
                                        rot = SymmOp.from_rotation_and_translation(op2.rotation_matrix,[0,0,0])
                                        relative_coords = rot.operate(relative_coords)
                                        new_vector = center1 + relative_coords
                                        new_vector -= np.floor(new_vector)
                                        final_coor.append(new_vector)
                                        final_site.append(site.specie.name)
                                        final_number.append(site.specie.number)

                        elif self.check_atomic_distances is True:
                            final_coor = deepcopy(atomic_coordinates_total)
                            final_site = deepcopy(atomic_sites_total)
                            final_number = list(Element(ele).z for ele in atomic_sites_total)
                            self.mol_generators = deepcopy(mol_generators_total)
                            """A list of mol_site objects which can be used
                            for generating the crystal."""

                        final_coor -= np.floor(final_coor)
                        if verify_distances(final_coor, final_site, final_lattice, factor=1.0) is True:
                            self.lattice = final_lattice
                            """A 3x3 matrix representing the lattice of the
                            unit cell."""  
                            self.coordinates = np.array(final_coor)
                            """The fractional coordinates for each molecule
                            in the final structure"""
                            self.sites = final_site
                            """The indices within self.molecules corresponding
                            to the type of molecule for each site in
                            self.coordinates."""              
                            self.struct = Structure(final_lattice, final_site, np.array(final_coor))
                            """A pymatgen.core.structure.Structure object for
                            the final generated crystal."""
                            self.spg_struct = (final_lattice, np.array(final_coor), final_number)
                            """A list of information describing the generated
                            crystal, which may be used by spglib for symmetry
                            analysis."""
                            self.valid = True
                            """Whether or not a valid crystal was generated."""
                            return
                        #else: print("Failed final distance check.")
        print("Couldn't generate crystal after max attempts.")
        if degrees == 0:
            print("Note: Wyckoff positions have no degrees of freedom.")
        self.struct = self.Msg2
        self.valid = False
        return self.Msg2

class molecular_crystal_2D():
    """
    A 2d counterpart to molecular_crystal. Given a layer group, list of
    molecule objects, molecular stoichiometry, and
    a volume factor, generates a molecular crystal consistent with the given
    constraints. This crystal is stored as a pymatgen struct via self.struct
    
    Args:
        number: the layer group number between 1 and 80. NOT equal to the
            international space group number, which is between 1 and 230
        molecules: a list of pymatgen.core.structure.Molecule objects for
            each type of molecule
        numMols: A list of the number of each type of molecule within the
            primitive cell (NOT the conventioal cell)
        thickness: the thickness, in Angstroms, of the unit cell in the 3rd
            dimension (the direction which is not repeated periodically)
        volume_factor: A volume factor used to generate a larger or smaller
            unit cell. Increasing this gives extra space between molecules
        allow_inversion: Whether or not to allow chiral molecules to be
            inverted. If True, the final crystal may contain mirror images of
            the original molecule. Unless the chemical properties of the mirror
            image are known, it is highly recommended to keep this value False
        orientations: Once a crystal with the same spacegroup and molecular
            stoichiometry has been generated, you may pass its
            valid_orientations attribute here to avoid repeating the
            calculation, but this is not required
        check_atomic_distances: If True, checks the inter-atomic distances
            after each Wyckoff position is added. This requires slightly more
            time, but vastly improves accuracy. For approximately spherical
            molecules, or for large inter-molecular distances, this may be
            turned off
    """
    def __init__(self, number, molecules, numMols, thickness, volume_factor, allow_inversion=False, orientations=None, check_atomic_distances=True):
        
        #Necessary input
        self.number = number
        """The layer group number of the crystal."""
        self.Msgs()
        """A list of warning messages to use during generation."""
        numMols = np.array(numMols) #must convert it to np.array
        self.factor = volume_factor
        """The supplied volume factor for the unit cell."""
        self.numMols0 = numMols
        self.lgp = Layergroup(number)
        """The number (between 1 and 80) for the crystal's layer group."""
        self.sg = self.lgp.sgnumber
        """The number (between 1 and 230) for the international spacegroup."""
        #Reorient the molecules along their principle axes
        oriented_molecules = []
        #Allow support for generating molecules from text via openbabel
        for i, mol in enumerate(molecules):
            if type(mol) == str:
                #Read strings into molecules, try collection first,
                #If string not in collection, use SMILES format
                try:
                    mo = molecule_collection[mol]
                except:
                    mo = ob_mol_from_string(mol)
                    mo = pmg_from_ob(mo)
                molecules[i] = mo
        for mol in molecules:
            pga = PointGroupAnalyzer(mol)
            mo = pga.symmetrize_molecule()['sym_mol']
            oriented_molecules.append(mo)
        self.molecules = oriented_molecules
        """A list of pymatgen.core.structure.Molecule objects, symmetrized and
        oriented along their symmetry axes."""
        self.thickness = thickness
        """the thickness, in Angstroms, of the unit cell in the 3rd
        dimension."""
        self.PBC = [1,2]
        """The periodic axes of the crystal."""
        self.PB = self.lgp.permutation[3:6] 
        #TODO: add docstring
        self.P = self.lgp.permutation[:3] 
        #TODO: add docstring
        self.boxes = []
        """A list of bounding boxes for each molecule. Used for estimating
        volume of the unit cell."""
        self.radii = []
        """A list of approximated radii for each molecule type. Used for
        checking inter-molecular distances."""
        for mol in self.molecules:
            self.boxes.append(get_box(reoriented_molecule(mol)[0]))
            max_r = 0
            for site in mol:
                radius = math.sqrt( site.x**2 + site.y**2 + site.z**2 )
                if radius > max_r: max_r = radius
            self.radii.append(max_r+1.0)
        self.numMols = numMols * cellsize(self.sg)
        """The number of each type of molecule in the CONVENTIONAL cell"""
        self.volume = estimate_volume_molecular(self.numMols, self.boxes, self.factor)
        """The volume of the generated unit cell"""
        self.wyckoffs = get_layer(self.number)
        """The Wyckoff positions for the crystal's spacegroup."""
        self.wyckoffs_organized = get_layer(self.number, organized=True)
        """The Wyckoff positions for the crystal's spacegroup. Sorted by
        multiplicity."""
        self.w_symm = get_layer_symmetry(self.number, molecular=True)
        """A list of site symmetry operations for the Wyckoff positions, obtained
            from get_wyckoff_symmetry."""
        self.wyckoff_generators = get_layer_generators(self.number)
        """A list of Wyckfof generators for the layer group"""
        self.check_atomic_distances = check_atomic_distances
        """Whether or not inter-atomic distances are checked at each step."""
        self.allow_inversion = allow_inversion
        """Whether or not to allow chiral molecules to be inverted."""
        #When generating multiple crystals of the same stoichiometry and sg,
        #allow the user to re-use the allowed orientations, to reduce time cost
        if orientations is None:
            self.get_orientations()
        else:
            self.valid_orientations = orientations
            """The valid orientations for each molecule and Wyckoff position.
            May be copied when generating a new molecular_crystal to save a
            small amount of time"""
        self.generate_crystal()


    def Msgs(self):
        self.Msg1 = 'Error: the stoichiometry is incompatible with the wyckoff sites choice'
        self.Msg2 = 'Error: failed in the cycle of generating structures'
        self.Msg3 = 'Warning: failed in the cycle of adding species'
        self.Msg4 = 'Warning: failed in the cycle of choosing wyckoff sites'
        self.Msg5 = 'Finishing: added the specie'
        self.Msg6 = 'Finishing: added the whole structure'

    def get_orientations(self):
        """
        Calculates the valid orientations for each Molecule and Wyckoff
        position. Returns a list with 4 indices:

        index 1: the molecular prototype's index within self.molecules

        index 2: the Wyckoff position's 1st index (based on multiplicity)

        index 3: the WP's 2nd index (within the group of equal multiplicity)

        index 4: the index of the valid orientation for the molecule/WP pair

        For example, self.valid_orientations[i][j][k] would be a list of valid
        orientations for self.molecules[i], in the Wyckoff position
        self.wyckoffs_organized[j][k]
        """
        self.valid_orientations = []
        for mol in self.molecules:
            self.valid_orientations.append([])
            wp_index = -1
            for i, x in enumerate(self.wyckoffs_organized):
                self.valid_orientations[-1].append([])
                for j, wp in enumerate(x):
                    wp_index += 1
                    allowed = orientation_in_wyckoff_position(mol, self.wyckoffs, self.w_symm, wp_index, already_oriented=True, allow_inversion=self.allow_inversion)
                    if allowed is not False:
                        self.valid_orientations[-1][-1].append(allowed)
                    else:
                        self.valid_orientations[-1][-1].append([])

    def check_compatible(self):
        """
        Checks if the number of molecules is compatible with the Wyckoff
        positions. Considers the number of degrees of freedom for each Wyckoff
        position, and makes sure at least one valid combination of WP's exists.
        """
        N_site = [len(x[0]) for x in self.wyckoffs_organized]
        has_freedom = False
        #remove WP's with no freedom once they are filled
        removed_wyckoffs = []
        for i, numMol in enumerate(self.numMols):
            #Check that the number of molecules is a multiple of the smallest Wyckoff position
            if numMol % N_site[-1] > 0:
                return False
            else:
                #Check if smallest WP has at least one degree of freedom
                op = self.wyckoffs_organized[-1][-1][0]
                if op.rotation_matrix.all() != 0.0:
                    if self.valid_orientations[i][-1][-1] != []:
                        has_freedom = True
                else:
                    #Subtract from the number of ions beginning with the smallest Wyckoff positions
                    remaining = numMol
                    for j, x in enumerate(self.wyckoffs_organized):
                        for k, wp in enumerate(x):
                            while remaining >= len(wp) and wp not in removed_wyckoffs:
                                if self.valid_orientations[i][j][k] != []:
                                    #Check if WP has at least one degree of freedom
                                    op = wp[0]
                                    remaining -= len(wp)
                                    if np.allclose(op.rotation_matrix, np.zeros([3,3])):
                                        if (len(self.valid_orientations[i][j][k]) > 1 or
                                            self.valid_orientations[i][j][k][0].degrees > 0):
                                            #NOTE: degrees of freedom may be inaccurate for linear molecules
                                            has_freedom = True
                                        else:
                                            removed_wyckoffs.append(wp)
                                    else:
                                        has_freedom = True
                                else:
                                    removed_wyckoffs.append(wp)
                    if remaining != 0:
                        return False
        if has_freedom:
            return True
        else:
            #print("Warning: Wyckoff Positions have no degrees of freedom.")
            return 0

        return True

    def generate_crystal(self, max1=max1, max2=max2, max3=max3):
        """
        The main code to generate a random molecular crystal. If successful,
        stores a pymatgen.core.structure object in self.struct and sets
        self.valid to True. If unsuccessful, sets self.valid to False and
        outputs an error message.

        Args:
            max1: the number of attempts for generating a lattice
            max2: the number of attempts for a given lattice
            max3: the number of attempts for a given Wyckoff position
        """
        #Check the minimum number of degrees of freedom within the Wyckoff positions
        degrees = self.check_compatible()
        if degrees is False:
            print(self.Msg1)
            self.struct = None
            self.valid = False
            return
        else:
            if degrees == 0:
                max1 = 10
                max2 = 10
                max3 = 10
            #Calculate a minimum vector length for generating a lattice
            #minvector = max(radius*2 for radius in self.radii)
            all_lengths = []
            for box in self.boxes:
                lengths = []
                for a, b in [[0,1],[2,3],[4,5]]:
                    lengths.append(abs(box[b]-box[a]))
                all_lengths.append(min(lengths))
            minvector = max(all_lengths)
            for cycle1 in range(max1):
                #1, Generate a lattice
                cell_para = generate_lattice_2D(self.sg, self.volume, thickness=self.thickness, minvec=minvector)
                if cell_para is None:
                    break
                else:
                    cell_matrix = para2matrix(cell_para)
                    if abs(self.volume - np.linalg.det(cell_matrix)) > 1.0: 
                        print('Error, volume is not equal to the estimated value: ', self.volume, ' -> ', np.linalg.det(cell_matrix))
                        print('cell_para:  ', cell_para)
                        sys.exit(0)

                    molecular_coordinates_total = [] #to store the added molecular coordinates
                    molecular_sites_total = []      #to store the corresponding molecular specie
                    atomic_coordinates_total = [] #to store the added atomic coordinates
                    atomic_sites_total = []      #to store the corresponding atomic specie
                    wps_total = []      #to store corresponding Wyckoff position indices
                    points_total = []   #to store the generating x,y,z points
                    mol_generators_total = []
                    good_structure = False

                    for cycle2 in range(max2):
                        molecular_coordinates_tmp = deepcopy(molecular_coordinates_total)
                        molecular_sites_tmp = deepcopy(molecular_sites_total)
                        atomic_coordinates_tmp = deepcopy(atomic_coordinates_total)
                        atomic_sites_tmp = deepcopy(atomic_sites_total)
                        wps_tmp = deepcopy(wps_total)
                        points_tmp = deepcopy(points_total)
                        mol_generators_tmp = []
                        
                	    #Add molecules specie by specie
                        for numMol, mol in zip(self.numMols, self.molecules):
                            i = self.molecules.index(mol)
                            numMol_added = 0

                            #Now we start to add the specie to the wyckoff position
                            for cycle3 in range(max3):
                                #Choose a random Wyckoff position for given multiplicity: 2a, 2b, 2c
                                #NOTE: The molecular version return wyckoff indices, not ops
                                indices = choose_wyckoff_molecular(self.wyckoffs_organized, numMol-numMol_added, self.valid_orientations[i])
                                if indices is not False:
                                    j, k = indices
                                    if self.valid_orientations[i][j][k] == []:
                                        print("Error: Failed to catch empty set...")
                                        print(i,j,k)
                	    	    #Generate a list of coords from ops
                                    ops = self.wyckoffs_organized[j][k]
                                    point = np.random.random(3)
                                    coords = np.array([op.operate(point) for op in ops])
                                    #merge coordinates if the atoms are close
                                    if self.check_atomic_distances is False:
                                        mtol = self.radii[i]*2
                                    elif self.check_atomic_distances is True:
                                        mtol = 3.0
                                    coords_toadd, good_merge = merge_coordinate_molecular(coords, cell_matrix, self.wyckoffs, self.w_symm, mtol, self.valid_orientations[i], PBC=self.PBC)
                                    if good_merge is not False:
                                        wp_index = good_merge

                                        #Check that coords_toadd are generated by point
                                        generators = self.wyckoffs[wp_index]
                                        point = find_generating_point(coords_toadd, generators, PBC=self.PBC)
                                        if point is None:
                                            print("Error: Could not generate merged coordinates from Wyckoff generators")
                                            print("Layer group: "+str(self.number))
                                            print("wp_index: "+str(wp_index))
                                            print("Coordinates:")
                                            for c in coords_toadd:
                                                print(c)
                                            print("Generators:")
                                            for g in generators:
                                                print(g.as_xyz_string())
                                            self.valid = False
                                            return
                                        #Scale the coordinates to [0,1], except in the no-PBC direction
                                        coords_toadd = filtered_coords(coords_toadd, PBC=self.PBC)

                                        #Check inter-molecular distances
                                        if self.check_atomic_distances is False:
                                            if check_distance_molecular(molecular_coordinates_tmp, coords_toadd, molecular_sites_tmp, i, cell_matrix, self.radii, PBC=self.PBC):
                                                molecular_coordinates_tmp.append(coords_toadd)
                                                molecular_sites_tmp.append(i)
                                                wps_tmp.append(wp_index)
                                                points_tmp.append(point)
                                                numMol_added += len(coords_toadd)
                                                if numMol_added == numMol:
                                                    molecular_coordinates_total = deepcopy(molecular_coordinates_tmp)
                                                    molecular_sites_total = deepcopy(molecular_sites_tmp)
                                                    wps_total = deepcopy(wps_tmp)
                                                    points_total = deepcopy(points_tmp)
                                                    break

                                        #Check inter-atomic distances
                                        elif self.check_atomic_distances is True:
                                            #Generate atomic coordinates from molecules
                                            mo = deepcopy(self.molecules[i])
                                            #get j, k from wp_index
                                            num = 0
                                            found = False
                                            j, k = jk_from_i(wp_index, self.wyckoffs_organized)
                                            op1 = choose(self.valid_orientations[i][j][k]).get_op()
                                            mo.apply_operation(op1)
                                            ms0 = mol_site(mo, point, self.wyckoffs[wp_index], self.wyckoff_generators[wp_index], cell_matrix)
                                            wp_atomic_sites = [] #The species for the Wyckoff position
                                            wp_atomic_coords = [] #The coords for the Wyckoff position
                                            flag1 = True
                                            for point_index, op2 in enumerate(self.wyckoff_generators[wp_index]):
                                                current_atomic_sites = []
                                                current_atomic_coords = []
                                                for site in mo:
                                                    #Place molecular coordinates in relative coordinates
                                                    #relative_coords = np.dot(np.linalg.inv(np.transpose(cell_matrix)), site.coords)
                                                    #relative_coords = np.dot(np.linalg.inv(cell_matrix), site.coords)
                                                    relative_coords = np.dot(site.coords, np.linalg.inv(cell_matrix))
                                                    center1 = op2.operate(point)
                                                    rot = SymmOp.from_rotation_and_translation(op2.rotation_matrix,[0,0,0])
                                                    relative_coords = rot.operate(relative_coords)
                                                    new_vector = center1 + relative_coords
                                                    new_vector = filtered_coords(new_vector, PBC=self.PBC)
                                                    current_atomic_sites.append(site.specie.name)
                                                    current_atomic_coords.append(new_vector)
                                                wp_atomic_sites.append(current_atomic_sites)
                                                wp_atomic_coords.append(current_atomic_coords)
                                                #Check distances between molecules in current WP
                                                if point_index == 1:
                                                    for a_index, specie2 in enumerate(current_atomic_sites):
                                                        flag1 = check_distance([wp_atomic_coords[0]], [current_atomic_coords[a_index]], wp_atomic_sites[0], specie2, cell_matrix, d_factor=2.0, PBC=self.PBC)
                                                        if flag1 is False:
                                                            break
                                                    if flag1 is False:
                                                        break
                                            if flag1 is True:
                                                #Check distances between current and previous molecular atoms
                                                a = []
                                                for x in wp_atomic_coords:
                                                    a += x
                                                b = []
                                                for x in wp_atomic_sites:
                                                    b += x
                                                flag2 = True
                                                for a_index, specie2 in enumerate(b):
                                                    flag2 = check_distance([atomic_coordinates_total], [a[a_index]], atomic_sites_total, specie2, cell_matrix, d_factor=2.0, PBC=self.PBC)
                                                    if flag2 is False: break
                                                if flag2 is True:
                                                    mol_generators_tmp.append(ms0)
                                                    molecular_coordinates_tmp.append(coords_toadd)
                                                    molecular_sites_tmp.append(i)
                                                    atomic_coordinates_tmp += a
                                                    atomic_sites_tmp += b
                                                    wps_tmp.append(wp_index)
                                                    points_tmp.append(point)
                                                    numMol_added += len(coords_toadd)
                                                    if numMol_added == numMol:
                                                        mol_generators_total = deepcopy(mol_generators_tmp)
                                                        molecular_coordinates_total = deepcopy(molecular_coordinates_tmp)
                                                        atomic_sites_total = deepcopy(atomic_sites_tmp)
                                                        atomic_coordinates_total = deepcopy(atomic_coordinates_tmp)
                                                        molecular_sites_total = deepcopy(molecular_sites_tmp)
                                                        wps_total = deepcopy(wps_tmp)
                                                        points_total = deepcopy(points_tmp)
                                                        break

                            if numMol_added != numMol:
                                break  #need to repeat from the 1st species

                        if numMol_added == numMol:
                            #print(self.Msg6)
                            good_structure = True
                            break
                        else: #reset the coordinates and sites
                            molecular_coordinates_total = []
                            molecular_sites_total = []
                            wps_total = []
                    #placing molecules here
                    if good_structure:
                        final_lattice = cell_matrix 
                        final_coor = []
                        final_site = []
                        final_number = []

                        if self.check_atomic_distances is False:
                            for center0, i, wp_index in zip(points_total, molecular_sites_total, wps_total):
                                mo = deepcopy(self.molecules[i])
                                #get j, k from wp_index
                                num = 0
                                found = False
                                j, k = jk_from_i(wp_index, self.wyckoffs_organized)
                                op1 = choose(self.valid_orientations[i][j][k]).get_op()
                                mo.apply_operation(op1)
                                ms0 = mol_site(mo, center0, self.wyckoffs[wp_index], self.wyckoff_generators[wp_index], cell_matrix)
                                mol_generators_total.append(ms0)
                                for index, op2 in enumerate(self.wyckoff_generators[wp_index]):
                                    for site in mo:
                                        #Place molecular coordinates in relative coordinates
                                        relative_coords = np.dot(site.coords, np.linalg.inv(cell_matrix))
                                        center1 = op2.operate(center0)
                                        rot = SymmOp.from_rotation_and_translation(op2.rotation_matrix,[0,0,0])
                                        relative_coords = rot.operate(relative_coords)
                                        new_vector = center1 + relative_coords
                                        new_vector = filtered_coords(new_vector, PBC=self.PBC)
                                        final_coor.append(new_vector)
                                        final_site.append(site.specie.name)
                                        final_number.append(site.specie.number)
                            self.mol_generators = deepcopy(mol_generators_total)

                        elif self.check_atomic_distances is True:
                            final_coor = deepcopy(atomic_coordinates_total)
                            final_site = deepcopy(atomic_sites_total)
                            final_number = list(Element(ele).z for ele in atomic_sites_total)
                            self.mol_generators = deepcopy(mol_generators_total)
                            """A list of mol_site objects which can be used
                            for generating the crystal."""

                        final_coor = filtered_coords(final_coor, PBC=self.PBC)
                        if verify_distances(final_coor, final_site, final_lattice, factor=1.0, PBC=self.PBC) is True:
                            final_coor = np.array(final_coor)
                            #final_lattice, final_coor = Permutation(final_lattice, final_coor, self.PB)
                            final_lattice, final_coor = Add_vacuum(final_lattice, final_coor)
                            self.lattice = final_lattice
                            """A 3x3 matrix representing the lattice of the
                            unit cell."""  
                            self.coordinates = np.array(final_coor)
                            """The fractional coordinates for each molecule
                            in the final structure"""
                            self.sites = final_site
                            """The indices within self.molecules corresponding
                            to the type of molecule for each site in
                            self.coordinates."""              
                            self.struct = Structure(final_lattice, final_site, np.array(final_coor))
                            """A pymatgen.core.structure.Structure object for
                            the final generated crystal."""
                            self.spg_struct = (final_lattice, np.array(final_coor), final_number)
                            """A list of information describing the generated
                            crystal, which may be used by spglib for symmetry
                            analysis."""
                            self.valid = True
                            """Whether or not a valid crystal was generated."""
                            return
        print("Couldn't generate crystal after max attempts.")
        if degrees == 0:
            print("Note: Wyckoff positions have no degrees of freedom.")
        self.struct = self.Msg2
        self.valid = False
        return self.Msg2

class molecular_crystal_1D():
    """
    A 1d counterpart to molecular_crystal. Given a Rod group, list of
    molecule objects, molecular stoichiometry, volume factor, and area,
    generates a molecular crystal consistent with the given constraints.
    The crystal is stored as a pymatgen struct via self.struct
    
    Args:
        number: the Rod group number between 1 and 80. NOT equal to the
            international space group number, which is between 1 and 230
        molecules: a list of pymatgen.core.structure.Molecule objects for
            each type of molecule
        numMols: A list of the number of each type of molecule within the
            primitive cell (NOT the conventioal cell)
        area: cross-sectional area of the unit cell in Angstroms squared
        volume_factor: A volume factor used to generate a larger or smaller
            unit cell. Increasing this gives extra space between molecules
        allow_inversion: Whether or not to allow chiral molecules to be
            inverted. If True, the final crystal may contain mirror images of
            the original molecule. Unless the chemical properties of the mirror
            image are known, it is highly recommended to keep this value False
        orientations: Once a crystal with the same spacegroup and molecular
            stoichiometry has been generated, you may pass its
            valid_orientations attribute here to avoid repeating the
            calculation, but this is not required
        check_atomic_distances: If True, checks the inter-atomic distances
            after each Wyckoff position is added. This requires slightly more
            time, but vastly improves accuracy. For approximately spherical
            molecules, or for large inter-molecular distances, this may be
            turned off
    """
    def __init__(self, number, molecules, numMols, area, volume_factor, allow_inversion=False, orientations=None, check_atomic_distances=True):
        
        #Necessary input
        self.number = number
        """The number (between 1 and 80) for the crystal's Rod group."""
        self.Msgs()
        """A list of warning messages to use during generation."""
        numMols = np.array(numMols) #must convert it to np.array
        self.factor = volume_factor
        """The supplied volume factor for the unit cell."""
        self.numMols0 = numMols
        """The number of molecules of each type in the primitive cell"""
        #Reorient the molecules along their principle axes
        oriented_molecules = []
        #Allow support for generating molecules from text via openbabel
        for i, mol in enumerate(molecules):
            if type(mol) == str:
                #Read strings into molecules, try collection first,
                #If string not in collection, use SMILES format
                try:
                    mo = molecule_collection[mol]
                except:
                    mo = ob_mol_from_string(mol)
                    mo = pmg_from_ob(mo)
                molecules[i] = mo
        for mol in molecules:
            pga = PointGroupAnalyzer(mol)
            mo = pga.symmetrize_molecule()['sym_mol']
            oriented_molecules.append(mo)
        self.molecules = oriented_molecules
        """A list of pymatgen.core.structure.Molecule objects, symmetrized and
        oriented along their symmetry axes."""
        self.area = area
        """the effective cross-sectional area, in Angstroms squared, of the
        unit cell."""
        self.PBC = [3]
        self.boxes = []
        """A list of bounding boxes for each molecule. Used for estimating
        volume of the unit cell."""
        self.radii = []
        """A list of approximated radii for each molecule type. Used for
        checking inter-molecular distances."""
        for mol in self.molecules:
            self.boxes.append(get_box(reoriented_molecule(mol)[0]))
            max_r = 0
            for site in mol:
                radius = math.sqrt( site.x**2 + site.y**2 + site.z**2 )
                if radius > max_r: max_r = radius
            self.radii.append(max_r+1.0)
        self.numMols = numMols
        """The number of each type of molecule in the CONVENTIONAL cell"""
        self.volume = estimate_volume_molecular(self.numMols, self.boxes, self.factor)
        """The volume of the generated unit cell"""
        self.wyckoffs = get_rod(self.number)
        """The Wyckoff positions for the crystal's Rod group."""
        self.wyckoffs_organized = get_rod(self.number, organized=True)
        """The Wyckoff positions for the crystal's Rod group. Sorted by
        multiplicity."""
        self.w_symm = get_rod_symmetry(self.number, molecular=True)
        """A list of site symmetry operations for the Wyckoff positions, obtained
            from get_wyckoff_symmetry."""
        self.wyckoff_generators = get_rod_generators(self.number)
        self.check_atomic_distances = check_atomic_distances
        """Whether or not inter-atomic distances are checked at each step."""
        self.allow_inversion = allow_inversion
        """Whether or not to allow chiral molecules to be inverted."""
        #When generating multiple crystals of the same stoichiometry and sg,
        #allow the user to re-use the allowed orientations, to reduce time cost
        if orientations is None:
            self.get_orientations()
        else:
            self.valid_orientations = orientations
            """The valid orientations for each molecule and Wyckoff position.
            May be copied when generating a new molecular_crystal to save a
            small amount of time"""
        self.generate_crystal()


    def Msgs(self):
        self.Msg1 = 'Error: the stoichiometry is incompatible with the wyckoff sites choice'
        self.Msg2 = 'Error: failed in the cycle of generating structures'
        self.Msg3 = 'Warning: failed in the cycle of adding species'
        self.Msg4 = 'Warning: failed in the cycle of choosing wyckoff sites'
        self.Msg5 = 'Finishing: added the specie'
        self.Msg6 = 'Finishing: added the whole structure'

    def get_orientations(self):
        """
        Calculates the valid orientations for each Molecule and Wyckoff
        position. Returns a list with 4 indices:

        index 1: the molecular prototype's index within self.molecules

        index 2: the Wyckoff position's 1st index (based on multiplicity)

        index 3: the WP's 2nd index (within the group of equal multiplicity)

        index 4: the index of the valid orientation for the molecule/WP pair

        For example, self.valid_orientations[i][j][k] would be a list of valid
        orientations for self.molecules[i], in the Wyckoff position
        self.wyckoffs_organized[j][k]
        """
        self.valid_orientations = []
        for mol in self.molecules:
            self.valid_orientations.append([])
            wp_index = -1
            for i, x in enumerate(self.wyckoffs_organized):
                self.valid_orientations[-1].append([])
                for j, wp in enumerate(x):
                    wp_index += 1
                    allowed = orientation_in_wyckoff_position(mol, self.wyckoffs, self.w_symm, wp_index, already_oriented=True, allow_inversion=self.allow_inversion)
                    if allowed is not False:
                        self.valid_orientations[-1][-1].append(allowed)
                    else:
                        self.valid_orientations[-1][-1].append([])

    def check_compatible(self):
        """
        Checks if the number of molecules is compatible with the Wyckoff
        positions. Considers the number of degrees of freedom for each Wyckoff
        position, and makes sure at least one valid combination of WP's exists.
        """
        N_site = [len(x[0]) for x in self.wyckoffs_organized]
        has_freedom = False
        #remove WP's with no freedom once they are filled
        removed_wyckoffs = []
        for i, numMol in enumerate(self.numMols):
            #Check that the number of molecules is a multiple of the smallest Wyckoff position
            if numMol % N_site[-1] > 0:
                return False
            else:
                #Check if smallest WP has at least one degree of freedom
                op = self.wyckoffs_organized[-1][-1][0]
                if op.rotation_matrix.all() != 0.0:
                    if self.valid_orientations[i][-1][-1] != []:
                        has_freedom = True
                else:
                    #Subtract from the number of ions beginning with the smallest Wyckoff positions
                    remaining = numMol
                    for j, x in enumerate(self.wyckoffs_organized):
                        for k, wp in enumerate(x):
                            while remaining >= len(wp) and wp not in removed_wyckoffs:
                                if self.valid_orientations[i][j][k] != []:
                                    #Check if WP has at least one degree of freedom
                                    op = wp[0]
                                    remaining -= len(wp)
                                    if np.allclose(op.rotation_matrix, np.zeros([3,3])):
                                        if (len(self.valid_orientations[i][j][k]) > 1 or
                                            self.valid_orientations[i][j][k][0].degrees > 0):
                                            #NOTE: degrees of freedom may be inaccurate for linear molecules
                                            has_freedom = True
                                        else:
                                            removed_wyckoffs.append(wp)
                                    else:
                                        has_freedom = True
                                else:
                                    removed_wyckoffs.append(wp)
                    if remaining != 0:
                        return False
        if has_freedom:
            return True
        else:
            #print("Warning: Wyckoff Positions have no degrees of freedom.")
            return 0

        return True

    def generate_crystal(self, max1=max1, max2=max2, max3=max3):
        """
        The main code to generate a random molecular crystal. If successful,
        stores a pymatgen.core.structure object in self.struct and sets
        self.valid to True. If unsuccessful, sets self.valid to False and
        outputs an error message.

        Args:
            max1: the number of attempts for generating a lattice
            max2: the number of attempts for a given lattice
            max3: the number of attempts for a given Wyckoff position
        """
        #Check the minimum number of degrees of freedom within the Wyckoff positions
        degrees = self.check_compatible()
        if degrees is False:
            print(self.Msg1)
            self.struct = None
            self.valid = False
            return
        else:
            if degrees == 0:
                max1 = 10
                max2 = 10
                max3 = 10
            #Calculate a minimum vector length for generating a lattice
            #minvector = max(radius*2 for radius in self.radii)
            all_lengths = []
            for box in self.boxes:
                lengths = []
                for a, b in [[0,1],[2,3],[4,5]]:
                    lengths.append(abs(box[b]-box[a]))
                all_lengths.append(min(lengths))
            minvector = max(max(all_lengths), tol_m)
            for cycle1 in range(max1):
                #1, Generate a lattice
                cell_para = generate_lattice_1D(self.number, self.volume, area=self.area, minvec=minvector)
                if cell_para is None:
                    break
                else:
                    cell_matrix = para2matrix(cell_para)
                    if abs(self.volume - np.linalg.det(cell_matrix)) > 1.0: 
                        print('Error, volume is not equal to the estimated value: ', self.volume, ' -> ', np.linalg.det(cell_matrix))
                        print('cell_para:  ', cell_para)
                        sys.exit(0)

                    molecular_coordinates_total = [] #to store the added molecular coordinates
                    molecular_sites_total = []      #to store the corresponding molecular specie
                    atomic_coordinates_total = [] #to store the added atomic coordinates
                    atomic_sites_total = []      #to store the corresponding atomic specie
                    wps_total = []      #to store corresponding Wyckoff position indices
                    points_total = []   #to store the generating x,y,z points
                    mol_generators_total = []
                    good_structure = False

                    for cycle2 in range(max2):
                        molecular_coordinates_tmp = deepcopy(molecular_coordinates_total)
                        molecular_sites_tmp = deepcopy(molecular_sites_total)
                        atomic_coordinates_tmp = deepcopy(atomic_coordinates_total)
                        atomic_sites_tmp = deepcopy(atomic_sites_total)
                        wps_tmp = deepcopy(wps_total)
                        points_tmp = deepcopy(points_total)
                        mol_generators_tmp = []
                        
                	    #Add molecules specie by specie
                        for numMol, mol in zip(self.numMols, self.molecules):
                            i = self.molecules.index(mol)
                            numMol_added = 0

                            #Now we start to add the specie to the wyckoff position
                            for cycle3 in range(max3):
                                #Choose a random Wyckoff position for given multiplicity: 2a, 2b, 2c
                                #NOTE: The molecular version return wyckoff indices, not ops
                                indices = choose_wyckoff_molecular(self.wyckoffs_organized, numMol-numMol_added, self.valid_orientations[i])
                                if indices is not False:
                                    j, k = indices
                                    if self.valid_orientations[i][j][k] == []:
                                        print("Error: Failed to catch empty set...")
                                        print(i,j,k)
                	    	    #Generate a list of coords from ops
                                    ops = self.wyckoffs_organized[j][k]
                                    point = np.random.random(3)
                                    coords = np.array([op.operate(point) for op in ops])
                                    #merge coordinates if the atoms are close
                                    if self.check_atomic_distances is False:
                                        mtol = self.radii[i]*2
                                    elif self.check_atomic_distances is True:
                                        mtol = 3.0
                                    coords_toadd, good_merge = merge_coordinate_molecular(coords, cell_matrix, self.wyckoffs, self.w_symm, mtol, self.valid_orientations[i], PBC=self.PBC)
                                    if good_merge is not False:
                                        wp_index = good_merge

                                        #Check that coords_toadd are generated by point
                                        generators = self.wyckoffs[wp_index]
                                        point = find_generating_point(coords_toadd, generators, PBC=self.PBC)
                                        if point is None:
                                            print("Error: Could not generate merged coordinates from Wyckoff generators")
                                            print("Layer group: "+str(self.number))
                                            print("wp_index: "+str(wp_index))
                                            print("Coordinates:")
                                            for c in coords_toadd:
                                                print(c)
                                            print("Generators:")
                                            for g in generators:
                                                print(g.as_xyz_string())
                                            self.valid = False
                                            return
                                        #Scale the coordinates to [0,1], except in the no-PBC direction
                                        coords_toadd = filtered_coords(coords_toadd, PBC=self.PBC)

                                        #Check inter-molecular distances
                                        if self.check_atomic_distances is False:
                                            if check_distance_molecular(molecular_coordinates_tmp, coords_toadd, molecular_sites_tmp, i, cell_matrix, self.radii, PBC=self.PBC):
                                                molecular_coordinates_tmp.append(coords_toadd)
                                                molecular_sites_tmp.append(i)
                                                wps_tmp.append(wp_index)
                                                points_tmp.append(point)
                                                numMol_added += len(coords_toadd)
                                                if numMol_added == numMol:
                                                    molecular_coordinates_total = deepcopy(molecular_coordinates_tmp)
                                                    molecular_sites_total = deepcopy(molecular_sites_tmp)
                                                    wps_total = deepcopy(wps_tmp)
                                                    points_total = deepcopy(points_tmp)
                                                    break

                                        #Check inter-atomic distances
                                        elif self.check_atomic_distances is True:
                                            #Generate atomic coordinates from molecules
                                            mo = deepcopy(self.molecules[i])
                                            #get j, k from wp_index
                                            num = 0
                                            found = False
                                            j, k = jk_from_i(wp_index, self.wyckoffs_organized)
                                            op1 = choose(self.valid_orientations[i][j][k]).get_op()
                                            mo.apply_operation(op1)
                                            ms0 = mol_site(mo, point, self.wyckoffs[wp_index], self.wyckoff_generators[wp_index], cell_matrix)
                                            wp_atomic_sites = [] #The species for the Wyckoff position
                                            wp_atomic_coords = [] #The coords for the Wyckoff position
                                            flag1 = True
                                            for point_index, op2 in enumerate(self.wyckoff_generators[wp_index]):
                                                current_atomic_sites = []
                                                current_atomic_coords = []
                                                for site in mo:
                                                    #Place molecular coordinates in relative coordinates
                                                    #relative_coords = np.dot(np.linalg.inv(np.transpose(cell_matrix)), site.coords)
                                                    #relative_coords = np.dot(np.linalg.inv(cell_matrix), site.coords)
                                                    relative_coords = np.dot(site.coords, np.linalg.inv(cell_matrix))
                                                    center1 = op2.operate(point)
                                                    rot = SymmOp.from_rotation_and_translation(op2.rotation_matrix,[0,0,0])
                                                    relative_coords = rot.operate(relative_coords)
                                                    new_vector = center1 + relative_coords
                                                    new_vector = filtered_coords(new_vector, PBC=self.PBC)
                                                    current_atomic_sites.append(site.specie.name)
                                                    current_atomic_coords.append(new_vector)
                                                wp_atomic_sites.append(current_atomic_sites)
                                                wp_atomic_coords.append(current_atomic_coords)
                                                #Check distances between molecules in current WP
                                                if point_index == 1:
                                                    for a_index, specie2 in enumerate(current_atomic_sites):
                                                        flag1 = check_distance([wp_atomic_coords[0]], [current_atomic_coords[a_index]], wp_atomic_sites[0], specie2, cell_matrix, d_factor=2.0, PBC=self.PBC)
                                                        if flag1 is False:
                                                            break
                                                    if flag1 is False:
                                                        break
                                            if flag1 is True:
                                                #Check distances between current and previous molecular atoms
                                                a = []
                                                for x in wp_atomic_coords:
                                                    a += x
                                                b = []
                                                for x in wp_atomic_sites:
                                                    b += x
                                                flag2 = True
                                                for a_index, specie2 in enumerate(b):
                                                    flag2 = check_distance([atomic_coordinates_total], [a[a_index]], atomic_sites_total, specie2, cell_matrix, d_factor=2.0, PBC=self.PBC)
                                                    if flag2 is False: break
                                                if flag2 is True:
                                                    mol_generators_tmp.append(ms0)
                                                    molecular_coordinates_tmp.append(coords_toadd)
                                                    molecular_sites_tmp.append(i)
                                                    atomic_coordinates_tmp += a
                                                    atomic_sites_tmp += b
                                                    wps_tmp.append(wp_index)
                                                    points_tmp.append(point)
                                                    numMol_added += len(coords_toadd)
                                                    if numMol_added == numMol:
                                                        mol_generators_total = deepcopy(mol_generators_tmp)
                                                        molecular_coordinates_total = deepcopy(molecular_coordinates_tmp)
                                                        atomic_sites_total = deepcopy(atomic_sites_tmp)
                                                        atomic_coordinates_total = deepcopy(atomic_coordinates_tmp)
                                                        molecular_sites_total = deepcopy(molecular_sites_tmp)
                                                        wps_total = deepcopy(wps_tmp)
                                                        points_total = deepcopy(points_tmp)
                                                        break

                            if numMol_added != numMol:
                                break  #need to repeat from the 1st species

                        if numMol_added == numMol:
                            #print(self.Msg6)
                            good_structure = True
                            break
                        else: #reset the coordinates and sites
                            molecular_coordinates_total = []
                            molecular_sites_total = []
                            wps_total = []
                    #placing molecules here
                    if good_structure:
                        final_lattice = cell_matrix 
                        final_coor = []
                        final_site = []
                        final_number = []

                        if self.check_atomic_distances is False:
                            for center0, i, wp_index in zip(points_total, molecular_sites_total, wps_total):
                                mo = deepcopy(self.molecules[i])
                                #get j, k from wp_index
                                num = 0
                                found = False
                                j, k = jk_from_i(wp_index, self.wyckoffs_organized)
                                op1 = choose(self.valid_orientations[i][j][k]).get_op()
                                mo.apply_operation(op1)
                                ms0 = mol_site(mo, center0, self.wyckoffs[wp_index], self.wyckoff_generators[wp_index], cell_matrix)
                                mol_generators_total.append(ms0)
                                for index, op2 in enumerate(self.wyckoff_generators[wp_index]):
                                    for site in mo:
                                        #Place molecular coordinates in relative coordinates
                                        relative_coords = np.dot(site.coords, np.linalg.inv(cell_matrix))
                                        center1 = op2.operate(center0)
                                        rot = SymmOp.from_rotation_and_translation(op2.rotation_matrix,[0,0,0])
                                        relative_coords = rot.operate(relative_coords)
                                        new_vector = center1 + relative_coords
                                        new_vector = filtered_coords(new_vector, PBC=self.PBC)
                                        final_coor.append(new_vector)
                                        final_site.append(site.specie.name)
                                        final_number.append(site.specie.number)
                            self.mol_generators = deepcopy(mol_generators_total)

                        elif self.check_atomic_distances is True:
                            final_coor = deepcopy(atomic_coordinates_total)
                            final_site = deepcopy(atomic_sites_total)
                            final_number = list(Element(ele).z for ele in atomic_sites_total)
                            self.mol_generators = deepcopy(mol_generators_total)
                            """A list of mol_site objects which can be used
                            for generating the crystal."""

                        final_coor = filtered_coords(final_coor, PBC=self.PBC)
                        if verify_distances(final_coor, final_site, final_lattice, factor=1.0, PBC=self.PBC) is True:
                            final_coor = np.array(final_coor)
                            #final_lattice, final_coor = Permutation(final_lattice, final_coor, self.PB)
                            final_lattice, final_coor = Add_vacuum(final_lattice, final_coor)
                            self.lattice = final_lattice
                            """A 3x3 matrix representing the lattice of the
                            unit cell."""  
                            self.coordinates = np.array(final_coor)
                            """The fractional coordinates for each molecule
                            in the final structure"""
                            self.sites = final_site
                            """The indices within self.molecules corresponding
                            to the type of molecule for each site in
                            self.coordinates."""              
                            self.struct = Structure(final_lattice, final_site, np.array(final_coor))
                            """A pymatgen.core.structure.Structure object for
                            the final generated crystal."""
                            self.spg_struct = (final_lattice, np.array(final_coor), final_number)
                            """A list of information describing the generated
                            crystal, which may be used by spglib for symmetry
                            analysis."""
                            self.valid = True
                            """Whether or not a valid crystal was generated."""
                            return
        print("Couldn't generate crystal after max attempts.")
        if degrees == 0:
            print("Note: Wyckoff positions have no degrees of freedom.")
        self.struct = self.Msg2
        self.valid = False
        return self.Msg2


if __name__ == "__main__":
    #-------------------------------- Options -------------------------
    import os

    parser = OptionParser()
    parser.add_option("-s", "--spacegroup", dest="sg", metavar='sg', default=36, type=int,
            help="desired space group number: 1-230, e.g., 36")
    parser.add_option("-e", "--molecule", dest="molecule", default='H2O', 
            help="desired molecules: e.g., H2O", metavar="molecule")
    parser.add_option("-n", "--numMols", dest="numMols", default=4, 
            help="desired numbers of molecules: 4", metavar="numMols")
    parser.add_option("-f", "--factor", dest="factor", default=3.0, type=float, 
            help="volume factor: default 3.0", metavar="factor")
    parser.add_option("-v", "--verbosity", dest="verbosity", default=0, type=int, help="verbosity: default 0; higher values print more information", metavar="verbosity")
    parser.add_option("-a", "--attempts", dest="attempts", default=1, type=int, 
            help="number of crystals to generate: default 1", metavar="attempts")
    parser.add_option("-o", "--outdir", dest="outdir", default="out", type=str, 
            help="Directory for storing output cif files: default 'out'", metavar="outdir")
    parser.add_option("-c", "--checkatoms", dest="checkatoms", default="True", type=str, 
            help="Whether to check inter-atomic distances at each step: default True", metavar="outdir")
    parser.add_option("-i", "--allowinversion", dest="allowinversion", default="False", type=str, 
            help="Whether to allow inversion of chiral molecules: default False", metavar="outdir")
    parser.add_option("-d", "--dimension", dest="dimension", metavar='dimension', default=3, type=int,
            help="desired dimension: (3 or 2 for 3d or 2d, respectively)")
    parser.add_option("-t", "--thickness", dest="thickness", metavar='thickness', default=2.0, type=float,
            help="Thickness, in Angstroms, of a 2D crystal: default 2.0")

    (options, args) = parser.parse_args()    
    molecule = options.molecule
    number = options.numMols
    verbosity = options.verbosity
    attempts = options.attempts
    outdir = options.outdir
    factor = options.factor
    dimension = options.dimension
    thickness = options.thickness

    if options.checkatoms == "True" or options.checkatoms == "False":
        checkatoms = eval(options.checkatoms)
    else:
        print("Invalid value for -c (--checkatoms): must be 'True' or 'False'.")
        checkatoms = True
    if options.allowinversion == "True" or options.allowinversion == "False":
        allowinversion = eval(options.allowinversion)
    else:
        print("Invalid value for -i (--allowinversion): must be 'True' or 'False'.")
        allowinversion = False
    
    numMols = []
    if molecule.find(',') > 0:
        strings = molecule.split(',')
        system = []
        for mol in strings:
            system.append(ob_mol_from_string(mol))
        for x in number.split(','):
            numMols.append(int(x))
    else:
        system = [ob_mol_from_string(molecule)]
        numMols = [int(number)]
    orientations = None

    try:
        os.mkdir(outdir)
    except: pass

    filecount = 1 #To check whether a file already exists
    for i in range(attempts):
        start = time()
        numMols0 = np.array(numMols)
        sg = options.sg
        if dimension == 3:
            rand_crystal = molecular_crystal(options.sg, system, numMols0, factor, check_atomic_distances=checkatoms, allow_inversion=allowinversion)
        elif dimension == 2:
            rand_crystal = molecular_crystal_2D(options.sg, system, numMols0, thickness, factor, allow_inversion=allowinversion, check_atomic_distances=checkatoms)
        end = time()
        timespent = np.around((end - start), decimals=2)
        if rand_crystal.valid:
            #Output a cif file
            written = False
            try:
                comp = str(rand_crystal.struct.composition)
                comp = comp.replace(" ", "")
                cifpath = outdir + '/' + comp + "_" + str(filecount) + '.cif'
                while os.path.isfile(cifpath):
                    filecount += 1
                    cifpath = outdir + '/' + comp + "_" + str(filecount) + '.cif'
                CifWriter(rand_crystal.struct, symprec=0.1).write_file(filename = cifpath)
                written = True
            except: pass

            #spglib style structure called cell
            ans = get_symmetry_dataset(rand_crystal.spg_struct, symprec=1e-1)['number']
            print('Space group requested: ', sg, 'generated', ans, 'vol: ', rand_crystal.volume)
            if written is True:
                print("    Output to "+cifpath)
            else:
                print("    Could not write cif file.")

            #Print additional information about the structure
            if verbosity > 0:
                print("Time required for generation: " + str(timespent) + "s")
                print("Molecular Wyckoff positions:")
                for ms in rand_crystal.mol_generators:
                    print(str(ms.mol.composition) + ": " + str(ms.multiplicity)+str(ms.letter)+" "+str(ms.position))
            if verbosity > 1:
                print(rand_crystal.struct)

        #If generation fails
        else: 
            print('something is wrong')
            print('Time spent during generation attempt: ' + str(timespent) + "s")
