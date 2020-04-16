#PyXtal imports
from pyxtal.crystal import *
from pyxtal.molecule import *

#Constants
#------------------------------
molecule_collection = Collection('molecules')
max1 = 40 #Attempts for generating lattices
max2 = 10 #Attempts for a given lattice
max3 = 10 #Attempts for a given Wyckoff position
max4 = 3 #Attempts for a given mol_site (changning orientation)

max_fast_mol_size = 30 #Max number of atoms per molecule before using fast distance check

#Define functions
#------------------------------
def check_intersection(ellipsoid1, ellipsoid2):
    """
    Given SymmOp's for 2 ellipsoids, checks whether or not they overlap

    Args:
        ellipsoid1: a SymmOp representing the first ellipsoid
        ellipsoid2: a SymmOp representing the second ellipsoid

    Returns:
        False if the ellipsoids overlap.
        True if they do not overlap.
    """
    #Transform so that one ellipsoid becomes a unit sphere at (0,0,0)
    Op = ellipsoid1.inverse * ellipsoid2
    #We define a new ellipsoid by moving the sphere around the old ellipsoid
    M = Op.rotation_matrix
    a = 1.0 /(1.0 / np.linalg.norm(M[0]) + 1)
    M[0] = M[0] / np.linalg.norm(M[0]) * a
    b = 1.0 / (1.0 / np.linalg.norm(M[1]) + 1)
    M[1] = M[1] / np.linalg.norm(M[1]) * b
    c = 1.0 / (1.0 / np.linalg.norm(M[2]) + 1)
    M[2] = M[2] / np.linalg.norm(M[2]) * c
    p = Op.translation_vector
    #Calculate the transformed distance from the sphere's center to the new ellipsoid
    dsq = np.dot(p, M[0])**2 + np.dot(p, M[1])**2 + np.dot(p, M[2])**2
    if dsq < 2:
        return False
    else:
        return True

def check_mol_sites(ms1, ms2, atomic=False, factor=1.0, tm=Tol_matrix(prototype="molecular")):
    """
    Checks whether or not the molecules of two mol sites overlap. Uses
    ellipsoid overlapping approximation to check. Takes PBC and lattice
    into consideration.

    Args:
        ms1: a mol_site object
        ms2: another mol_site object
        atomic: if True, checks inter-atomic distances. If False, checks
            overlap between molecular ellipsoids
        factor: the distance factor to pass to check_distances. (only for
            inter-atomic distance checking)
        tm: a Tol_matrix object (or prototype string) for distance checking

    Returns:
        False if the Wyckoff positions overlap. True otherwise
    """
    if atomic is False:
        es0 = ms1.get_ellipsoids()
        PBC_vectors = np.dot(create_matrix(PBC=ms1.PBC), ms1.lattice)
        PBC_ops = [SymmOp.from_rotation_and_translation(Euclidean_lattice, v) for v in PBC_vectors]
        es1 = []
        for op in PBC_ops:
            es1.append(np.dot(es0, op))
        es1 = np.squeeze(es1)
        truth_values = np.vectorize(check_intersection)(es1, ms2.get_ellipsoid())
        if np.sum(truth_values) < len(truth_values):
            return False
        else:
            return True

    elif atomic is True:
        #Get coordinates for both mol_sites
        c1, s1 = ms1.get_coords_and_species()
        c2, s2 = ms2.get_coords_and_species()

        #Calculate which distance matrix is smaller/faster
        m_length1 = len(ms1.mol.species)
        m_length2 = len(ms2.mol.species)
        wp_length1 = len(c1)
        wp_length2 = len(c2)
        size1 = m_length1 * wp_length2
        size2 = m_length2 * wp_length1

        #Case 1
        if size1 <= size2:
            coords_mol = c1[:m_length1]
            #Calculate tol matrix for species pairs
            tols = np.zeros((m_length1,m_length2))
            for i1, specie1 in enumerate(s1[:m_length1]):
                for i2, specie2 in enumerate(s2[:m_length2]):
                    tols[i1][i2] = tm.get_tol(specie1, specie2)
            tols = np.repeat(tols, ms2.multiplicity, axis=1)
            d = distance_matrix(coords_mol, c2, ms1.lattice, PBC=ms1.PBC)

        #Case 2
        elif size1 > size2:
            coords_mol = c2[:m_length2]
            #Calculate tol matrix for species pairs
            tols = np.zeros((m_length2,m_length1))
            for i1, specie1 in enumerate(s2[:m_length2]):
                for i2, specie2 in enumerate(s1[:m_length1]):
                    tols[i1][i2] = tm.get_tol(specie1, specie2)
            tols = np.repeat(tols, ms1.multiplicity, axis=1)
            d = distance_matrix(coords_mol, c1, ms1.lattice, PBC=ms1.PBC)

        #Check if distances are smaller than tolerances
        if (d<tols).any():
            return False
        return True

def estimate_volume_molecular(molecules, numMols, factor=2.0, boxes=None):
    """
    Given the molecular stoichiometry, estimate the volume needed for a unit cell.

    Args:
        molecules: a list of Pymatgen Molecule objects
        numMols: a list with the number of each type of molecule
        factor: a factor to multiply the final result by. Used to increase space
        between molecules
        boxes: a list of Box objects for each molecule. Obtained from get_box
            if None, boxes are calculated automatically.

    Returns:
        the estimated volume (in cubic Angstroms) needed for the unit cell
    """
    if boxes is None:
        boxes = []
        for mol in molecules:
            boxes.append(get_box(reoriented_molecule(mol)[0]))
    volume = 0
    for numMol, box in zip(numMols, boxes):
        volume += numMol*box.volume
    return abs(factor*volume)

def get_group_orientations(mol, group, allow_inversion=False):
    """
    Calculate the valid orientations for each Molecule and Wyckoff position.
    Returns a list with 3 indices:

        index 1: the Wyckoff position's 1st index (based on multiplicity)  

        index 2: the WP's 2nd index (within the group of equal multiplicity)  

        index 3: the index of the valid orientation for the molecule/WP pair

    For example, self.valid_orientations[i][j] would be a list of valid
    orientations for self.molecules[i], in the Wyckoff position
    self.group.wyckoffs_organized[i][j]

    Args:
        mol: a pymatgen Molecule object.
        group: a Group object
        allow_inversion: whether or not to allow inversion operations for chiral
            molecules

    Returns:
        a list of operations orientation objects for each Wyckoff position. 1st
            and 2nd indices correspond to the Wyckoff position
    """
    wyckoffs = group.wyckoffs_organized
    w_symm_all = group.w_symm_m
    valid_orientations = []
    wp_index = -1
    for i, x in enumerate(wyckoffs):
        valid_orientations.append([])
        for j, wp in enumerate(x):
            wp_index += 1
            allowed = orientation_in_wyckoff_position(mol, wp, allow_inversion=allow_inversion)
            if allowed is not False:
                valid_orientations[-1].append(allowed)
            else:
                valid_orientations[-1].append([])
    return valid_orientations

class Box():
    """
    Class for storing the binding box for a molecule. Box is oriented along the x, y, and
    z axes.

    Args:
        minx: the minimum x value
        maxx: the maximum x value
        miny: the minimum y value
        maxy: the maximum y value
        minz: the minimum z value
        maxz: the maximum z value
    """
    def __init__(self, minx, maxx, miny, maxy, minz, maxz):
        self.minx = float(minx)
        self.maxx = float(maxx)
        self.miny = float(miny)
        self.maxy = float(maxy)
        self.minz = float(minz)
        self.maxz = float(maxz)

        self.width = float(abs(maxx - minx))
        self.length = float(abs(maxy - miny))
        self.height = float(abs(maxz - minz))

        self.minl = min(self.width, self.length, self.height)
        self.maxl = max(self.width, self.length, self.height)
        for x in (self.width, self.length, self.height):
            if x <= self.maxl and x >= self.minl:
                self.midl = x

        self.volume = float(self.width * self.length * self.height)

def get_box(mol):
    """
    Given a molecule, find a minimum orthorhombic box containing it.
    Size is calculated using min and max x, y, and z values, plus the padding defined by the vdw radius
    For best results, call oriented_molecule first.
    
    Args:
        mol: a pymatgen Molecule object. Should be oriented along its principle axes.

    Returns:
        a Box object
    """
    minx, miny, minz, maxx, maxy, maxz = 0.,0.,0.,0.,0.,0.
    #for p in mol:
    for p in mol:
        x, y, z = p.coords
        r = Element(p.species_string).vdw_radius
        if x-r < minx: minx = x-r
        if y-r < miny: miny = y-r
        if z-r < minz: minz = z-r
        if x+r > maxx: maxx = x+r
        if y+r > maxy: maxy = y+r
        if z+r > maxz: maxz = z+r
    return Box(minx,maxx,miny,maxy,minz,maxz)

def check_distance_molecular(coord1, coord2, indices1, index2, lattice, radii, d_factor=1.0, PBC=[1,1,1]):
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
        PBC: A periodic boundary condition list, where 1 means periodic, 0 means not periodic.
            Ex: [1,1,1] -> full 3d periodicity, [0,0,1] -> periodicity along the z axis

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

def merge_coordinate_molecular(coor, lattice, group, tol, orientations):
    """
    Given a list of fractional coordinates, merges them within a given
    tolerance, and checks if the merged coordinates satisfy a Wyckoff
    position. Used for merging general Wyckoff positions into special Wyckoff
    positions within the random_crystal (and its derivative) classes.

    Args:
        coor: a list of fractional coordinates
        lattice: a 3x3 matrix representing the unit cell
        group: a pyxtal.symmetry.Group object
        tol: the cutoff distance for merging coordinates
        orientations: the valid orientations for a given molecule. Obtained
            from get_sg_orientations, which is called within molecular_crystal

    Returns:
        coor, index, point: (coor) is the new list of fractional coordinates after
        merging, and index is a single index of the Wyckoff position within
        the spacegroup. If merging is unsuccesful, or no index is found,
        returns the original coordinates and False. point is a 3-vector which can
        be plugged into the Wyckoff position to generate the rest of the points
    """
    #Get index of current Wyckoff position. If not one, return False
    index, point = check_wyckoff_position(coor, group)
    if index is False:
        return coor, False, None
    if point is None:
        return coor, False, None
    PBC = group.PBC
    #Main loop for merging multiple times
    while True:
        #Check distances of current WP. If too small, merge
        dm = distance_matrix([coor[0]], coor, lattice, PBC=PBC)
        passed_distance_check = True
        for i, x in enumerate(dm):
            for j, y in enumerate(x):
                if i != j and y < tol:
                    passed_distance_check = False
                    break
            if passed_distance_check is False:
                break
        if check_images([coor[0]], ['C'], lattice, PBC=PBC, tol=tol) is False:
            passed_distance_check = False
        if passed_distance_check is False:
            mult1 = group[index].multiplicity
            #Find possible wp's to merge into
            possible = []
            for i, wp in enumerate(group):
                mult2 = wp.multiplicity
                #Check that a valid orientation exists
                j, k = jk_from_i(i, orientations)
                if orientations[j][k] == []: continue
                #Only allow smaller WP's that are an integer factor of the current one
                if (mult2 < mult1) and (mult1 % mult2 == 0):
                    possible.append(i)
            if possible == []:
                return coor, False, None
            #Calculate minimum separation for each WP
            distances = []
            for i in possible:
                wp = group[i]
                projected_point = project_point(point, wp[0], lattice=lattice, PBC=PBC)
                new_coor = apply_ops(projected_point, wp)
                d = distance(point - projected_point, lattice, PBC=PBC)
                distances.append(np.min(d))
            #Choose wp with shortest translation for generating point
            tmpindex = np.argmin(distances)
            index = possible[tmpindex]
            newwp = group[index]
            projected_point = project_point(point, newwp[0], lattice=lattice, PBC=PBC)
            coor = apply_ops(projected_point, newwp)
            point = coor[0]
            index = newwp.index
        #Distances were not too small; return True
        else:
            return coor, index, point

def choose_wyckoff_molecular(group, number, orientations, general_site_only=True):
    """
    Choose a Wyckoff position to fill based on the current number of molecules
    needed to be placed within a unit cell

    Rules:

        1) The new position's multiplicity is equal/less than (number).
        2) We prefer positions with large multiplicity.
        3) The site must admit valid orientations for the desired molecule.

    Args:
        group: a pyxtal.symmetry.Group object
        number: the number of molecules still needed in the unit cell
        orientations: the valid orientations for a given molecule. Obtained
            from get_sg_orientations, which is called within molecular_crystal

    Returns:
        a single index for the Wyckoff position. If no position is found,
        returns False
    """
    wyckoffs = group.wyckoffs_organized
    
    if general_site_only or np.random.random()>0.5: #choose from high to low
        for j, wyckoff in enumerate(wyckoffs):
            if len(wyckoff[0]) <= number:
                good_wyckoff = []
                for k, w in enumerate(wyckoff):
                    if orientations[j][k] != []:
                        good_wyckoff.append(w)
                if len(good_wyckoff) > 0:
                    return random.choice(good_wyckoff)
        return False
    else:
        good_wyckoff = []
        for j, wyckoff in enumerate(wyckoffs):
            if len(wyckoff[0]) <= number:
                for k, w in enumerate(wyckoff):
                    if orientations[j][k] != []:
                        good_wyckoff.append(w)
        if len(good_wyckoff) > 0:
            return random.choice(good_wyckoff)
        else:
            return False

class mol_site():
    """
    Class for storing molecular Wyckoff positions and orientations within
    the molecular_crystal class. Each mol_site object represenents an
    entire Wyckoff position, not necessarily a single molecule. This is the
    molecular version of Wyckoff_site

    Args:
        mol: a Pymatgen Molecule object
        position: the fractional 3-vector representing the generating molecule's position
        orientation: an Orientation object for the generating molecule
        wyckoff_position: a Wyckoff_position object
        lattice: a Lattice object for the crystal
        ellipsoid: an optional binding Ellipsoid object for checking distances.
        tm: a Tol_matrix object for distance checking
    """
    def __init__(self, mol, position, orientation, wyckoff_position, lattice, ellipsoid=None, tm=Tol_matrix(prototype="molecular")):
        self.mol = mol
        """A Pymatgen molecule object"""
        self.position = wyckoff_position[0].operate(position)
        """Relative coordinates of the molecule's center within the unit cell"""
        self.orientation = orientation
        """The orientation object of the Mol in the first point in the WP"""
        self.ellipsoid = ellipsoid
        """A SymmOp representing the minimal ellipsoid for the molecule"""
        self.wp = wyckoff_position
        self.lattice = lattice
        """The crystal lattice in which the molecule resides"""
        self.multiplicity = self.wp.multiplicity
        """The multiplicity of the molecule's Wyckoff position"""
        self.PBC = wyckoff_position.PBC
        """The periodic axes"""
        self.tol_matrix = tm
        self.tols_matrix = self.get_tols_matrix()

    def __str__(self):
        s = str(self.mol.formula)+": "+str(self.position)+" "+str(self.wp.multiplicity)+self.wp.letter+", site symmetry "+ss_string_from_ops(self.wp.symmetry_m[0], self.wp.number, dim=self.wp.dim)
        phi, theta, psi = euler_from_matrix(self.orientation.matrix, radians=False)
        s += "\n    phi: "+str(phi)
        s += "\n    theta: "+str(theta)
        s += "\n    psi: "+str(psi)
        return s

    def get_tols_matrix(self):
        """
        Returns: a 2D matrix which is used internally for distance checking.
        """
        species = self.mol.species
        #Create tolerance matrix from subset of tm
        tm = self.tol_matrix
        tols = np.zeros((len(species),len(species)))
        for i1, specie1 in enumerate(species):
            for i2, specie2 in enumerate(species):
                tols[i1][i2] = tm.get_tol(specie1, specie2)
        return tols

    def get_ellipsoid(self):
        """
        Returns the bounding ellipsoid for the molecule. Applies the orientation
        transformation first.

        Returns:
            a re-orientated SymmOp representing the molecule's bounding ellipsoid
        """
        if self.ellipsoid == None:
            self.ellipsoid = find_ellipsoid(self.mol)
        e = self.ellipsoid
        #Appy orientation
        m = np.dot(e.rotation_matrix, self.orientation.get_matrix(angle=0))
        return SymmOp.from_rotation_and_translation(m, e.translation_vector)

    def get_ellipsoids(self):
        """
        Returns the bounding ellipsoids for the molecules in the WP. Includes the correct
        molecular centers and orientations.

        Returns:
            an array of re-orientated SymmOp's representing the molecule's bounding ellipsoids
        """
        #Get molecular centers
        centers0 = apply_ops(self.position, self.wp.generators)
        centers1 = np.dot(centers0, self.lattice)
        #Rotate ellipsoids
        e1 = self.get_ellipsoid()
        es = np.dot(self.wp.generators_m, e1)
        #Add centers to ellipsoids
        center_ops = [SymmOp.from_rotation_and_translation(Euclidean_lattice, c) for c in centers1]
        es_final = []
        for e, c in zip(es, center_ops):
            es_final.append(e*c)
        return np.array(es_final)

    def _get_coords_and_species(self, absolute=False, add_PBC=False, first=False):
        """
        Used to generate coords and species for get_coords_and_species

        Args:
            absolute: whether or not to return absolute (Euclidean)
                coordinates. If false, return relative coordinates instead
            add_PBC: whether or not to add coordinates in neighboring unit cells, used for
                distance checking
            first: whether or not to extract the information from only the first site
        
        Returns:
            atomic coords: a numpy array of fractional coordinates for the atoms in the site
            species: a list of atomic species for the atomic coords
        """
        mo = deepcopy(self.mol)
        mo.apply_operation(self.orientation.get_op(angle=0))
        wp_atomic_sites = []
        wp_atomic_coords = []
        for point_index, op2 in enumerate(self.wp.generators):
            current_atomic_sites = []
            current_atomic_coords = []

            #Rotate the molecule (Euclidean metric)
            op2_m = self.wp.generators_m[point_index]
            mo2 = deepcopy(mo)
            mo2.apply_operation(op2_m)
            #Obtain the center in absolute coords
            center_relative = op2.operate(self.position)
            center_absolute = np.dot(center_relative, self.lattice)
            #Add absolute center to molecule
            mo2.apply_operation(SymmOp.from_rotation_and_translation(np.identity(3),center_absolute))

            for site in mo2:
                #Place molecular coordinates in relative coordinates
                relative_coords = np.dot(site.coords, np.linalg.inv(self.lattice))
                #Do not filter: interferes with periodic image check
                #relative_coords = filtered_coords(relative_coords, PBC=self.PBC)
                current_atomic_sites.append(site.specie.name)
                current_atomic_coords.append(relative_coords)
            for s in current_atomic_sites:
                wp_atomic_sites.append(s)
            for c in current_atomic_coords:
                wp_atomic_coords.append(c)

            if first:
                break

        if add_PBC is True:
            #Filter PBC of wp_atomic_coords
            wp_atomic_coords = filtered_coords(wp_atomic_coords, PBC=self.PBC)
            #Add PBC copies of coords
            m = create_matrix(PBC=self.PBC)
            #Move [0,0,0] PBC vector to first position in array
            m2 = [[0,0,0]]
            v0 = np.array([0.,0.,0.])
            for v in m:
                if not (v==v0).all():
                    m2.append(v)
            new_coords = np.vstack([wp_atomic_coords + v for v in m2])
            wp_atomic_coords = new_coords

        return np.array(wp_atomic_coords), wp_atomic_sites

    def get_coords_and_species(self, absolute=False, add_PBC=False):
        """
        Lazily generates and returns the atomic coordinate and species for the
        Wyckoff position. Plugs the molecule into the provided orientation
        (with angle=0), and calculates the new positions.

        Args:
            absolute: whether or not to return absolute (Euclidean)
                coordinates. If false, return relative coordinates instead
            add_PBC: whether or not to add coordinates in neighboring unit cells, used for
                distance checking
        
        Returns:
            coords, species: coords is an np array of 3-vectors. species is
                a list of atomic species names, for example
                ['H', 'H', 'O', 'H', 'H', 'O']
        """
        if absolute is False:
            try:
                return self.relative_coords, self.species
            except:
                self.relative_coords, self.species = self._get_coords_and_species(absolute=absolute, add_PBC=add_PBC)
            return self.relative_coords, self.species
        else:
            try:
                return self.absolute_coords, self.species
            except:
                self.absolute_coords, self.species = self._get_coords_and_species(absolute=absolute, add_PBC=add_PBC)
            return self.absolute_coords, self.species

    def get_centers(self, absolute=False):
        """
        Returns the fractional coordinates for the center of mass for each molecule in
        the Wyckoff position

        Returns:
            A numpy array of fractional 3-vectors
        """
        centers = apply_ops(self.position, self.wp.generators)
        #centers1 = filtered_coords(centers0, self.PBC)
        if absolute is False:
            return centers
        else:
            return np.dot(centers, self.lattice)

    def get_radius(self):
        try:
            return self.radius
        except:
            r_max = 0
            for site in self.mol:
                radius = math.sqrt(site.x**2 + site.y**2 + site.z**2) + self.tol_matrix.get_tol(site.specie,site.specie)*0.5
                if radius > r_max:
                    r_max = radius
            self.radius = r_max
            return self.radius

    def check_distances(self, factor=1.0, atomic=True):
        """
        Checks if the atoms in the Wyckoff position are too close to each other
        or not. Does not check distances between atoms in the same molecule. Uses
        crystal.check_distance as the base code.
        
        Args:
            factor: the tolerance factor to use. A higher value means atoms must
                be farther apart
            atomic: if True, checks inter-atomic distances. If False, checks ellipsoid
                overlap between molecules instead
        
        Returns:
            True if the atoms are not too close together, False otherwise
        """
        if atomic is True:
            m_length = len(self.mol.species)
            #TODO: Use tm instead of tols lists
            #Get coords of WP with PBC
            coords, species = self._get_coords_and_species()

            #Get coords of the generating molecule
            coords_mol = coords[:m_length]
            #Remove generating molecule's coords from large array
            coords = coords[m_length:]

            if self.PBC != [0,0,0]:
                #Check periodic images
                m = create_matrix(PBC=self.PBC)
                #Remove original coordinates
                m2 = []
                v0 = np.array([0.,0.,0.])
                for v in m:
                    if not (v==v0).all():
                        m2.append(v)
                coords_PBC = np.vstack([coords_mol + v for v in m2])
                d = distance_matrix(coords_mol, coords_PBC, self.lattice, PBC=[0,0,0])
                tols = np.repeat(self.tols_matrix, len(m2), axis=1)
                if (d<tols).any():
                    return False

            if self.multiplicity > 1:
                #Check inter-atomic distances
                d = distance_matrix(coords_mol, coords, self.lattice, PBC=self.PBC)
                tols = np.repeat(self.tols_matrix, self.multiplicity-1, axis=1)
                if (d<tols).any():
                    return False

            return True



            """New method - only checks some atoms/molecules"""
            #Store length of molecule
            m_length = len(self.mol)
            #Get coordinates of center molecule and Wyckoff position
            coords, species = self._get_coords_and_species(absolute=True)
            coords_mol = coords[:m_length]

            if self.PBC == [0,0,0]:
                #Check non-periodic Wyckoff positions
                if self.multiplicity == 1:
                    return True
                coords_other = coords[m_length:]
                tols = np.repeat(self.tols_matrix, self.multiplicity-1, axis=1)
                d = cdist(coords_mol, coords_other)
                if (d<tols).any():
                    return False
                else:
                    return True

            #Create PBC vectors
            m = create_matrix(PBC=self.PBC)
            ml = np.dot(m, self.lattice)
            #Store the index of the (0,0,0) vector within ml
            mid_index = len(ml) // 2

            if self.multiplicity == 1:
                #Only check periodic images
                #Remove original coordinates
                m2 = []
                v0 = np.array([0.,0.,0.])
                for v in ml:
                    if not (v==v0).all():
                        m2.append(v)
                coords_PBC = np.vstack([coords_mol + v for v in m2])
                d = distance_matrix(coords_mol, coords_PBC, None, PBC=[0,0,0])
                tols = np.repeat(self.tols_matrix, len(m2), axis=1)
                if (d<tols).any():
                    return False
                else:
                    return True

            #Generate centers of all molecules
            centers = self.get_centers(absolute=True)
            vectors = np.repeat(centers, len(ml), axis=0) + np.tile(ml,(len(centers),1)) - np.dot(self.position, self.lattice)
            #Calculate distances between centers
            distances = np.linalg.norm(vectors, axis=-1)
            #Find which molecules need to be checked
            indices_mol = np.where(distances < self.get_radius()*2)[0]
            #Get indices of Wyckoff positions and PBC vectors
            indices_wp = []
            indices_pbc = []
            indices_vector = []
            for index in indices_mol:
                i_wp, i_pbc = divmod(index, len(ml))
                #Omit original center molecule
                if not (i_wp == 0 and i_pbc == mid_index):
                    indices_wp.append(i_wp)
                    indices_pbc.append(i_pbc)
                    indices_vector.append(index)

            if indices_wp == []:
                return True

            #Get atomic positions of molecules with small separation vectors
            original_coords = np.vstack([coords[index_wp*m_length:index_wp*m_length+m_length] for index_wp in indices_wp])
            pbc_toadd = np.repeat(ml[indices_pbc], m_length, axis=0)
            atomic_coords = original_coords + pbc_toadd
            #Get inter-atomic tolerances
            tols = np.tile(self.get_tols_matrix(), len(indices_wp))
            if m_length <= max_fast_mol_size:
                #Check all atomic pairs
                d = cdist(coords_mol, atomic_coords)

                """
                print("~~~~~~~~~~~~~~~~~~~~~~~")
                print("ml:", ml.shape)
                print(ml)
                print("centers:", centers.shape)
                print(centers)
                print("vectors:", vectors.shape)
                print(vectors)
                print("radius*2: ", self.get_radius()*2)
                print("distances:", distances.shape)
                print(distances)
                print("indices_mol:", len(indices_mol))
                print(indices_mol)
                print("indices_wp:", len(indices_wp))
                print(indices_wp)
                print("indices_pbc:", len(indices_pbc))
                print(indices_pbc)
                print("indices_vector:", len(indices_vector))
                print(indices_vector)
                print("coords_mol:", coords_mol.shape)
                print(coords_mol)
                print("coords:", coords.shape)
                print(coords)
                print("original_coords:", original_coords.shape)
                print(original_coords)
                print("pbc_toadd:", pbc_toadd.shape)
                print(pbc_toadd[:12])
                print("atomic_coords: ", atomic_coords.shape)
                print(atomic_coords)
                print("d:", d.shape)
                print(d)
                print("tols_matrix:", self.get_tols_matrix().shape)
                print(self.get_tols_matrix())
                print("tols:", tols.shape)
                print(tols)
                """


                if (d<tols).any():
                    return False
                else:
                    return True

            elif m_length > max_fast_mol_size:
                #Get corresponding separation vectors
                new_vectors = np.repeat(vectors[indices_vector], m_length, axis=0)
                #Get atomic coordinates relative to molecular centers
                relative_atomic_coords = atomic_coords - new_vectors
                #Dot atomic coordinates with inter-molecular separation vectors
                dots = np.einsum('...j,...j', new_vectors, relative_atomic_coords)
                #Find where relative vectors point towards the original molecule
                new_indices = np.where(dots<0)[0]
                #Get new coordinates and tolerances for distance matrix
                new_atomic_coords = atomic_coords[new_indices]
                d = cdist(coords_mol, new_atomic_coords)
                tols2 = tols[:,new_indices]
                if (d<tols2).any():
                    return False
                else:
                    return True    
            


        elif atomic is False:
            #Check molecular ellipsoid overlap
            if self.multiplicity == 1:
                return True
            es0 = self.get_ellipsoids()[1:]
            PBC_vectors = np.dot(create_matrix(PBC=self.PBC), self.lattice)
            PBC_ops = [SymmOp.from_rotation_and_translation(Euclidean_lattice, v) for v in PBC_vectors]
            es1 = []
            for op in PBC_ops:
                es1.append(np.dot(es0, op))
            es1 = np.squeeze(es1)
            truth_values = np.vectorize(check_intersection)(es1, self.get_ellipsoid())
            if np.sum(truth_values) < len(truth_values):
                return False
            else:
                return True

class molecular_crystal():
    """
    Class for storing and generating molecular crystals based on symmetry
    constraints. Based on the crystal.random_crystal class for atomic crystals.
    Given a spacegroup, list of molecule objects, molecular stoichiometry, and
    a volume factor, generates a molecular crystal consistent with the given
    constraints. This crystal is stored as a pymatgen struct via self.struct
    
    Args:
        group: The international spacegroup number
            OR, a pyxtal.symmetry.Group object
        molecules: a list of pymatgen.core.structure.Molecule objects for
            each type of molecule. Alternatively, you may supply a file path,
            or give a string to convert (in which case fmt must be defined)
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
        fmt: Optional value for the input molecule string format. Used only
            when molecule values are strings
        lattice: an optional Lattice object to use for the unit cell
        tm: the Tol_matrix object used to generate the crystal
    """

    def init_common(self, molecules, numMols, volume_factor, select_high, allow_inversion, orientations, check_atomic_distances, group, lattice, tm):
        """
        init functionality which is shared by 3D, 2D, and 1D crystals
        """
        self.numattempts = 0
        """The number of attempts needed to generate the crystal."""
        if type(group) == Group:
            self.group = group
            """A pyxtal.symmetry.Group object storing information about the space/layer
            /Rod/point group, and its Wyckoff positions."""
        else:
            self.group = Group(group, dim=self.dim)
        self.number = self.group.number
        """The international group number of the crystal:
        1-230 for 3D space groups
        1-80 for 2D layer groups
        1-75 for 1D Rod groups
        1-32 for crystallographic point groups
        None otherwise
        """
        self.Msgs()
        """A list of warning messages to use during generation."""
        self.factor = volume_factor
        """The supplied volume factor for the unit cell."""
        numMols = np.array(numMols) #must convert it to np.array
        self.numMols0 = numMols
        """The number of each type of molecule in the PRIMITIVE cell"""
        self.numMols = self.numMols0 * cellsize(self.group)
        """The number of each type of molecule in the CONVENTIONAL cell"""
        oriented_molecules = []
        #Allow support for generating molecules from text via openbable
        for i, mol in enumerate(molecules):
            # Parsing the molecular input
            if type(mol) == str:
                # Read strings into molecules, try collection first,
                # If string not in collection, use pymatgen format
                tmp = mol.split('.')
                if len(tmp) > 1:
                    print('\nLoad the molecule from the given file {:s}'.format(mol))
                    if tmp[-1] in ['xyz', 'gjf', 'g03', 'json']:
                        mo = mol_from_file(mol)
                    else:
                        raise NameError('{:s} is not a supported format'.format(tmp[-1]))
                else:
                    print('\nLoad the molecule {:s} from the built_in collections'.format(mol))
                    mo = molecule_collection[mol]
                if mo is not None:
                    molecules[i] = mo
                else:
                    raise NameError("Could not create molecules from given input: {:s}".format(mol))
        for mol in molecules:
            if len(mol) > 1:
                props = mol.site_properties
                pga = PointGroupAnalyzer(mol)
                mo = pga.symmetrize_molecule()['sym_mol']
                if len(props)>0:
                    for key in props.keys():
                        mo.add_site_property(key, props[key])
                #print(mo.site_properties)
                #import sys
                #sys.exit()
                oriented_molecules.append(mo)
            else:
                oriented_molecules.append(mol)
        self.molecules = oriented_molecules
        """A list of pymatgen.core.structure.Molecule objects, symmetrized and
        oriented along their symmetry axes."""
        self.boxes = []
        """A list of bounding boxes for each molecule. Used for estimating
        volume of the unit cell."""
        self.radii = []
        """A list of approximated radii for each molecule type. Used for
        checking inter-molecular distances."""
        #Calculate boxes and radii for each molecule
        for mol in self.molecules:
            self.boxes.append(get_box(reoriented_molecule(mol)[0]))
            max_r = 0
            for site in mol:
                radius = math.sqrt( site.x**2 + site.y**2 + site.z**2 )
                if radius > max_r: max_r = radius
            self.radii.append(max_r+1.0)
        """The volume of the generated unit cell"""
        self.check_atomic_distances = check_atomic_distances
        """Whether or not inter-atomic distances are checked at each step."""
        self.allow_inversion = allow_inversion
        """Whether or not to allow chiral molecules to be inverted."""
        self.select_high = select_high
        """Whether or not to select high multi Wycoff sites only."""
        #When generating multiple crystals of the same stoichiometry and sg,
        #allow the user to re-use the allowed orientations, to reduce time cost
        if orientations is None:
            self.get_orientations()
        else:
            self.valid_orientations = orientations
            """The valid orientations for each molecule and Wyckoff position.
            May be copied when generating a new molecular_crystal to save a
            small amount of time"""
        if lattice is not None:
            #Use the provided lattice
            self.lattice = lattice
            self.volume = lattice.volume
            #Make sure the custom lattice PBC axes are correct.
            if lattice.PBC != self.PBC:
                self.lattice.PBC = self.PBC
                printx("\n  Warning: converting custom lattice PBC to "+str(self.PBC))
        elif lattice is None:
            #Determine the unique axis
            if self.dim == 2:
                if self.number in range(3, 8):
                    unique_axis = "c"
                else:
                    unique_axis = "a"
            elif self.dim == 1:
                if self.number in range(3, 8):
                    unique_axis = "a"
                else:
                    unique_axis = "c"
            else:
                unique_axis = "c"
            #Generate a Lattice instance
            self.volume = estimate_volume_molecular(self.molecules, self.numMols, self.factor, boxes=self.boxes)
            """The volume of the generated unit cell."""

            #Calculate the minimum, middle, and maximum box lengths for the unit cell.
            #Used to make sure at least one non-overlapping orientation exists for each molecule
            minls = []
            midls = []
            maxls = []
            for box in self.boxes:
                minls.append(box.minl)
                midls.append(box.midl)
                maxls.append(box.maxl)

            if self.dim == 3 or self.dim == 0:
                self.lattice = Lattice(self.group.lattice_type, self.volume, PBC=self.PBC, unique_axis=unique_axis, min_l=max(minls), mid_l=max(midls), max_l=max(maxls))
            elif self.dim == 2:
                self.lattice = Lattice(self.group.lattice_type, self.volume, PBC=self.PBC, unique_axis=unique_axis, min_l=max(minls), mid_l=max(midls), max_l=max(maxls), thickness=self.thickness)
            elif self.dim == 1:
                self.lattice = Lattice(self.group.lattice_type, self.volume, PBC=self.PBC, unique_axis=unique_axis, min_l=max(minls), mid_l=max(midls), max_l=max(maxls), area=self.area)
            """The Lattice object used to generate lattice matrices for the structure."""
        #Set the tolerance matrix
        if type(tm) == Tol_matrix:
            self.tol_matrix = tm
            """The Tol_matrix object used for checking inter-atomic distances within the structure."""
        else:
            try:
                self.tol_matrix = Tol_matrix(prototype=tm)
            except:
                printx("Error: tm must either be a Tol_matrix object or a prototype string for initializing one.", priority=1)
                self.valid = False
                self.struct = None
                return
        self.generate_crystal()

    def __init__(self, group, molecules, numMols, volume_factor, select_high=True, allow_inversion=False, orientations=None, check_atomic_distances=True, fmt="xyz", lattice=None, tm=Tol_matrix(prototype="molecular")):
        self.dim = 3
        """The number of periodic dimensions of the crystal"""
        #Necessary input
        self.PBC = [1,1,1]
        """The periodic axes of the crystal"""
        if type(group) != Group:
            group = Group(group, self.dim)
        self.sg = group.number
        self.selec_high = select_high
        """The international spacegroup number of the crystal."""
        self.init_common(molecules, numMols, volume_factor, select_high, allow_inversion, orientations, check_atomic_distances, group, lattice, tm)

    def Msgs(self):
        self.Msg1 = 'Error: the stoichiometry is incompatible with the wyckoff sites choice'
        self.Msg2 = 'Error: failed in the cycle of generating structures'
        self.Msg3 = 'Warning: failed in the cycle of adding species'
        self.Msg4 = 'Warning: failed in the cycle of choosing wyckoff sites'
        self.Msg5 = 'Finishing: added the specie'
        self.Msg6 = 'Finishing: added the whole structure'
        self.Msg7 = 'Error: invalid paramaters for initialization'

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
        self.group.wyckoffs_organized[j][k]
        """
        self.valid_orientations = []
        for mol in self.molecules:
            self.valid_orientations.append([])
            wp_index = -1
            for i, x in enumerate(self.group.wyckoffs_organized):
                self.valid_orientations[-1].append([])
                for j, wp in enumerate(x):
                    wp_index += 1
                    allowed = orientation_in_wyckoff_position(mol, wp, already_oriented=True, allow_inversion=self.allow_inversion)
                    if allowed is not False:
                        self.valid_orientations[-1][-1].append(allowed)
                    else:
                        self.valid_orientations[-1][-1].append([])

    def check_compatible(group, numMols, valid_orientations):
        """
        Checks if the number of molecules is compatible with the Wyckoff
        positions. Considers the number of degrees of freedom for each Wyckoff
        position, and makes sure at least one valid combination of WP's exists.
        """
        #Store whether or not at least one degree of freedom exists
        has_freedom = False
        #Store the wp's already used that don't have any freedom
        used_indices = []
        #Loop over species
        for i_mol, numIon in enumerate(numMols):
            #Get lists of multiplicity, maxn and freedom
            l_mult0 = []
            l_maxn0 = []
            l_free0 = []
            indices0 = []
            for i_wp, wp in enumerate(group):
                #Check that at least one valid orientation exists
                j, k = jk_from_i(i_wp, group.wyckoffs_organized)
                if len(valid_orientations[i_mol][j][k]) > 0:
                    indices0.append(i_wp)
                    l_mult0.append(len(wp))
                    l_maxn0.append(numIon // len(wp))
                    if np.allclose(wp[0].rotation_matrix, np.zeros([3,3])):
                        l_free0.append(False)
                    else:
                        l_free0.append(True)
            #Remove redundant multiplicities:
            l_mult = []
            l_maxn = []
            l_free = []
            indices = []
            for mult, maxn, free, i_wp in zip(l_mult0, l_maxn0, l_free0, indices0):
                if free is True:
                    if mult not in l_mult:
                        l_mult.append(mult)
                        l_maxn.append(maxn)
                        l_free.append(True)
                        indices.append(i_wp)
                elif free is False and i_wp not in used_indices:
                    l_mult.append(mult)
                    l_maxn.append(1)
                    l_free.append(False)
                    indices.append(i_wp)
            #Loop over possible combinations
            #Create pointer variable to move through lists
            p = 0
            #Store the number of each WP, used across possible WP combinations
            n0 = [0]*len(l_mult)
            n = deepcopy(n0)
            for i, mult in enumerate(l_mult):
                if l_maxn[i] != 0:
                    p = i
                    n[i] = l_maxn[i]
                    break
            p2 = p
            if n == n0:
                return False
            while True:
                num = np.dot(n, l_mult)
                dobackwards = False
                #The combination works: move to next species
                if num == numIon:
                    #Check if at least one degree of freedom exists
                    for val, free, i_wp in zip(n, l_free, indices):
                        if val > 0:
                            if free is True:
                                has_freedom = True
                            elif free is False:
                                indices.append(i_wp)
                    break
                #All combinations failed: return False
                if n == n0 and p >= len(l_mult) - 1:
                    return False
                #Too few atoms
                if num < numIon:
                    #Forwards routine
                    #Move p to the right and max out
                    if p < len(l_mult) - 1:
                        p += 1
                        n[p] = min((numIon - num) // l_mult[p], l_maxn[p])
                    else:
                        #p is already at last position: trigger backwards routine
                        dobackwards = True
                #Too many atoms
                if num > numIon or dobackwards is True:
                    #Backwards routine
                    #Set n[p] to 0, move p backwards to non-zero, and decrease by 1
                    n[p] = 0
                    while p > 0 and p > p2:
                        p -= 1
                        if n[p] != 0:
                            n[p] -= 1
                            if n[p] == 0 and p == p2:
                                p2 = p + 1
                            break
        #All species passed: return True
        if has_freedom is True:
            return True
        #All species passed, but no degrees of freedom: return 0
        elif has_freedom is False:
            return 0

    def to_file(self, fmt=None, filename=None):
        """
        Creates a file with the given filename and file type to store the structure.
        By default, creates cif files for crystals and xyz files for clusters.
        By default, the filename is based on the stoichiometry.

        Args:
            fmt: the file type ('cif', 'xyz', etc.)
            filename: the file path

        Returns:
            Nothing. Creates a file at the specified path
        """
        if filename == None:
            given = False
        else:
            given = True
        if self.valid:
            if fmt == None:
                fmt = "cif"
            if filename == None:
                filename = str(self.struct.formula).replace(" ","") + "." + fmt
            #Check if filename already exists
            #If it does, add a new number to end of filename
            if os.path.exists(filename):
                if given is False:
                    filename = filename[:(-len(fmt)-1)]
                i = 1
                while True:
                    outdir = filename + "_" + str(i)
                    if given is False:
                        outdir += "." + fmt
                    if not os.path.exists(outdir):
                        break
                    i += 1
                    if i > 10000:
                        return "Could not create file: too many files already created."
            else:
                outdir = filename
            self.struct.to(fmt=fmt, filename=outdir)
            return "Output file to " + outdir
        elif self.valid:
            printx("Cannot create file: structure did not generate.", priority=1)

    def print_all(self):
        print("--Molecular Crystal--")
        print("Dimension: "+str(self.dim))
        print("Group: "+self.group.symbol)
        print("Volume factor: "+str(self.factor))
        if self.valid:
            print("Wyckoff sites:")
            for x in self.mol_generators:
                print("  "+str(x))
            print("Pymatgen Structure:")
            print(self.struct)

    def generate_crystal(self, max1=max1, max2=max2, max3=max3, max4=max4):
        """
        The main code to generate a random molecular crystal. If successful,
        stores a pymatgen.core.structure object in self.struct and sets
        self.valid to True. If unsuccessful, sets self.valid to False and
        outputs an error message.

        Args:
            max1: the number of attempts for generating a lattice
            max2: the number of attempts for a given lattice
            max3: the number of attempts for a given Wyckoff position
            max4: the number of attempts for changing the molecular orientation
        """
        #Check the minimum number of degrees of freedom within the Wyckoff positions
        degrees = molecular_crystal.check_compatible(self.group, self.numMols, self.valid_orientations)
        if degrees is False:
            self.struct = None
            self.valid = False
            return
        else:
            if degrees == 0:
                max1 = 20
                max2 = 3
                max3 = 1
                max4 = 1
            #Calculate a minimum vector length for generating a lattice
            all_lengths = []
            for box in self.boxes:
                all_lengths.append(box.minl)
            minvector = max(all_lengths) * 1.2 #Require slightly larger min vector for lattices
            for cycle1 in range(max1):
                self.cycle1 = cycle1
                #1, Generate a lattice
                if self.lattice.allow_volume_reset is True:
                    self.volume = estimate_volume_molecular(self.molecules, self.numMols, self.factor, boxes=self.boxes)
                    self.lattice.volume = self.volume
                self.lattice.reset_matrix()
                try:
                    cell_matrix = self.lattice.get_matrix
                except:
                    continue
                cell_para = self.lattice.get_para()

                if cell_para is None:
                    break
                else:
                    cell_matrix = para2matrix(cell_para)
                    if abs(self.volume - np.linalg.det(cell_matrix)) > 1.0: 
                        printx('Error, volume is not equal to the estimated value: ', self.volume, ' -> ', np.linalg.det(cell_matrix), priority=0)
                        printx('cell_para:  ', cell_para, priority=0)
                        sys.exit(0)

                    for cycle2 in range(max2):
                        molecular_coordinates_total = [] #to store the added molecular coordinates
                        molecular_sites_total = []      #to store the corresponding molecular specie
                        coordinates_total = [] #to store the added atomic coordinates
                        species_total = []      #to store the corresponding atomic specie
                        wps_total = []      #to store corresponding Wyckoff position indices
                        points_total = []   #to store the generating x,y,z points
                        mol_generators_total = []
                        good_structure = False


                        self.cycle2 = cycle2
                        molecular_coordinates_tmp = deepcopy(molecular_coordinates_total)
                        molecular_sites_tmp = deepcopy(molecular_sites_total)
                        coordinates_tmp = deepcopy(coordinates_total)
                        species_tmp = deepcopy(species_total)
                        wps_tmp = deepcopy(wps_total)
                        points_tmp = deepcopy(points_total)
                        mol_generators_tmp = []
                        
                        #Add molecules specie by specie
                        for i, numMol in enumerate(self.numMols):
                            mol = self.molecules[i]
                            numMol_added = 0

                            #Now we start to add the specie to the wyckoff position
                            for cycle3 in range(max3):
                                self.cycle3 = cycle3
                                self.cycle4 = 0
                                self.numattempts += 1
                                #Choose a random Wyckoff position for given multiplicity: 2a, 2b, 2c
                                #NOTE: The molecular version return wyckoff indices, not ops
                                wp = choose_wyckoff_molecular(self.group, numMol-numMol_added, self.valid_orientations[i], self.select_high)
                                if wp is not False:
                                    #Generate a list of coords from the wyckoff position
                                    point = self.lattice.generate_point()
                                    projected_point = project_point(point, wp[0], lattice=cell_matrix, PBC=self.PBC)
                                    coords = apply_ops(projected_point, wp)
                                    #merge coordinates if the atoms are close
                                    if self.check_atomic_distances is False:
                                        mtol = self.radii[i]*2
                                    elif self.check_atomic_distances is True:
                                        mtol = self.radii[i]*0.5
                                    coords_toadd, good_merge, point = merge_coordinate_molecular(coords, cell_matrix, self.group, mtol, self.valid_orientations[i])
                                    if good_merge is not False:
                                        wp_index = good_merge
                                        coords_toadd = filtered_coords(coords_toadd, PBC=self.PBC) #scale the coordinates to [0,1], very important!

                                        #Create a mol_site object
                                        mo = deepcopy(self.molecules[i])
                                        j, k = jk_from_i(wp_index, self.group.wyckoffs_organized)
                                        ori = random.choice(self.valid_orientations[i][j][k]).random_orientation()
                                        ms0 = mol_site(mo, point, ori, self.group[wp_index], cell_matrix, tm=self.tol_matrix)
                                        #Check distances within the WP
                                        if ms0.check_distances(atomic=self.check_atomic_distances) is False: #continue
                                            #Check distance between centers
                                            d = distance_matrix(ms0.get_centers(), ms0.get_centers(), ms0.lattice, PBC=ms0.PBC)
                                            min_box_l = self.boxes[i].minl
                                            xys = np.where(d < min_box_l)
                                            passed_center = True
                                            for i_y, x in enumerate(xys[0]):
                                                y = xys[1][i_y]
                                                val = d[x][y]
                                                #Ignore self-distances
                                                if x == y:
                                                    continue
                                                else:
                                                    passed_center = False
                                            if not passed_center: continue
                                            #If centers are farther apart than min box length, allow multiple orientation attempts
                                            passed_ori = False
                                            for cycle4 in range(max4):
                                                self.cycle4 = cycle4
                                                ori = random.choice(self.valid_orientations[i][j][k]).random_orientation()
                                                ms0 = mol_site(mo, point, ori, self.group[wp_index], cell_matrix, tm=self.tol_matrix)
                                                if ms0.check_distances(atomic=self.check_atomic_distances):
                                                    passed_ori = True
                                                    break
                                        else:
                                            passed_ori = True
                                        if passed_ori is False: continue
                                        #Check distances with other WP's
                                        coords_toadd, species_toadd = ms0.get_coords_and_species()
                                        passed = True
                                        for ms1 in mol_generators_tmp:
                                            if check_mol_sites(ms0, ms1, atomic=self.check_atomic_distances, tm=self.tol_matrix) is False:
                                                passed = False
                                                break
                                        if passed is False: continue
                                        elif passed is True:
                                            #Distance checks passed; store the new Wyckoff position
                                            mol_generators_tmp.append(ms0)
                                            if coordinates_tmp == []:
                                                coordinates_tmp = coords_toadd
                                            else:
                                                coordinates_tmp = np.vstack([coordinates_tmp, coords_toadd])
                                            species_tmp += species_toadd
                                            numMol_added += len(coords_toadd)/len(mo)
                                            if numMol_added == numMol:
                                                #We have enough molecules of the current type
                                                mol_generators_total = deepcopy(mol_generators_tmp)
                                                coordinates_total = deepcopy(coordinates_tmp)
                                                species_total = deepcopy(species_tmp)
                                                break

                            if numMol_added != numMol:
                                break  #need to repeat from the 1st species

                        if numMol_added == numMol:
                            printx(self.Msg6, priority=3)
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

                        final_coor = deepcopy(coordinates_total)
                        final_site = deepcopy(species_total)
                        final_number = list(Element(ele).z for ele in species_total)
                        self.mol_generators = deepcopy(mol_generators_total)
                        """A list of mol_site objects which can be used
                        for generating the crystal."""

                        final_coor = filtered_coords(final_coor, PBC=self.PBC)
                        final_lattice, final_coor = Add_vacuum(final_lattice, final_coor, PBC=self.PBC)
                        #if verify_distances(final_coor, final_site, final_lattice, factor=0.75, PBC=self.PBC):
                        self.lattice_matrix = final_lattice
                        """A 3x3 matrix representing the lattice of the
                        unit cell."""
                        self.frac_coords = np.array(final_coor)
                        """The fractional coordinates for each molecule
                        in the final structure"""
                        self.cart_coords = np.dot(final_coor, final_lattice)
                        """The absolute coordinates for each molecule
                        in the final structure"""
                        self.sites = final_site
                        """The indices within self.molecules corresponding
                        to the type of molecule for each site in
                        self.cart_coords."""              
                        self.struct = Structure(final_lattice, self.sites, self.frac_coords)
                        """A pymatgen.core.structure.Structure object for
                        the final generated crystal."""
                        self.spg_struct = (final_lattice, self.frac_coords, final_number)
                        """A list of information describing the generated
                        crystal, which may be used by spglib for symmetry
                        analysis."""
                        self.valid = True
                        """Whether or not a valid crystal was generated."""
                        return
                    else: printx("Failed final distance check.", priority=3)
        printx("Couldn't generate crystal after max attempts.", priority=1)
        if degrees == 0:
            printx("Note: Wyckoff positions have no degrees of freedom.", priority=2)
        self.struct = None
        self.valid = False
        return

class molecular_crystal_2D(molecular_crystal):
    """
    A 2d counterpart to molecular_crystal. Given a layer group, list of
    molecule objects, molecular stoichiometry, and
    a volume factor, generates a molecular crystal consistent with the given
    constraints. This crystal is stored as a pymatgen struct via self.struct
    
    Args:
        group: the layer group number between 1 and 80. NOT equal to the
            international space group number, which is between 1 and 230
            OR, a pyxtal.symmetry.Group object
        molecules: a list of pymatgen.core.structure.Molecule objects for
            each type of molecule. Alternatively, you may supply a file path,
            or give a string to convert (in which case fmt must be defined)
        numMols: A list of the number of each type of molecule within the
            primitive cell (NOT the conventioal cell)
        thickness: the thickness, in Angstroms, of the unit cell in the 3rd
            dimension (the direction which is not repeated periodically). A
            value of None causes a thickness to be chosen automatically. Note
            that this constraint applies only to the molecular centers; some
            atomic coordinates may lie outside of this range
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
        fmt: Optional value for the input molecule string format. Used only
            when molecule values are strings
        lattice: an optional Lattice object to use for the unit cell
        tm: the Tol_matrix object used to generate the crystal
    """
    def __init__(self, group, molecules, numMols, volume_factor, select_high=True, allow_inversion=False, orientations=None, check_atomic_distances=True, fmt='xyz', thickness=None, lattice=None, tm=Tol_matrix(prototype="molecular")):
        self.dim = 2
        """The number of periodic dimensions of the crystal"""
        self.numattempts = 0
        """The number of attempts needed to generate the crystal."""
        #Necessary input
        if type(group) != Group:
            group = Group(group, self.dim)
        number = group.number
        """The layer group number of the crystal."""
        self.thickness = thickness
        """the thickness, in Angstroms, of the unit cell in the 3rd
        dimension."""
        self.PBC = [1,1,0]
        """The periodic axes of the crystal."""
        self.init_common(molecules, numMols, volume_factor, select_high, allow_inversion, orientations, check_atomic_distances, group, lattice, tm)

class molecular_crystal_1D(molecular_crystal):
    """
    A 1d counterpart to molecular_crystal. Given a Rod group, list of
    molecule objects, molecular stoichiometry, volume factor, and area,
    generates a molecular crystal consistent with the given constraints.
    The crystal is stored as a pymatgen struct via self.struct
    
    Args:
        group: the Rod group number between 1 and 80. NOT equal to the
            international space group number, which is between 1 and 230
            OR, a pyxtal.symmetry.Group object
        molecules: a list of pymatgen.core.structure.Molecule objects for
            each type of molecule. Alternatively, you may supply a file path,
            or give a string to convert (in which case fmt must be defined)
        numMols: A list of the number of each type of molecule within the
            primitive cell (NOT the conventioal cell)
        area: cross-sectional area of the unit cell in Angstroms squared. A
            value of None causes an area to be chosen automatically. Note that
            this constraint applies only to the molecular centers; some atomic
            coordinates may lie outside of this range
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
        fmt: Optional value for the input molecule string format. Used only
            when molecule values are strings
        lattice: an optional Lattice object to use for the unit cell
        tm: the Tol_matrix object used to generate the crystal
    """
    def __init__(self, group, molecules, numMols, volume_factor, select_high=True, allow_inversion=False, orientations=None, check_atomic_distances=True, fmt='xyz', area=None, lattice=None, tm=Tol_matrix(prototype="molecular")):
        self.dim = 1
        """The number of periodic dimensions of the crystal"""
        #Necessary input
        self.area = area
        """the effective cross-sectional area, in Angstroms squared, of the
        unit cell."""
        self.PBC = [0,0,1]
        """The periodic axes of the crystal (1,2,3)->(x,y,z)."""
        self.sg = None
        """The international space group number (there is not a 1-1 correspondence
        with Rod groups)."""
        self.init_common(molecules, numMols, volume_factor, select_high, allow_inversion, orientations, check_atomic_distances, group, lattice, tm)


if __name__ == "__main__":
    #-------------------------------- Options -------------------------
    import os
    from time import time
    from spglib import get_symmetry_dataset

    parser = OptionParser()
    parser.add_option("-s", "--spacegroup", dest="sg", metavar='sg', default=36, type=int,
            help="desired space group number: 1-230, e.g., 36")
    parser.add_option("-e", "--molecule", dest="molecule", default='H2O', 
            help="desired molecules: e.g., H2O", metavar="molecule")
    parser.add_option("-n", "--numMols", dest="numMols", default=4, 
            help="desired numbers of molecules: 4", metavar="numMols")
    parser.add_option("-f", "--factor", dest="factor", default=1.0, type=float, 
            help="volume factor: default 1.0", metavar="factor")
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
            help="desired dimension: (3 or 2 for 3d or 2d, respectively): default 3")
    parser.add_option("-t", "--thickness", dest="thickness", metavar='thickness', default=None, type=float,
            help="Thickness, in Angstroms, of a 2D crystal, or area of a 1D crystal, None generates a value automatically: default None")

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
            system.append(mol)
        for x in number.split(','):
            numMols.append(int(x))
    else:
        system = [molecule]
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
