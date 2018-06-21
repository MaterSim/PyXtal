from structure import *
from molecule import *

def estimate_volume_molecular(numMols, boxes, factor=2.0):
    '''
    Estimate the volume needed for a molecular crystal unit cell.
    args:
        numMols: A list with the number of each type of molecule
        boxes: A list of bounding boxes for each molecule. Obtained from get_box
    '''
    volume = 0
    for numMol, box in zip(numMols, boxes):
        volume += (box[1]-box[0])*(box[3]-box[2])*(box[5]-box[4])
    return factor*volume

def get_sg_orientations(mol, sg, allow_inversion=False):
    """
    Calculate the valid orientations for each Molecule and Wyckoff position.
    Returns a list with 3 indices:
    index 1: the Wyckoff position's 1st index (based on multiplicity)
    index 2: the WP's 2nd index (within the group of equal multiplicity)
    index 3: the index of the valid orientation for the molecule/WP pair
    For example, self.valid_orientations[i][j] would be a list of valid
        orientations for self.molecules[i],
        in the Wyckoff position self.wyckoffs[i][j]
    """
    valid_orientations = []
    wyckoffs = get_wyckoffs(sg, organized=True)
    wp_index = -1
    for i, x in enumerate(wyckoffs):
        valid_orientations.append([])
        for j, wp in enumerate(x):
            wp_index += 1
            allowed = orientation_in_wyckoff_position(mol, sg, wp_index, already_oriented=True, allow_inversion=allow_inversion)
            if allowed is not False:
                valid_orientations[-1].append(allowed)
            else:
                valid_orientations[-1].append([])
    return valid_orientations

def get_box(mol, padding=1.5):
    '''
    Given a molecule, find a minimum orthorhombic box containing it.
    Size is calculated using min and max x, y, and z values.
    Returns a list [x1,x2,y1,y2,z1,z2] where x1 is the relative displacement in
    the negative x direction, x2 is the displacement in the positive x
    direction, and so on. For best results, call oriented_molecule first.
    args:
        mol: a pymatgen molecule object
        padding: the extra space to be added in each direction. Double this
            amount will be added to each of the x, y, and z directions.
    '''
    minx, miny, minz, maxx, maxy, maxz = 0.,0.,0.,0.,0.,0.
    for i, p in enumerate(mol):
        x, y, z = p.coords
        if x < minx: minx = x
        if y < minx: minx = y
        if z < minx: minx = z
        if x > minx: maxx = x
        if y > minx: maxx = y
        if z > minx: maxx = z
    return [minx-padding,maxx+padding,miny-padding,maxy+padding,minz-padding,minz+padding]

def check_distance_molecular(coord1, coord2, indices1, index2, lattice, radii):
    #NOTE: Currently does not depend on molecular orientations
    """
    check the distances between two set of molecules
    Args:
    coord1: multiple list of positions e.g. [[0,0,0],[1,1,1]]
    indices1: the corresponding molecular indices of coord1, e.g. [1, 3]
    coord2: a list of new positions: [0.5, 0.5 0.5]
    index2: the molecular index for coord2: 4
    lattice: cell matrix
    """
    #add PBC
    coord2s = []
    matrix = create_matrix()
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

            #print(d_min, tol)
            if d_min < tol:
                return False
        return True
    else:
        return True

def check_wyckoff_position_molecular(points, sg, orientations, wyckoffs=None, exact_translation=False):
    '''
    Given a list of points, return index of Wyckoff position in space group.
    If no match found, returns False.

    Args:
        points: a list of 3d coordinates or SymmOps to check
        sg: the international space group number to check
        wyckoffs: a list of Wyckoff positions obtained from get_wyckoffs.
        exact_translation: whether we require two SymmOps to have exactly equal
            translational components. If false, translations related by +-1
            are considered equal
    '''
    points = np.array(points)
    points = np.around((points*1e+10))/1e+10

    if wyckoffs == None:
        wyckoffs = get_wyckoffs(sg)
        gen_pos = wyckoffs[0]
    else:
        gen_pos = wyckoffs[0][0]
    new_points = []
    #
    if exact_translation == False:
        for p in points:
            new_points.append(p - np.floor(p))
        points = new_points
    w_symm_all = get_wyckoff_symmetry(sg)
    p_symm = []
    #If exact_translation is false, store WP's which might be a match
    possible = []
    for x in points:
        p_symm.append(site_symm(x, gen_pos))
    for i, wp in enumerate(wyckoffs):
        w_symm = w_symm_all[i]
        if len(p_symm) == len(w_symm):
            temp = deepcopy(w_symm)
            for p in p_symm:
                for w in temp:
                    if exact_translation:
                        if p == w:
                            temp.remove(w)
                    elif not exact_translation:
                        temp2 = deepcopy(w)
                        for op_p in p:
                            for op_w in w:
                                #Check that SymmOp's are equal up to some integer translation
                                if are_equal(op_w, op_p, allow_pbc=True):
                                    temp2.remove(op_w)
                        if temp2 == []:
                            temp.remove(w)
            if temp == []:
                #If we find a match with exact translations
                if exact_translation:
                    #Check orientations
                    j, k = jk_from_i(i, orientations)
                    if orientations[j][k] != []:
                        return i
                    else:
                        return False
                elif not exact_translation:
                    possible.append(i)
        #If no matching WP's are found
    if len(possible) == 0:
        return False
    #If exactly one matching WP is found
    elif len(possible) == 1:
        j, k = jk_from_i(possible[0], orientations)
        if orientations[j][k] != []:
            return possible[0]
        else:
            return False
    #If multiple WP's are found
    elif len(possible) > 1:
        #TODO: add a way to differentiate between possible WP's
        #print("Warning: multiple Wyckoff positions found")
        new = []
        for i in possible:
            j, k = jk_from_i(i, orientations)
            if orientations[j][k] != []:
                new.append(i)
        if len(new) == 0:
            return False
        elif len(new) == 1:
            return new[0]
        elif len(new) > 1:
            return choose(new)

def merge_coordinate_molecular(coor, lattice, wyckoff, sg, tol, orientations):
    while True:
        pairs, graph = find_short_dist(coor, lattice, tol)
        index = None
        valid = True
        if len(pairs)>0 and valid is True:
            if len(coor) > len(wyckoff[-1][0]):
                merged = []
                groups = connected_components(graph)
                for group in groups:
                    merged.append(get_center(coor[group], lattice))
                merged = np.array(merged)
                #if check_wyckoff_position(merged, sg, wyckoff) is not False:
                index = check_wyckoff_position_molecular(merged, sg, orientations, exact_translation=False)
                if index is False:
                    return coor, False
                elif index is None:
                    valid = False
                else:
                    #Check each possible merged Wyckoff position for orientaitons
                    coor = merged

            else:#no way to merge
                #print('no way to Merge, FFFFFFFFFFFFFFFFFFFFFFF----------------')
                return coor, False
        else:
            if index is None:
                index = check_wyckoff_position_molecular(coor, sg, orientations, exact_translation=False)
            return coor, index

def choose_wyckoff_molecular(wyckoffs, number, orientations):
    """
    choose the wyckoff sites based on the current number of atoms
    rules 
    1, the newly added sites is equal/less than the required number.
    2, prefer the sites with large multiplicity
    orientations: the valid orientations --for a given molecule--
    """
    if np.random.random()>0.5: #choose from high to low
        for j, wyckoff in enumerate(wyckoffs):
            if len(wyckoff[0]) <= number:
                good_wyckoff = []
                for k, w in enumerate(wyckoff):
                    if orientations[j][k] != []:
                        good_wyckoff.append([j,k])
                if len(good_wyckoff) > 0:
                    for indices in good_wyckoff:
                        if orientations[indices[0]][indices[1]] == []:
                            print(str(j)+", "+str(k)+str(" X"))
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
            for indices in good_wyckoff:
                j, k = indices
                if orientations[j][k] == []:
                    print(str(j)+", "+str(k)+str(" Y"))
            return choose(good_wyckoff)
        else:
            return False

class mol_site():
    '''
    Class for storing molecular Wyckoff positions and orientations within
    the molecular_crystal class. 
    '''
    def __init__(self, mol, position, sg, wp_index, lattice):
        #Pymatgen molecule object
        self.mol = mol
        #Relative coordinates within the unit cell
        self.position = position
        #Spacegroup number
        self.sg = sg
        #single index of the Wyckoff position within the spacegroup
        self.wp_index = wp_index
        #letter of the Wyckoff position
        self.multiplicity = len(get_wyckoffs(sg)[wp_index])
        self.letter = letter_from_index(wp_index, sg)

class molecular_crystal():
    '''
    Class for storing and generating molecular crystals based on symmetry
    constraints. Based on the random_crystal class for atomic crystals.
    Given a spacegroup, list of molecule objects, molecular stoichiometry, and
    a volume factor, generates a molecular crystal consistent with the given
    constraints. This crystal is stored as a pymatgen struct via self.struct
    '''
    def __init__(self, sg, molecules, numMols, factor, allow_inversion=False, orientations=None):
        
        #Necessary input
        numMols = np.array(numMols) #must convert it to np.array
        self.factor = factor
        self.numMols0 = numMols
        self.sg = sg
        #Reorient the molecules along their principle axes
        oriented_molecules = []
        for mol in molecules:
            pga = PointGroupAnalyzer(mol)
            mo = pga.symmetrize_molecule()['sym_mol']
            oriented_molecules.append(mo)
        self.molecules = oriented_molecules
        self.boxes = []
        #Calculate binding boxes for each molecule
        for mol in self.molecules:
            self.boxes.append(get_box(reoriented_molecule(mol)[0]))
        self.radii = []
        self.minlen = []
        self.maxlen = []
        for box in self.boxes:
            self.radii.append(math.sqrt( max(box[1],box[0])**2 + max(box[3],box[2])**2 + max(box[5],box[4])**2 ))
            lens = [box[1]-box[0], box[3]-box[2], box[5]-box[4]]
            self.minlen.append(min(lens))
            self.maxlen.append(min(lens))
        self.Msgs()
        self.numMols = numMols * cellsize(self.sg)
        self.volume = estimate_volume_molecular(self.numMols, self.boxes, self.factor)
        self.wyckoffs = get_wyckoffs(self.sg, organized=True) #2D Array of Wyckoff positions organized by multiplicity
        #Whether or not to allow chiral molecules to be flipped
        self.allow_inversion = allow_inversion
        #When generating multiple crystals of the same stoichiometry and sg,
        #allow the user to re-use the allowed orientations, to reduce time cost
        if orientations is None:
            self.get_orientations()
        else:
            self.valid_orientations = orientations
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
        Calculate the valid orientations for each Molecule and Wyckoff position.
        Returns a list with 4 indices:
        index 1: the molecular prototype's index
        index 2: the Wyckoff position's 1st index (based on multiplicity)
        index 3: the WP's 2nd index (within the group of equal multiplicity)
        index 4: the index of the valid orientation for the molecule/WP pair
        For example, self.valid_orientations[i][j][k] would be a list of valid
            orientations for self.molecules[i],
            in the Wyckoff position self.wyckoffs[j][k]
        """
        self.valid_orientations = []
        for mol in self.molecules:
            self.valid_orientations.append([])
            wp_index = -1
            for i, x in enumerate(self.wyckoffs):
                self.valid_orientations[-1].append([])
                for j, wp in enumerate(x):
                    wp_index += 1
                    allowed = orientation_in_wyckoff_position(mol, self.sg, wp_index, already_oriented=True, allow_inversion=self.allow_inversion)
                    if allowed is not False:
                        self.valid_orientations[-1][-1].append(allowed)
                    else:
                        self.valid_orientations[-1][-1].append([])

    def check_compatible(self):
        '''
        check if the number of molecules is compatible with the
        wyckoff positions
        needs to improve later
        '''
        N_site = [len(x[0]) for x in self.wyckoffs]
        has_freedom = False
        #remove WP's with no freedom once they are filled
        removed_wyckoffs = []
        for i, numMol in enumerate(self.numMols):
            #Check that the number of molecules is a multiple of the smallest Wyckoff position
            if numMol % N_site[-1] > 0:
                return False
            else:
                #Check if smallest WP has at least one degree of freedom
                op = self.wyckoffs[-1][-1][0]
                if op.rotation_matrix.all() != 0.0:
                    if self.valid_orientations[i][-1][-1] != []:
                        has_freedom = True
                else:
                    #Subtract from the number of ions beginning with the smallest Wyckoff positions
                    remaining = numMol
                    for j, x in enumerate(self.wyckoffs):
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
        """the main code to generate random crystal """
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
            minvector = max(lens for lens in self.minlen)
            #print(self.radii, minvector)
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

                    coordinates_total = [] #to store the added coordinates
                    sites_total = []      #to store the corresponding molecular specie
                    wps_total = []      #to store corresponding Wyckoff position indices
                    points_total = []   #to store the generating x,y,z points
                    good_structure = False

                    for cycle2 in range(max2):
                        coordinates_tmp = deepcopy(coordinates_total)
                        sites_tmp = deepcopy(sites_total)
                        wps_tmp = deepcopy(wps_total)
                        points_tmp = deepcopy(points_total)
                        
                	    #Add molecules specie by specie
                        for numMol, mol in zip(self.numMols, self.molecules):
                            i = self.molecules.index(mol)
                            numMol_added = 0

                            #Now we start to add the specie to the wyckoff position
                            for cycle3 in range(max3):
                                #Choose a random Wyckoff position for given multiplicity: 2a, 2b, 2c
                                #NOTE: The molecular version return wyckoff indices, not ops
                                indices = choose_wyckoff_molecular(self.wyckoffs, numMol-numMol_added, self.valid_orientations[i])
                                if indices is not False:
                                    j, k = indices
                                    if self.valid_orientations[i][j][k] == []:
                                        print("Error: Failed to catch empty set...")
                                        print(i,j,k)
                	    	    #Generate a list of coords from ops
                                    ops = self.wyckoffs[j][k]
                                    point = np.random.random(3)
                                    #print('generating new points:', point)
                                    coords = np.array([op.operate(point) for op in ops])
                                    #merge_coordinate if the atoms are close
                                    coords_toadd, good_merge = merge_coordinate_molecular(coords, cell_matrix, 
                                            self.wyckoffs, self.sg, 3.0, self.valid_orientations[i])
                                    if good_merge is not False:
                                        wp_index = good_merge
                                        point = coords_toadd[0]
                                        coords_toadd -= np.floor(coords_toadd) #scale the coordinates to [0,1], very important!
                                        if check_distance_molecular(coordinates_tmp, coords_toadd, sites_tmp, i, cell_matrix, self.radii):
                                            coordinates_tmp.append(coords_toadd)
                                            sites_tmp.append(i)
                                            wps_tmp.append(wp_index)
                                            points_tmp.append(point)
                                            numMol_added += len(coords_toadd)
                                        if numMol_added == numMol:
                                            coordinates_total = deepcopy(coordinates_tmp)
                                            sites_total = deepcopy(sites_tmp)
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
                            coordinates_total = []
                            sites_total = []
                            wps_total = []
                    #placing molecules here
                    if good_structure:
                        final_lattice = cell_matrix 
                        final_coor = []
                        final_site = []
                        final_number = []
                        wyckoffs = get_wyckoffs(sg)
                        self.mol_generators = []
                        for center0, i, wp_index in zip(points_total, sites_total, wps_total):
                            mol = self.molecules[i]
                            mo = deepcopy(mol)
                            #get j, k from wp_index
                            num = 0
                            found = False
                            j, k = jk_from_i(wp_index, self.wyckoffs)
                            op1 = choose(self.valid_orientations[i][j][k]).get_op()
                            mo.apply_operation(op1)
                            ms0 = mol_site(mo, center0, self.sg, wp_index, cell_matrix)
                            self.mol_generators.append(ms0)
                            for index, op2 in enumerate(get_wyckoff_generators(self.sg)[wp_index]):
                                for site in mol:
                                    #Place molecular coordinates in relative coordinates
                                    relative_coords = np.dot(np.linalg.inv(cell_matrix), site.coords)
                                    raw_vector = center0 + relative_coords
                                    new_vector = op2.operate(raw_vector)
                                    new_vector -= np.floor(new_vector)
                                    final_coor.append(new_vector)
                                    final_site.append(site.species_string)
                                    final_number.append(site.specie.number)

                        final_coor -= np.floor(final_coor)
                        if verify_distances(final_coor, final_site, final_lattice) is True:      
                            self.lattice = final_lattice  
                            self.coordinates = np.array(final_coor)
                            self.sites = final_site              
                            self.struct = Structure(final_lattice, final_site, np.array(final_coor))
                            self.spg_struct = (final_lattice, np.array(final_coor), final_number)
                            self.valid = True
                            return
                        #else: print("Failed final distance check.")
        print("Couldn't generate crystal after max attempts.")
        if degrees == 0:
            print("Note: Wyckoff positions have no degrees of freedom.")
        self.struct = self.Msg2
        self.valid = False
        return self.Msg2


if __name__ == "__main__":
    #-------------------------------- Options -------------------------
    from ase.build import molecule as ase_molecule
    from pymatgen import Molecule
    def get_ase_mol(molname):
        """convert ase molecule to pymatgen style"""
        ase_mol = ase_molecule(molname)
        pos = ase_mol.get_positions()
        symbols = ase_mol.get_chemical_symbols()
        return(Molecule(symbols, pos))
    
    #-------------------------------- Options -------------------------
    parser = OptionParser()
    parser.add_option("-s", "--spacegroup", dest="sg", metavar='sg', default=36, type=int,
            help="desired space group number: 1-230, e.g., 194")
    parser.add_option("-e", "--molecule", dest="molecule", default='H2O', 
            help="desired molecules: e.g., H2O", metavar="molecule")
    parser.add_option("-n", "--numMols", dest="numMols", default=2, 
            help="desired numbers of molecules: 4", metavar="numMols")
    parser.add_option("-f", "--factor", dest="factor", default=3.0, type=float, 
            help="volume factor: default 3.0", metavar="factor")
    parser.add_option("-v", "--verbosity", dest="verbosity", default=0, type=int, help="verbosity: default 0; higher values print more information", metavar="verbosity")
    parser.add_option("-a", "--attempts", dest="attempts", default=1, type=int, 
            help="number of crystals to generate: default 1", metavar="attempts")

    (options, args) = parser.parse_args()    
    molecule = options.molecule
    number = options.numMols
    verbosity = options.verbosity
    attempts = options.attempts
    
    numMols = []
    if molecule.find(',') > 0:
        strings = molecule.split(',')
        system = []
        for mol in strings:
            system.append(get_ase_mol(mol))
        for x in number.split(','):
            numMols.append(int(x))
    else:
        system = [get_ase_mol(molecule)]
        numMols = [int(number)]
    orientations = None
    for i in range(attempts):
        numMols0 = np.array(numMols)
        sg = options.sg
        rand_crystal = molecular_crystal(options.sg, system, numMols0, options.factor, orientations=orientations)

        if rand_crystal.valid:
            #orientations = rand_crystal.valid_orientations
            #pymatgen style
            #print("Generated number "+str(i+1))
            CifWriter(rand_crystal.struct, symprec=0.1).write_file(filename = "out/"+str(i+1)+'.cif')

            #spglib style structure called cell
            ans = get_symmetry_dataset(rand_crystal.spg_struct, symprec=1e-1)['number']
            print('Space group requested: ', sg, 'generated', ans, 'vol: ', rand_crystal.volume)
            if verbosity > 0:
                for ms in rand_crystal.mol_generators:
                    print(str(ms.multiplicity)+str(ms.letter)+" "+str(ms.position))

            #print(CifWriter(new_struct, symprec=0.1).__str__())
            #print('Space group:', finder.get_space_group_symbol(), 'tolerance:', tol)
            #output wyckoff sites only

        else: 
            print('something is wrong')
            #print(len(rand_crystal.coordinates))
            break
            #print(new_struct)
