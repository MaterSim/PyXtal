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

def get_box(mol, padding=1):
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
            oriented_molecules.append(reoriented_molecule(mol)[0])
        self.molecules = oriented_molecules
        self.boxes = []
        #Calculate binding boxes for each molecule
        for mol in self.molecules:
            self.boxes.append(get_box(mol))
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
                    elif allowed is False:
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
        for numMol in self.numMols:
            #Check that the number of molecules is a multiple of the smallest Wyckoff position
            if numMol % N_site[-1] > 0:
                return False
            else:
                #Check if smallest WP has at least one degree of freedom
                op = self.wyckoffs[-1][-1][0]
                if op.rotation_matrix.all() != 0.0:
                    has_freedom = True
                else:
                    #Subtract from the number of ions beginning with the smallest Wyckoff positions
                    remaining = numMol
                    for x in self.wyckoffs:
                        for wp in x:
                            removed = False
                            while remaining >= len(wp) and wp not in removed_wyckoffs:
                                #Check if WP has at least one degree of freedom
                                op = wp[0]
                                remaining -= len(wp)
                                if np.allclose(op.rotation_matrix, np.zeros([3,3])):
                                    removed_wyckoffs.append(wp)
                                    removed = True
                                else:
                                    has_freedom = True
                    if remaining != 0:
                        return False
        if has_freedom:
            return True
        else:
            #print("Warning: Wyckoff Positions have no degrees of freedom.")
            return 0

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
            if degrees is 0:
                max1 = 5
                max2 = 5
                max3 = 5
            #Calculate a minimum vector length for generating a lattice
            minvector = max(max(2.0*radius for radius in self.radii), tol_m)
            for cycle1 in range(max1):
                #1, Generate a lattice
                cell_para = generate_lattice(self.sg, self.volume, minvec=minvector)
                cell_matrix = para2matrix(cell_para)
                if abs(self.volume - np.linalg.det(cell_matrix)) > 1.0: 
                    print('Error, volume is not equal to the estimated value: ', self.volume, ' -> ', np.linalg.det(cell_matrix))
                    print('cell_para:  ', cell_para)
                    sys.exit(0)

                coordinates_total = [] #to store the added coordinates
                sites_total = []      #to store the corresponding specie
                good_structure = False

                for cycle2 in range(max2):
                    coordinates_tmp = deepcopy(coordinates_total)
                    sites_tmp = deepcopy(sites_total)
                    
            	    #Add specie by specie
                    for numIon, radius in zip(self.numIons, self.radii):
                        numIon_added = 0
                        tol = max(0.5*radius, tol_m)

                        #Now we start to add the specie to the wyckoff position
                        for cycle3 in range(max3):
                            #Choose a random Wyckoff position for given multiplicity: 2a, 2b, 2c
                            ops = choose_wyckoff(self.wyckoffs, numIon-numIon_added) 
                            if ops is not False:
            	    	    #Generate a list of coords from ops
                                point = np.random.random(3)
                                #print('generating new points:', point)
                                coords = np.array([op.operate(point) for op in ops])
                                #merge_coordinate if the atoms are close
                                coords_toadd, good_merge = merge_coordinate(coords, cell_matrix, self.wyckoffs, self.sg, tol)
                                if good_merge:
                                    coords_toadd -= np.floor(coords_toadd) #scale the coordinates to [0,1], very important!
                                    #print('existing: ', coordinates_tmp)
                                    if check_distance(coordinates_tmp, coords_toadd, sites_tmp, specie, cell_matrix):
                                        coordinates_tmp.append(coords_toadd)
                                        sites_tmp.append(specie)
                                        numIon_added += len(coords_toadd)
                                    if numIon_added == numIon:
                                        coordinates_total = deepcopy(coordinates_tmp)
                                        sites_total = deepcopy(sites_tmp)
                                        break
                        if numIon_added != numIon:
                            break  #need to repeat from the 1st species

                    if numIon_added == numIon:
                        #print(self.Msg6)
                        good_structure = True
                        break
                    else: #reset the coordinates and sites
                        coordinates_total = []
                        sites_total = []

                if good_structure:
                    final_coor = []
                    final_site = []
                    final_number = []
                    final_lattice = cell_matrix
                    for coor, ele in zip(coordinates_total, sites_total):
                        for x in coor:
                            final_coor.append(x)
                            final_site.append(ele)
                            final_number.append(Element(ele).z)

                    self.lattice = final_lattice                    
                    self.coordinates = np.array(final_coor)
                    self.sites = final_site                    
                    self.struct = Structure(final_lattice, final_site, np.array(final_coor))
                    self.spg_struct = (final_lattice, np.array(final_coor), final_number)
                    self.valid = True
                    return
        if degrees == 0: print("Wyckoff positions have no degrees of freedom.")
        self.struct = self.Msg2
        self.valid = False
        return self.Msg2
