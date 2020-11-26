# Standard Libraries
from copy import deepcopy
from random import choice
import numpy as np

# PyXtal imports #avoid *
from pyxtal.version import __version__
from pyxtal.molecular_crystal import (
    molecular_crystal,
    molecular_crystal_2D,
    molecular_crystal_1D,
)
from pyxtal.crystal import (
    random_cluster,
    random_crystal,
    random_crystal_1D,
    random_crystal_2D,
)
from pyxtal.symmetry import Group, Wyckoff_position
from pyxtal.wyckoff_site import atom_site, mol_site
from pyxtal.wyckoff_split import wyckoff_split
from pyxtal.lattice import Lattice
from pyxtal.operations import apply_ops
from pyxtal.tolerance import Tol_matrix
from pyxtal.io import write_cif, structure_from_ext
from pyxtal.XRD import XRD

name = "pyxtal"

def print_logo():
    """
    print out the logo and version information
    """

    print(
        """
             ______       _    _          _   
            (_____ \     \ \  / /        | |   
             _____) )   _ \ \/ / |_  ____| |  
            |  ____/ | | | )  (|  _)/ _  | | 
            | |    | |_| |/ /\ \ |_( (_| | |___
            |_|     \__  /_/  \_\___)__|_|_____)
                   (____/      """
    )
    print("\n")
    print("----------------------(version", __version__, ")----------------------\n")
    print("A Python package for random crystal generation")
    print("The source code is available at https://github.com/qzhu2017/pyxtal")
    print("Developed by Zhu's group at University of Nevada Las Vegas\n\n")

class pyxtal:
    """
    Class for handling atomic crystals based on symmetry constraints. 

    Examples
    --------
    >>> from pyxtal import pyxtal
    >>> struc = pyxtal()
    >>> struc.from_random(3, 227, ['C'], [8])
    >>> struc
    ------Crystal from random------
    Dimension: 3
    Composition: C8
    Group: Fd-3m (227)
    cubic lattice:   4.3529   4.3529   4.3529  90.0000  90.0000  90.0000
    Wyckoff sites:
	 C @ [0.1250 0.1250 0.1250], WP:  8a, Site symmetry: -4 3 m

    The structure object can be easily manipulated via `apply_perturbtation`
    or `subgroup` function

    >>> struc2 = struc.subgroup(H=141, once=True)
    >>> struc2
    ------Crystal from Wyckoff Split------
    Dimension: 3
    Composition: C8
    Group: I41/amd (141)
    tetragonal lattice:   3.3535   3.3535   4.6461  90.0000  90.0000  90.0000
    Wyckoff sites:
    	 C @ [0.0000 0.2500 0.3750], WP:  4b, Site symmetry: -4 m 2

    Alternatively, one can also easily compute its XRD via the `pyxtal.XRD` class

    >>> xrd = struc.get_XRD()
    >>> xrd
      2theta     d_hkl     hkl       Intensity  Multi
      32.706     2.738   [ 1  1  1]   100.00        8
      54.745     1.677   [ 2  2  0]    40.95       12
      65.249     1.430   [ 3  1  1]    20.65       24
      81.116     1.186   [ 4  0  0]     5.15        6
      90.236     1.088   [ 3  3  1]     8.24       24
     105.566     0.968   [ 4  2  2]    14.44       24
     115.271     0.913   [ 5  1  1]    10.03       24
     133.720     0.838   [ 4  4  0]     9.80       12
     148.177     0.802   [ 5  3  1]    28.27       48

    Finally, the structure can be saved to different formats

    >>> struc.to_file('my.cif')
    >>> struc.to_file('my_poscar', fmt='poscar')

    or to Pymatgen/ASE structure object

    >>> pmg_struc = struc.to_pymatgen()
    >>> ase_struc = struc.to_ase()
    """

    def __init__(self, molecular=False):
        self.valid = False
        self.molecular = molecular
        self.diag = False
        self.numIons = None
        self.numMols = None

    def __str__(self):
        if self.valid:
            s = "------Crystal from {:s}------".format(self.source)
            s += "\nDimension: {}".format(self.dim)
            s += "\nComposition: {}".format(self.formula)
            if self.group.number in [7, 14, 15] and self.diag:
                symbol = self.group.alias
            else:
                symbol = self.group.symbol
            s += "\nGroup: {} ({})".format(symbol, self.group.number)
            s += "\n{}".format(self.lattice)
            s += "\nWyckoff sites:"
            if self.molecular:
                for wyc in self.mol_sites:
                    s += "\n\t{}".format(wyc)
            else:
                for wyc in self.atom_sites:
                    s += "\n\t{}".format(wyc)
        else:
            s = "\nStructure not available."
        return s

    def __repr__(self):
        return str(self)

    def from_random(
        self,
        dim = 3,
        group=None,
        species=None,
        numIons=None,
        factor=1.1,
        thickness = None,
        area = None,
        lattice=None,
        sites = None,
        conventional = True,
        diag = False,
        t_factor = 1.0,
        max_count = 10,
        force_pass = False,
    ):
        if self.molecular:
            prototype = "molecular"
        else:
            prototype = "atomic"
        tm = Tol_matrix(prototype=prototype, factor=t_factor)

        count = 0
        quit = False

        while True:
            count += 1
            if self.molecular:
                if dim == 3:
                    struc = molecular_crystal(group, species, numIons, factor, 
                    lattice=lattice, sites=sites, conventional=conventional, diag=diag, tm=tm)
                elif dim == 2:
                    struc = molecular_crystal_2D(group, species, numIons, factor, 
                    thickness=thickness, sites=sites, conventional=conventional, tm=tm)
                elif dim == 1:
                    struc = molecular_crystal_1D(group, species, numIons, factor, 
                    area=area, sites=sites, conventional=conventional, tm=tm)
            else:
                if dim == 3:
                    struc = random_crystal(group, species, numIons, factor, 
                            lattice, sites, conventional, tm)
                elif dim == 2:
                    struc = random_crystal_2D(group, species, numIons, factor, 
                            thickness, lattice, sites, conventional, tm)
                elif dim == 1:
                    struc = random_crystal_1D(group, species, numIons, factor, 
                            area, lattice, sites, conventional, tm)
                else:
                    struc = random_cluster(group, species, numIons, factor, 
                            lattice, sites, tm)
            if force_pass:
                quit = True
                break
            elif struc.valid:
                quit = True
                break

            if count >= max_count:
                raise RuntimeError("It takes long time to generate the structure, check inputs")

        if quit:
            self.valid = struc.valid
            self.dim = dim
            try:
                self.lattice = struc.lattice
                if self.molecular:
                    self.numMols = struc.numMols
                    self.molecules = struc.molecules
                    self.mol_sites = struc.mol_sites
                    self.diag = struc.diag
                else:
                    self.numIons = struc.numIons
                    self.species = struc.species
                    self.atom_sites = struc.atom_sites
                self.group = struc.group
                self.PBC = struc.PBC
                self.source = 'random'
                self.factor = struc.factor
                self.number = struc.number
                self.get_formula()
            except:
                pass


    def from_seed(self, seed, molecule=None, relax_h=False):
        """
        Load the seed structure from Pymatgen/ASE/POSCAR/CIFs
        Internally they will be handled by Pymatgen
        """
        from ase import Atoms
        from pymatgen import Structure

        if self.molecular:
            from pyxtal.molecule import pyxtal_molecule
            pmol = pyxtal_molecule(molecule).mol
            struc = structure_from_ext(seed, pmol, relax_h=relax_h)
            if struc.match():
                self.mol_sites = [struc.make_mol_site()]
                self.group = Group(struc.wyc.number)
                self.lattice = struc.lattice
                self.molecules = [pyxtal_molecule(struc.molecule, symmetrize=False)]
                self.numMols = struc.numMols
                self.diag = struc.diag
                self.valid = True # Need to add a check function
            else:
                raise ValueError("Cannot extract the molecular crystal from cif")
        else:
            if isinstance(seed, dict):
                self.from_dict()
            elif isinstance(seed, Atoms): #ASE atoms
                from pymatgen.io.ase import AseAtomsAdaptor
                pmg_struc = AseAtomsAdaptor.get_structure(seed)
                self.from_pymatgen(pmg_struc)
            elif isinstance(seed, Structure): #Pymatgen
                self.from_pymatgen(seed)
            elif isinstance(seed, str):
                pmg_struc = Structure.from_file(seed)
                self.from_pymatgen(pmg_struc)

        self.factor = 1.0
        self.number = self.group.number
        self.source = 'Seed'
        self.dim = 3
        self.PBC = [1, 1, 1]
        self.get_formula()

    def get_formula(self):
        formula = ""
        if self.molecular:
            numspecies = self.numMols
            species = [str(mol) for mol in self.molecules]
        else:
            numspecies = self.numIons
            species = self.species
        for i, s in zip(numspecies, species):
            formula += "{:s}{:d}".format(s, int(i))
        self.formula = formula

    def from_pymatgen(self, structure):
        """
        Load the seed structure from Pymatgen/ASE/POSCAR/CIFs
        """
        from pymatgen.symmetry.analyzer import SpacegroupAnalyzer as sga

        self.valid = True
        try:
            # needs to do it twice in order to get the conventional cell
            s = sga(structure)
            structure = s.get_refined_structure()
            s = sga(structure)
            sym_struc = s.get_symmetrized_structure()
            number = s.get_space_group_number()
        except:
            print("Failed to load the Pymatgen structure")
            self.valid = False

        if self.valid:
            d = sym_struc.composition.as_dict()
            species = [key for key in d.keys()]
            numIons = []
            for ele in species:
                numIons.append(int(d[ele]))
            self.numIons = numIons
            self.species = species
            self.group = Group(number)
            atom_sites = []
            for i, site in enumerate(sym_struc.equivalent_sites):
                pos = site[0].frac_coords
                wp = Wyckoff_position.from_group_and_index(number, sym_struc.wyckoff_symbols[i])
                specie = site[0].specie.number
                atom_sites.append(atom_site(wp, pos, specie))
            self.atom_sites = atom_sites
            matrix, ltype = sym_struc.lattice.matrix, self.group.lattice_type
            self.lattice = Lattice.from_matrix(matrix, ltype=ltype)

    def check_short_distances(self, r=0.7, exclude_H = True):
        """
        A function to check short distance pairs
        Mainly used for debug, powered by pymatgen

        Args:
            r: the given cutoff distances
            exclude_H: whether or not exclude the H atoms

        Returns:
            list of pairs within the cutoff
        """
        if self.dim > 0:
            pairs = []
            pmg_struc = self.to_pymatgen()
            if exclude_H:
                pmg_struc.remove_species('H')
            res = pmg_struc.get_all_neighbors(r)
            for i, neighs in enumerate(res):
                for n in neighs:
                    pairs.append([pmg_struc.sites[i].specie, n.specie, n.nn_distance])
        else:
            raise NotImplementedError("Does not support cluster for now")
        return pairs

    def to_file(self, filename=None, fmt=None, permission='w', sym_num=None):
        """
        Creates a file with the given filename and file type to store the structure.
        By default, creates cif files for crystals and xyz files for clusters.
        For other formats, Pymatgen is used

        Args:
            filename: the file path
            fmt: the file type (`cif`, `xyz`, etc.)
            permission: `w` or `a+`

        Returns:
            Nothing. Creates a file at the specified path
        """
        if self.valid:
            if fmt is None:
                if self.dim == 0:
                    fmt = 'xyz'
                else:
                    fmt = 'cif'

            if fmt == "cif":
                if self.dim == 3:
                    return write_cif(self, filename, "from_pyxtal", permission, sym_num=sym_num)
                else:
                    pmg_struc = self.to_pymatgen()
                    if self.molecular:
                        pmg_struc.sort()
                return pmg_struc.to(fmt=fmt, filename=filename)
            else:
                pmg_struc = self.to_pymatgen()
                if self.molecular:
                    pmg_struc.sort()
                return pmg_struc.to(fmt=fmt, filename=filename)
        else:
            raise RuntimeError("Cannot create file: structure did not generate")


    def subgroup(self, H=None, eps=0.05, idx=None, once=False, group_type='t'):
        """
        generate a structure with lower symmetry

        Args:
            H: space group number (int)
            eps: pertubation term (float)
            idx: list
            once: generate only one structure, otherwise output all

        Returns:
            a list of pyxtal structures with lower symmetries
        """

        #randomly choose a subgroup from the available list
        if group_type == 't':
            Hs = self.group.get_max_t_subgroup()['subgroup']
        else:
            Hs = self.group.get_max_k_subgroup()['subgroup']

        if idx is None:
            idx = range(len(Hs))
        else:
            for id in idx:
                if id >= len(Hs):
                    raise ValueError("The idx exceeds the number of possible splits")
        if H is not None: 
            idx = [id for id in idx if Hs[id] == H]

        if len(idx) == 0:
            raise ValueError("The space group H is incompatible with idx")

        sites = [str(site.wp.multiplicity)+site.wp.letter for site in self.atom_sites]
        valid_splitters = []
        bad_splitters = []
        for id in idx:
            splitter = wyckoff_split(G=self.group.number, wp1=sites, idx=id, group_type=group_type)
            if splitter.valid_split:
                valid_splitters.append(splitter)
            else:
                bad_splitters.append(splitter)

        if len(valid_splitters) == 0:
            # do one more step
            new_strucs = []
            for splitter in bad_splitters:
                trail_struc = self.subgroup_by_splitter(splitter)
                new_strucs.append(trail_struc.subgroup(once=True, group_type=group_type))
            return new_strucs
        else:
            if once:
                return self.subgroup_by_splitter(choice(valid_splitters), eps=eps)
            else:
                new_strucs = []
                for splitter in valid_splitters:
                    new_strucs.append(self.subgroup_by_splitter(splitter, eps=eps))
            return new_strucs

    def subgroup_by_splitter(self, splitter, eps=0.05):
        lat1 = np.dot(splitter.R[:3,:3].T, self.lattice.matrix)
        multiples = np.linalg.det(splitter.R[:3,:3])
        split_sites = []
        for i, site in enumerate(self.atom_sites):
            pos = site.position
            for ops1, ops2 in zip(splitter.G2_orbits[i], splitter.H_orbits[i]):
                pos0 = apply_ops(pos, ops1)[0]
                pos0 -= np.floor(pos0)
                pos0 += eps*(np.random.sample(3) - 0.5)
                wp, _ = Wyckoff_position.from_symops(ops2, group=splitter.H.number, permutation=False)
                split_sites.append(atom_site(wp, pos0, site.specie))
        new_struc = deepcopy(self)
        new_struc.group = splitter.H
        lattice = Lattice.from_matrix(lat1, ltype=new_struc.group.lattice_type)
        new_struc.lattice=lattice.mutate(degree=eps, frozen=True)
        new_struc.atom_sites = split_sites
        new_struc.numIons = [int(multiples*numIon) for numIon in self.numIons]
        new_struc.source = 'Wyckoff Split'
        
        return new_struc

    def apply_perturbation(self, d_lat=0.05, d_coor=0.05, d_rot=1):
        """
        perturb the structure without breaking the symmetry

        Args:
            d_coor: magnitude of perturbation on atomic coordinates (in A)
            d_lat: magnitude of perturbation on lattice (in percentage)
        """

        self.lattice = self.lattice.mutate(degree=d_lat)

        if self.molecular:
            for i, site in enumerate(self.mol_sites):
                site.perturbate(lattice=self.lattice.matrix, trans=d_coor, rot=d_rot)
        else:
            for i, site in enumerate(self.atom_sites):
                site.perturbate(lattice=self.lattice.matrix, magnitude=d_coor)

        self.source = 'Perturbation'

    def copy(self):
        """
        simply copy the structure
        """
        return deepcopy(self)

    def _get_coords_and_species(self, absolute=False, unitcell=True):
        """
        extract the coordinates and species information 

        Args:
            abosulte: if True, return the cartesian coords otherwise fractional

        Returns:
            total_coords (N*3 numpy array) and the list of species
        """
        species = []
        total_coords = None
        if self.molecular:
            for site in self.mol_sites:
                coords, site_species = site.get_coords_and_species(absolute, unitcell=unitcell)
                species.extend(site_species)
                if total_coords is None:
                    total_coords = coords
                else:
                    total_coords = np.append(total_coords, coords, axis=0)
        else:
            for site in self.atom_sites:
                species.extend([site.specie]*site.multiplicity)
                if total_coords is None:
                    total_coords = site.coords
                else:
                    total_coords = np.append(total_coords, site.coords, axis=0)

            if absolute:
                total_coords = total_coords.dot(self.lattice.matrix)
        
        return total_coords, species

    def to_ase(self, resort=True):
        """
        export to ase Atoms object
        """
        from ase import Atoms
        if self.valid:
            if self.dim > 0:
                lattice = self.lattice.copy()
                if self.molecular:
                    coords, species = self._get_coords_and_species(True)
                    latt, coords = lattice.add_vacuum(coords, frac=False, PBC=self.PBC)
                    atoms = Atoms(species, positions=coords, cell=latt, pbc=self.PBC)
                    if resort:
                        permutation = np.argsort(atoms.numbers)
                        atoms = atoms[permutation]
                    return atoms
                else:
                    coords, species = self._get_coords_and_species()
                    latt, coords = lattice.add_vacuum(coords, PBC=self.PBC)
                    return Atoms(species, scaled_positions=coords, cell=latt, pbc=self.PBC)
            else:
                coords, species = self._get_coords_and_species(True)
                return Atoms(species, positions=coords)
        else:
            raise RuntimeError("No valid structure can be converted to ase.")

    def to_pymatgen(self):
        """
        export to Pymatgen structure object
        """
        from pymatgen.core.structure import Structure, Molecule

        if self.valid:
            if self.dim > 0:
                lattice = self.lattice.copy()
                coords, species = self._get_coords_and_species()
                # Add space above and below a 2D or 1D crystals
                latt, coords = lattice.add_vacuum(coords, PBC=self.PBC)
                return Structure(latt, species, coords)
            else:
                # Clusters are handled as large molecules
                coords, species = self._get_coords_and_species(True)
                return Molecule(species, coords)
        else:
            raise RuntimeError("No valid structure can be converted to pymatgen.")

    def get_XRD(self, **kwargs):
        """
        compute the PXRD object
        ** kwargs include
            - wavelength (1.54184)
            - thetas [0, 180]
            - preferred_orientation: False
            - march_parameter: None
        """

        return XRD(self.to_ase(), **kwargs)

    def show(self, **kwargs):
        """
        display the crystal structure
        """
        if self.molecular:
            from pyxtal.viz import display_molecular
            return display_molecular(self, **kwargs)
        else:
            from pyxtal.viz import display_atomic
            return display_atomic(self, **kwargs)

    def optimize_lattice(self):
        """
        optimize the lattice if the cell has a bad inclination angles
        """
        if self.molecular:
            for i in range(5):
                lattice, trans, opt = self.lattice.optimize()
                if opt:
                    for site in self.mol_sites:
                        pos_abs = np.dot(site.position, self.lattice.matrix)
                        pos_frac = pos_abs.dot(lattice.inv_matrix)
                        site.position = pos_frac - np.floor(pos_frac)
                        site.lattice = lattice
                        # for P21/c, Pc, C2/c, check if opt the inclination angle
                        if self.group.number in [7, 14, 15]:
                            for j, op in enumerate(site.wp.ops):
                                vec = op.translation_vector.dot(trans)
                                vec -= np.floor(vec)
                                op1 = op.from_rotation_and_translation(op.rotation_matrix, vec)
                                site.wp.ops[j] = op1
                    _, perm = Wyckoff_position.from_symops(site.wp.ops, self.group.number)            
                    if not isinstance(perm, list): 
                        self.diag = True
                    else:
                        self.diag = False
                    self.lattice = lattice
                else:
                    break


#    def save(self, filename=None):
#        """
#        Save the structure
#        Args:
#            filename: the file to save txt information
#        """
#        dict0 = self.save_dict(db_filename)
#        with open(filename, "w") as fp:
#            json.dump(dict0, fp)
#
#        print("save the GP model to", filename, "and database to", db_filename)
#
#    def load(self, filename=None):
#        """
#        load the structure
#        Args:
#            filename: the file to save txt information
#        """
#        #print(filename)
#        with open(filename, "r") as fp:
#            dict0 = json.load(fp)
#        self.load_from_dict(dict0, N_max=N_max, device=device)
#        self.fit(opt=opt)
#        print("load the GP model from ", filename)
#
    def save_dict(self):
        """
        save the model as a dictionary
        """
        sites = []
        if self.molecular:
            for site in self.mol_sites:
                sites.append(site.save_dict())
        else:
            for site in self.atom_sites:
                sites.append(site.save_dict())

        dict0 = {"lattice": self.lattice.matrix,
                 "sites": sites,
                 "group": self.group.number,
                 "molecular": self.molecular,
                 "numIons": self.numIons,
                 "numMols": self.numMols,
                 "factor": self.factor,
                 "PBC": self.PBC,
                 "formula": self.formula,
                 "source": self.source,
                 "dim": self.dim,
                 "valid": self.valid,
                }
       
        return dict0

    def load_dict(self, dict0):
        """
        load the structure from a dictionary
        """
        self.group = Group(dict0["group"])
        self.lattice = Lattice.from_matrix(dict0["lattice"], ltype=self.group.lattice_type)
        self.molecular = dict0["molecular"]
        self.number = self.group.number
        self.factor = dict0["factor"]
        self.source = dict0["source"]
        self.dim = dict0["dim"]
        self.PBC = dict0["PBC"]
        self.numIons = dict0["numIons"]
        self.numMols = dict0["numMols"]
        self.valid = dict0["valid"]
        self.formula = dict0["formula"]
        sites = []
        if dict0["molecular"]:
            for site in dict0["sites"]:
                sites.append(mol_site.load_dict(site))
            self.mol_sites = sites
        else:
            for site in dict0["sites"]:
                sites.append(atom_site.load_dict(site))
            self.atom_sites = sites
