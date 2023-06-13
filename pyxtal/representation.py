import os
import numpy as np
from pyxtal.symmetry import Group
from pyxtal.lattice import Lattice
from pyxtal.wyckoff_site import mol_site, atom_site
from pyxtal.molecule import find_rotor_from_smile
 
class representation_atom():
    """
    A class to handle the 1D representation of atomic crystal   
    Works for Zprime > 1

    Args:
        x: a list of [cell, site_1, site_2, ...]
    """

    def __init__(self, x):
        self.x = x

    def __str__(self):
        return self.to_string()

    @classmethod
    def from_pyxtal(cls, struc, standard=False):
        """
        Initialize 1D rep. from the pyxtal object

        Args:
            struc: pyxtal object
        """
        if standard and not struc.standard_setting:
            pmg = struc.to_pymatgen()
            struc.from_seed(pmg, standard=True)
        symmetry = [struc.atom_sites[0].wp.hall_number]
        lat = struc.lattice.encode()
        vector = [symmetry + lat]
        smiles = []
        for site in struc.atom_sites:
            vector.append(site.encode())
        x = vector
        return cls(x)
    
    def to_standard_setting(self):
        xtal = self.to_pyxtal()
        self.x = representation.from_pyxtal(xtal, standard=True).x
 
    def to_pyxtal(self):
        """
        Export the pyxtal structure

        Args:
            smiles: list of smiles
            compoisition: list of composition
        """
        from pyxtal import pyxtal

        # symmetry
        v = self.x[0]
        struc = pyxtal()
        struc.group, number = Group(v[0], use_hall=True), v[0]
    
        # lattice
        ltype = struc.group.lattice_type
        if ltype == 'triclinic':
            a, b, c, alpha, beta, gamma = v[1], v[2], v[3], v[4], v[5], v[6]
        elif ltype == 'monoclinic':
            a, b, c, alpha, beta, gamma = v[1], v[2], v[3], 90, v[4], 90
        elif ltype == 'orthorhombic':
            a, b, c, alpha, beta, gamma = v[1], v[2], v[3], 90, 90, 90
        elif ltype == 'tetragonal':
            a, b, c, alpha, beta, gamma = v[1], v[1], v[2], 90, 90, 90
        elif ltype == 'hexagonal':
            a, b, c, alpha, beta, gamma = v[1], v[1], v[2], 90, 90, 120
        else:
            a, b, c, alpha, beta, gamma = v[1], v[1], v[1], 90, 90, 90
        try:
            struc.lattice = Lattice.from_para(a, b, c, alpha, beta, gamma, ltype=ltype)
        except:
            print(a, b, c, alpha, beta, gamma, ltype)
            raise ValueError("Problem in Lattice")
    
        # sites
        struc.numIons = [0] * len(smiles) 
        struc.atom_sites = [] 

        count = 1
        for i, comp in enumerate(composition): 
            for j in range(comp):
                v = self.x[count]
                dicts = {}
                dicts['type'] = i
                dicts['dim'] = 3
                dicts['PBC'] = [1, 1, 1]
                dicts['hn'] = struc.group.hall_number
                dicts['index'] = 0
                dicts['lattice'] = struc.lattice.matrix
                dicts['lattice_type'] = ltype
                site = atom_site.from_1D_dicts(dicts)
                site.type = i
                struc.atom_sites.append(site)
                struc.numIons[i] += site.wp.multiplicity
                #move to next rep
                count += 1
            struc.species.append(site.specie)

        struc._get_formula()
        struc.source = '1D rep.'
        struc.valid = True
        struc.standard_setting = site.wp.is_standard_setting()

        return struc
    
    def to_array(self):
        """
        Export only varibles to a 1d numpy array
        """
        cells, xyzs = self.x[0][1:], self.x[1:]
        x = cells
        for xyz in xyzs: x = np.hstack((x, xyz[2:]))
        return x
        
    def to_string(self, time=None, eng=None, tag=None):
        """
        Export string representation

        Args:
            time: float
            eng: float
            tag: string
        """
        x = self.x
        strs = "{:3d} ".format(int(x[0][0]))

        # data for cell
        if x[0][0] <= 348:
            num = 4
        elif x[0][0] <= 488:
            num = 3
        else: #cubic
            num = 2

        for c in x[0][1:num]:
            strs += "{:5.2f} ".format(c)
        for c in x[0][num:]:
            strs += "{:5.1f} ".format(c)
        
        # data for atoms
        strs += "{:d} ".format(len(x)-1)  # Number of sites
        for i in range(1, len(x)):
            strs += "{:s} ".format(x[i][0])
            strs += "{:d} ".format(x[i][1])
            for v in x[i][2:]:
                strs += "{:4.2f} ".format(v)      

        if time is not None:
            strs += "{:5.2f}".format(time)

        if eng is not None:
            strs += "{:11.3f}".format(eng)
    
        if tag is not None:
            strs += " {:s}".format(tag)
    
        return strs


class representation():
    """
    A class to handle the 1D representation of molecular crystal   
    Works for Zprime > 1

    Args:
        x: a list of [cell, site_1, site_2, ...]
        smiles: a list of [smile_1, smile_2, ...]
    """

    def __init__(self, x, smiles=None):
        if smiles is not None:
            self.smiles = []
            for i, smile in enumerate(smiles):
                if smile.endswith('.smi'): 
                    smile = smile[:-4]
                self.smiles.append(smile)
        else:
            self.smiles = None
        self.x = x

    def __str__(self):
        return self.to_string()

    @classmethod
    def from_pyxtal(cls, struc, standard=False):
        """
        Initialize 1D rep. from the pyxtal object

        Args:
            struc: pyxtal object
        """
        if standard and not struc.standard_setting:
            #struc.optimize_lattice(standard=True)
            pmg = struc.to_pymatgen()
            struc.from_seed(pmg, molecules=struc.molecules, standard=True)
        symmetry = [struc.mol_sites[0].wp.hall_number]
        lat = struc.lattice.encode()
        vector = [symmetry + lat]
        smiles = []
        for site in struc.mol_sites:
            vector.append(site.encode())
            smiles.append(site.molecule.smile)
        x = vector
        return cls(x, smiles)
    
    @classmethod
    def from_string(cls, inputs, smiles, composition=None):
        """
        Initialize 1D rep. from the string

        Args:
            inputs: input string 
            smiles: list of smiles
        """
        #parse the cell
        if composition is None:
            composition = [1] * len(smiles)

        inputs = [float(tmp) for tmp in inputs.split()]
        hn = int(inputs[0])
        if hn <= 2:
            n_cell = 8
        elif hn <= 107:
            n_cell = 6
        elif hn <= 348:
            n_cell = 5
        elif hn <= 488:
            n_cell = 4
        else:
            n_cell = 3 #cubic
        cell = [hn] + inputs[1:n_cell-1]
        
        x = [cell]
        n_site = int(inputs[n_cell-1])
        if n_site != sum(composition):
            msg = "Composition is inconsistent: {:d}/{:d}\n".format(sum(composition), n_site)
            msg += str(inputs)
            raise ValueError(msg)
        n_cell += 1

        for i, smile in enumerate(smiles):
            if smile.endswith('.smi'): 
                smile=smile[:-4]
            for c in range(composition[i]):
                if smile in ["Cl-"]:
                    n_mol = 4
                else:
                    n_torsion = len(find_rotor_from_smile(smile))
                    n_mol = 7 + n_torsion
                #inversion
                inputs[n_cell+n_mol-2] = int(inputs[n_cell+n_mol-2])
                x.append(inputs[n_cell-1:n_cell+n_mol-1])
                n_cell += n_mol
        return cls(x, smiles)

    def to_standard_setting(self):
        xtal = self.to_pyxtal()
        self.x = representation.from_pyxtal(xtal, standard=True).x
 
    def to_pyxtal(self, smiles=None, composition=None):
        """
        Export the pyxtal structure

        Args:
            smiles: list of smiles
            compoisition: list of composition
        """
        from pyxtal import pyxtal
        if smiles is None:
            smiles = self.smiles

        if composition is None:
            composition = [1] * len(smiles)

        if sum(composition) + 1 != len(self.x):
            msg = "Composition is inconsistent:\n"
            msg += str(composition) + "\n"
            msg += self.to_string()
            raise ValueError(msg)

        # symmetry
        v = self.x[0]
        struc = pyxtal(molecular=True)
        struc.group, number = Group(v[0], use_hall=True), v[0]
    
        # lattice
        ltype = struc.group.lattice_type
        if ltype == 'triclinic':
            a, b, c, alpha, beta, gamma = v[1], v[2], v[3], v[4], v[5], v[6]
        elif ltype == 'monoclinic':
            a, b, c, alpha, beta, gamma = v[1], v[2], v[3], 90, v[4], 90
        elif ltype == 'orthorhombic':
            a, b, c, alpha, beta, gamma = v[1], v[2], v[3], 90, 90, 90
        elif ltype == 'tetragonal':
            a, b, c, alpha, beta, gamma = v[1], v[1], v[2], 90, 90, 90
        elif ltype == 'hexagonal':
            a, b, c, alpha, beta, gamma = v[1], v[1], v[2], 90, 90, 120
        else:
            a, b, c, alpha, beta, gamma = v[1], v[1], v[1], 90, 90, 90
        try:
            struc.lattice = Lattice.from_para(a, b, c, alpha, beta, gamma, ltype=ltype)
        except:
            print(a, b, c, alpha, beta, gamma, ltype)
            raise ValueError("Problem in Lattice")
    
        # sites
        struc.numMols = [0] * len(smiles) 
        struc.molecules = []
        struc.mol_sites = [] 

        count = 1
        for i, comp in enumerate(composition): 
            smile = smiles[i]
            if smile.endswith('.smi'): smile=smile[:-4]
            for j in range(comp):
                v = self.x[count]
                dicts = {}
                dicts['smile'] = smile
                dicts['type'] = i
                dicts['dim'] = 3
                dicts['PBC'] = [1, 1, 1]
                #dicts['number'] = number
                dicts['hn'] = struc.group.hall_number
                dicts['index'] = 0
                dicts['lattice'] = struc.lattice.matrix
                dicts['lattice_type'] = ltype
                dicts['center'] = v[:3]
                if smile not in ["Cl-"]:
                    dicts['orientation'] = np.array(v[3:6])
                    dicts['rotor'] = v[6:-1]
                    dicts['reflect'] = int(v[-1])
                site = mol_site.from_1D_dicts(dicts)
                site.type = i
                struc.mol_sites.append(site)
                struc.numMols[i] += site.wp.multiplicity
                #move to next rep
                count += 1
            struc.molecules.append(site.molecule)

        struc._get_formula()
        struc.source = '1D rep.'
        struc.valid = True
        struc.standard_setting = site.wp.is_standard_setting()

        return struc
    
    def to_string(self, time=None, eng=None, tag=None):
        """
        Export string representation

        Args:
            time: float
            eng: float
            tag: string
        """
        x = self.x
        strs = "{:3d} ".format(int(x[0][0]))

        # data for cell
        if x[0][0] <= 348:
            num = 4
        elif x[0][0] <= 488:
            num = 3
        else: #cubic
            num = 2

        for c in x[0][1:num]:
            strs += "{:5.2f} ".format(c)
        for c in x[0][num:]:
            strs += "{:5.1f} ".format(c)
        
        # data for molecule
        strs += "{:d} ".format(len(x)-1)
        for i in range(1, len(x)):
            for v in x[i][:3]:
                strs += "{:4.2f} ".format(v)      
            for v in x[i][3:-1]:
                strs += "{:6.1f} ".format(v)      
            strs += "{:d} ".format(int(x[i][-1]))

        if time is not None:
            strs += "{:5.2f}".format(time)

        if eng is not None:
            strs += "{:11.3f}".format(eng)
    
        if tag is not None:
            strs += " {:s}".format(tag)
    
        return strs

    def same_smiles(self, smiles):
        if len(self.smiles) == smiles:
            for s1, s2 in zip(self.smiles, smiles):
                if s2 != s2:
                    return False
            return True
        else:
            return False
            
    def get_dist(self, rep):
        """
        get distance with the other rep1
        Now only supports Z'=1
        """
        from pyxtal.symmetry import Wyckoff_position as WP
        if self.same_smiles(rep.smiles):
            msg = 'different smiles'
            print(msg)
            return None
        elif len(self.x) != len(rep.x):
            msg = 'different number of sites'
            print(msg)
            return None
        elif self.x[0][0] != rep.x[0][0]:
            msg = 'different space group numbers'
            print(msg)
            return None
        else:
            diffs = []
            wp = WP.from_group_and_index(self.x[0][0], 0, use_hall=True)
            for i in range(len(self.x)):
                diff = np.zeros(len(self.x[i]))
                tmp1 = np.array(self.x[i])
                tmp2 = np.array(rep.x[i])
                # cell difference
                if i == 0:
                    diff_cell = tmp2 - tmp1
                    diffs.extend(diff_cell)
                # site difference
                else:
                    # symmmetry variation
                    xyzs = wp.apply_ops(tmp2[:3])
                    diff_xyzs = xyzs - tmp1[:3]
                    diff_xyzs -= np.round(diff_xyzs)
                    id = np.argmin(np.linalg.norm(diff_xyzs, axis=1))
                    diff_xyz = diff_xyzs[id]
                    diff_ori = tmp2[3:6] - tmp1[3:6]
                    diff_ori /= [360.0, 180.0, 360.0]
                    diff_ori -= np.round(diff_ori)
                    diff_ori *= [360.0, 180.0, 360.0]
                    diff_tor = tmp2[6:] - tmp1[6:]
                    diff_tor /= 360.0
                    diff_tor -= np.round(diff_tor)
                    diff_tor *= 360.0
                    diffs.extend(diff_xyz)
                    diffs.extend(diff_ori)
                    diffs.extend(diff_tor)
            return np.array(diffs)

if __name__ == "__main__":

    #aspirin
    smiles = ['CC(=O)OC1=CC=CC=C1C(=O)O']
    x = [[81,11.43,6.49,11.19,83.31],[0.77,0.57,0.53,48.55,24.31,145.94,-77.85,-4.40,170.86,False]]
    #rep0 = representation(x, smiles)
    #print(rep0.to_string())
    rep1 = representation(x, smiles)
    xtal = rep1.to_pyxtal()
    print(xtal)
    rep2 = representation.from_pyxtal(xtal)
    print(rep2.to_pyxtal())
    print(rep2.to_string())
    string = "82 11.43  6.49 11.19 83.31 1  0.77  0.57  0.53 48.55 24.31 145.9 -77.85 -4.40 170.9 0"
    rep3 = representation.from_string(string, smiles)
    print(rep3.to_string())
    print(rep3.to_pyxtal())
    rep3.to_standard_setting()
    print(rep3.to_pyxtal())
    print(rep3.to_string())

    string1 = "81 14.08  6.36 25.31  83.9 1 0.83 0.40 0.63  136.6  -21.6 -151.1 -101.1 -131.2  154.7 -176.4 -147.8  178.2 -179.1  -53.3 0"
    string2 = "81 14.08  6.36 25.31  83.9 1 0.03 0.84 0.89  149.1   -8.0  -37.8  -39.9 -104.2  176.2 -179.6  137.8 -178.5 -173.3 -103.6 0"
    smiles = ['CC1=CC=C(C=C1)S(=O)(=O)C2=C(N=C(S2)C3=CC=C(C=C3)NC(=O)OCC4=CC=CC=C4)C']
    rep4 = representation.from_string(string1, smiles)
    rep5 = representation.from_string(string2, smiles)
    print(string1)
    print(string2)
    print(rep4.get_dist(rep5))

    #strings = [
    #"83 14.08  6.36 25.31  83.9 1 0.72 0.40 0.27  131.6  -17.0 -120.0  -83.8 -134.1 -174.5 -175.7 -168.8  173.9  178.0 -157.4 0",
    #"81 14.08  6.36 25.31  83.9 1 0.59 0.81 0.39 -117.8  -50.1  -95.3  -25.8  -80.6  164.7  155.9 -124.9 -159.2  178.6 -154.7 0",
    #"81 14.08  6.36 25.31  83.9 1 0.75 0.09 0.01  133.8  -19.5  -55.1  -86.7  -91.7 -175.0 -170.4 -176.8  173.3 -164.8  -58.4 0",
    #"81 14.08  6.36 25.31  83.9 1 0.72 0.44 0.01  135.2   27.5   97.2 -101.1 -105.1  -29.7 -169.7  -50.1  172.2 -173.1  131.6 0",
    #"82 14.00  6.34 25.26  83.6 1 0.21 0.08 0.54  146.0  -12.0   50.2  108.0  112.3 -166.3 -158.7  -35.5  172.3 -168.7  133.0 0",
    #"81 14.08  6.36 25.31  83.9 1 0.05 0.30 0.89  -68.2   41.2  148.8  -66.9  -85.0 -167.4  172.3 -166.2 -178.3  166.4  -45.9 0",
    #]

    #import pymatgen.analysis.structure_matcher as sm
    #matcher = sm.StructureMatcher(ltol=0.3, stol=0.3, angle_tol=10)
    #for i, string in enumerate(strings):
    #    print(str(i) + '  ' +string)
    #    rep4 = representation.from_string(string, smiles)
    #    pmg1 = rep4.to_pyxtal().to_pymatgen(); pmg1.remove_species('H')
    #    rep4.to_standard_setting()
    #    pmg2 = rep4.to_pyxtal().to_pymatgen(); pmg2.remove_species('H')
    #    print(i, rep4.to_string(), matcher.fit(pmg1, pmg2))
