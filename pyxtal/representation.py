import os
import numpy as np
from pyxtal.symmetry import Group
from pyxtal.lattice import Lattice
from pyxtal.wyckoff_site import mol_site
from pyxtal.molecule import find_id_from_smile
 
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
    def from_pyxtal(cls, struc):
        """
        Initialize 1D rep. from the pyxtal object

        Args:
            struc: pyxtal object
        """
        symmetry = [struc.group.number, struc.diag]
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
        g, diag = int(inputs[0]), int(inputs[1])
        if g <= 2:
            n_cell = 8
        elif g <= 15:
            n_cell = 6
        elif g <= 74:
            n_cell = 5
        elif g <= 194:
            n_cell = 4
        else:
            n_cell = 3 #cubic
        cell = [g, diag] + inputs[2:n_cell]
        
        x = [cell]
        n_site = int(inputs[n_cell])
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
                    n_torsion = len(find_id_from_smile(smile))
                    n_mol = 7 + n_torsion
                inputs[n_cell+n_mol-1] = int(inputs[n_cell+n_mol-1])
                x.append(inputs[n_cell:n_cell+n_mol])
                n_cell += n_mol
        return cls(x, smiles)
 
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
        struc.group, number, struc.diag = Group(v[0]), v[0], v[1]
    
        # lattice
        ltype = struc.group.lattice_type
        if ltype == 'triclinic':
            a, b, c, alpha, beta, gamma = v[2], v[3], v[4], v[5], v[6], v[7]
        elif ltype == 'monoclinic':
            a, b, c, alpha, beta, gamma = v[2], v[3], v[4], 90, v[5], 90
        elif ltype == 'orthorhombic':
            a, b, c, alpha, beta, gamma = v[2], v[3], v[4], 90, 90, 90
        elif ltype == 'tetragonal':
            a, b, c, alpha, beta, gamma = v[2], v[2], v[3], 90, 90, 90
        elif ltype == 'hexagonal':
            a, b, c, alpha, beta, gamma = v[2], v[2], v[3], 90, 90, 120
        else:
            a, b, c, alpha, beta, gamma = v[2], v[2], v[2], 90, 90, 90
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
                dicts['number'] = number
                dicts['diag'] = struc.diag
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
        strs = "{:3d} {:d} ".format(int(x[0][0]), int(x[0][1]))

        # data for cell
        if x[0][0] <= 74:
            num = 5
        elif x[0][0] <=194:
            num = 4
        else:
            num = 3

        for c in x[0][2:num]:
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



if __name__ == "__main__":

    #aspirin
    smiles = ['CC(=O)OC1=CC=CC=C1C(=O)O']
    x = [[14,False,11.43,6.49,11.19,83.31],[0.77,0.57,0.53,48.55,24.31,145.94,-77.85,-4.40,170.86,False]]
    #rep0 = representation(x, smiles)
    #print(rep0.to_string())
    rep1 = representation(x, smiles)
    xtal = rep1.to_pyxtal()
    print(xtal)
    rep2 = representation.from_pyxtal(xtal)
    print(rep2.to_pyxtal())
    print(rep2.to_string())
    string = "14 0 11.43  6.49 11.19 83.31 1  0.77  0.57  0.53 48.55 24.31 145.9 -77.85 -4.40 170.9 0"
    rep3 = representation.from_string(string, smiles)
    print(rep3.to_pyxtal())
