import os
import numpy as np
from pyxtal import pyxtal
from pyxtal.symmetry import Group
from pyxtal.lattice import Lattice
from pyxtal.wyckoff_site import mol_site
 
class representation():
    """
    A class to handle the 1D representation of molecular crystal   
    Only works for Zprime = 1
    """

    def __init__(self, x, smiles=None):
        self.smiles = smiles
        self.x = x

    @classmethod
    def from_pyxtal(cls, struc):
        """
        Convert the crystal into 1D vector
        """
        symmetry = [struc.group.number, struc.diag]
        lat = struc.lattice.encode()
        vector = symmetry + lat
        for site in struc.mol_sites:
            vector += site.encode()
            smiles = [site.molecule.smile]
        x = vector
        return cls(x, smiles)
    
    def to_pyxtal(self, smiles=None):
        """
        convert to the pyxtal structure
        """
        if smiles is None:
            smiles = self.smiles

        if smiles is None:
            raise ValueError("no smiles")
        else:
            v = self.x
            if type(v) != list:
                v = list(v)
            v[0] = int(v[0])
            v[1] = bool(v[1])
            v[-1] = bool(v[-1])
            # symmetry
            struc = pyxtal(molecular=True)
            struc.group, struc.diag = Group(v[0]), v[1]
    
            # lattice
            count = 2
            ltype = struc.group.lattice_type
            if ltype == 'triclinic':
                a, b, c, alpha, beta, gamma = v[2], v[3], v[4], v[5], v[6], v[7]
                count += 6
            elif ltype == 'monoclinic':
                a, b, c, alpha, beta, gamma = v[2], v[3], v[4], 90, v[5], 90
                count += 4
            elif ltype == 'orthorhombic':
                a, b, c, alpha, beta, gamma = v[2], v[3], v[4], 90, 90, 90
                count += 3
            struc.lattice = Lattice.from_para(a, b, c, alpha, beta, gamma, ltype=ltype)
    
            # sites
            struc.numMols = [] 
            struc.molecules = []
            struc.mol_sites = [] 

            for i, smile in enumerate(smiles):
                dicts = {}
                dicts['smile'] = smile
                dicts['dim'] = 3
                dicts['PBC'] = [1, 1, 1]
                dicts['number'] = v[0]
                dicts['diag'] = int(v[1])
                dicts['index'] = 0
                dicts['lattice'] = struc.lattice.matrix
                dicts['lattice_type'] = ltype
                dicts['center'] = v[count:count+3]
                dicts['orientation'] = np.array(v[count+3:count+6])
                dicts['rotor'] = v[count+6:-1]
                dicts['reflect'] = int(v[-1])
                site = mol_site.from_1D_dicts(dicts)
                struc.mol_sites.append(site)
                struc.molecules.append(site.molecule)
                struc.numMols.append(site.wp.multiplicity)

            struc._get_formula()
            struc.source = '1D rep.'
            struc.valid = True

            return struc
    
    def to_string(self, time=None, eng=None, tag=None):
        """
        string representation
        """
        x = self.x
        strs = "{:3d} {:d} ".format(int(x[0]), int(x[1]))
        if int(x[0]) <=2:
            num = 11
        elif int(x[0]) <=15:
            num = 9
        elif int(x[0]) <=74:
            num = 8
            
        for v in x[2:num]:
            if v>99.9999:
                strs += "{:5.1f} ".format(v)      
            else:
                strs += "{:5.2f} ".format(v)      
        for v in x[num:-1]:
            strs += "{:6.1f} ".format(v)      
        strs += "{:d} ".format(int(x[-1]))      

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
    x = [14,False,11.43,6.49,11.19,83.31,0.77,0.57,0.53,48.55,24.31,145.94,-77.85,-4.40,170.86,False]
    rep0 = representation(x)
    print(rep0.to_string())
    rep1 = representation(x, smiles)
    xtal = rep1.to_pyxtal()
    print(xtal)
    rep2 = representation.from_pyxtal(xtal)
    print(rep2.to_pyxtal())
    print(rep2.to_string())

