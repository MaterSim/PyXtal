import numpy as np
import pyxtal.symmetry as sym
from copy import deepcopy
from pymatgen.core.operations import SymmOp
from random import choice

class wyckoff_split:
    """
    Class for performing wyckoff split between two space groups.
    Essentially, this code is to look for the database from the
    internationally crystallographic table and find the group-subgroup
    relations

    Args:
        G: int (1-230), number of super space group
        w1: string ("4a") or integer (1)
    """


    def __init__(self, G=197, idx=None, wp1=[0, 1]):
        
        self.G = sym.Group(G)  # Group object
        self.wyc = self.G.get_max_t_subgroup()
        id_lists = []
        for wp in wp1:
            if type(wp) == int:
                id_lists.append(wp)
            else:
                id_lists.append(sym.index_from_letter(wp[-1], self.G))
        self.wp1_indices = id_lists
        self.wp1_lists = [self.G[id] for id in id_lists] # a WP object

        # choose a random spliting option if idx is not specified
        if idx is None:
            ids = [id for id in range(len(self.wyc['subgroup']))]
            idx = choice(ids)

        H = self.wyc['subgroup'][idx]
        self.H = sym.Group(H)  # Group object

        self.parse_wp2(idx)

        if (self.G.lattice_type == self.H.lattice_type):
            self.valid_split = False
            for wps in self.wp2_lists:
                for wp in wps:
                    rotation = np.array(wp[0].as_dict()['matrix'])[:3,:3]
                    if np.linalg.matrix_rank(rotation) > 0:
                        self.valid_split = True
                        break
        else:
            self.valid_split = True
        
        #if self.valid_split:
        self.G1_orbits = []
        self.G2_orbits = []
        self.H_orbits = []
        for i, wp1 in enumerate(self.wp1_lists):
            self.H_orbits.append([wp2.ops for wp2 in self.wp2_lists[i]])
            G1_orbits, G2_orbits = self.split(wp1, self.wp2_lists[i])
            self.G1_orbits.append(G1_orbits)
            self.G2_orbits.append(G2_orbits)


    def parse_wp2(self, idx):
        """
        query the wp2 and transformation matrix from the given {G, H, wp1}
        """
        #print("trans", idx)
        #print(self.wyc['transformation'])
        trans = self.wyc['transformation'][idx]
        subgroup_relations = self.wyc['relations'][idx]
        subgroup_relations = [ele for ele in reversed(subgroup_relations)] 

        self.R = np.zeros([4,4])
        self.R[:3,:3] += trans[:3,:3]
        self.R[3,3] = 1
        self.inv_R = np.linalg.inv(self.R)
        inv_t = np.dot(self.inv_R[:3,:3], trans[:,3].T)
        self.inv_R[:3,3] = -inv_t.T
        self.R[:3,3] = trans[:3,3]


        wp2_lists = []
        for wp1_index in self.wp1_indices:
            wp2_list = []
            for letter in subgroup_relations[wp1_index]:
                id = sym.index_from_letter(letter[-1], self.H)
                wp2_list.append(self.H[id])
            wp2_lists.append(wp2_list)
        self.wp2_lists = wp2_lists
        #import sys; sys.exit()
    
    def split(self, wp1, wp2_lists):
        """
        split the generators in w1 to different w2s
        """
        #print(wp1)
        # wyckoff objects
        wp1_generators_visited = []
        wp1_generators = [np.array(wp.as_dict()['matrix']) for wp in wp1]
        
        # convert them to numpy array
        for generator in wp1_generators:
            generator = np.array(generator)

        G1_orbits = []
        G2_orbits = []
        factor = max([1,np.linalg.det(self.R)])

        for wp2 in wp2_lists:
            #print(wp2)
            #import sys; sys.exit()
            # try all generators here
            for gen in wp1_generators:
                good_generator = False
                #QZ: temporary solution, Needs to be fixed later
                if gen[0,3] == 1/4 and gen[1,3] == 3/4:
                    gen[0,3] -=1
                trans_generator = np.matmul(self.inv_R, gen)
                #print(trans_generator)
                
                g1_orbits = []
                g2_orbits = []
                strs = []
                for i, wp in enumerate(wp2):
                    new_basis_orbit = np.matmul(wp.as_dict()['matrix'], trans_generator)
                    #print(wp.as_dict()['matrix'])
                    #print(new_basis_orbit)
                    #import sys; sys.exit()
                    old_basis_orbit = np.matmul(self.R, new_basis_orbit).round(3)
                    #old_basis_orbit[3,:] = [0, 0, 0, 1]
                    tmp = deepcopy(old_basis_orbit)
                    tmp[3,:] = [0, 0, 0, 1]
                    if i==0: 
                        #print("wp1", SymmOp(gen).as_xyz_string())
                        #print("sgb", SymmOp(new_basis_orbit).as_xyz_string())
                        #print("gb", SymmOp(tmp).as_xyz_string())
                        #for w in wp1_generators:
                        #    print(SymmOp(w).as_xyz_string())
                        #print(in_lists(tmp, wp1_generators_visited), in_lists(tmp, wp1_generators))
                        if not in_lists(tmp, wp1_generators_visited) and in_lists(tmp, wp1_generators):
                            good_generator = True
                            #print("good_gener")
                        else:
                            break
                    # to consider PBC
                    g1_orbits.append(old_basis_orbit)        
                    g2_orbits.append(new_basis_orbit)        
                #print(g1_orbits)
                if good_generator:
                    temp=[]
                    # remove duplicates due to peridoic boundary conditions
                    for gen in g1_orbits:
                        if not in_lists(gen, temp):
                            temp.append(gen)
                    #print("adding unique generators", len(temp))
                    if int(len(temp)*factor) == len(wp2):           
                        wp1_generators_visited.extend(temp)
                        g1_orbits = [SymmOp(orbit) for orbit in g1_orbits]
                        g2_orbits = [SymmOp(orbit) for orbit in g2_orbits]
                        G1_orbits.append(g1_orbits)
                        G2_orbits.append(g2_orbits)
                        break

            if len(g1_orbits) < len(wp2):
                s1 = str(wp1.multiplicity)+wp1.letter
                s2 = ""
                for wp2 in wp2_lists:
                    s2 += str(wp2.multiplicity)+wp2.letter
                    s2 += ', '
                print("Error in split between {:d}[{:s}] -> {:d}[{:s}]".format(self.G.number, s1, self.H.number, s2))
                print(self.R)
                print(g1_orbits)
                raise ValueError("Cannot find the generator for wp2")
        return G1_orbits, G2_orbits
        
    def __str__(self):
        s = "Wycokff split from {:d} to {:d}\n".format(self.G.number, self.H.number)
        for i, wp1 in enumerate(self.wp1_lists):
            s += "Origin: {:d}{:s}\n".format(wp1.multiplicity, wp1.letter)
            s += "After spliting\n"
        
            for j, wp2 in enumerate(self.wp2_lists[i]):
                s += "{:d}{:s}\n".format(wp2.multiplicity, wp2.letter)
                for g1_orbit, g2_orbit, h_orbit in zip(self.G1_orbits[i][j], self.G2_orbits[i][j], self.H_orbits[i][j]):
                    s += "{:30s} -> {:30s} -> {:30s}\n".format(g1_orbit.as_xyz_string(), \
                                                               g2_orbit.as_xyz_string(), \
                                                               h_orbit.as_xyz_string())
        return s
        
    def __repr__(self):
        return str(self)


def in_lists(mat1, mat2, eps=1e-4):
    if len(mat2) == 0:
        return False
    else:
        for mat in mat2:
            if np.array_equal(mat[:3,:3], mat1[:3,:3]):
                diffs = np.abs(mat[:3,3] - mat1[:3,3])
                diffs -= np.floor(diffs)
                #print("diffs", diffs)
                if (diffs*diffs).sum() < 1e-2:
                    return True
        return False

        
if __name__ == "__main__":

    from pyxtal.crystal import random_crystal
    from pyxtal.operations import apply_ops
    from ase import Atoms
    import pymatgen.analysis.structure_matcher as sm
    from spglib import get_symmetry_dataset
    from pymatgen.io.ase import AseAtomsAdaptor
    #sites = ['24f','6b']
    #G, H, fac = 197, 23, 2
    sites = ['32e']
    #sites = ['8a']
    G, H, fac = 227, 166, 4
    numIons = int(sum([int(i[:-1]) for i in sites])/fac)
    C = random_crystal(G, ['C'], [numIons], sites=[sites])
    spg1 = get_symmetry_dataset(C.to_ase(), symprec=1e-4)['international']

    splitter = wyckoff_split(G=G,H=H,wp1=sites)
    print(splitter)
    lat1 = np.dot(C.lattice.matrix, splitter.R[:3,:3].T)
    pos1 = None
    for i, site in enumerate(C.atom_sites):
        pos = site.position
        for ops1, ops2 in zip(splitter.G2_orbits[i], splitter.H_orbits[i]):
            pos0 = pos + 0.05*(np.random.sample(3) - 0.5)
            #print(pos0)
            pos_tmp = apply_ops(pos0, ops1)
            pos_tmp = apply_ops(pos_tmp[0], ops2)
            if pos1 is None:
                pos1 = pos_tmp
            else:
                pos1 = np.vstack((pos1, pos_tmp))

    C1 = Atoms(["C"]*len(pos1), scaled_positions=pos1, cell=lat1, pbc=[1,1,1])
    C1.write("1.vasp", format='vasp', vasp5=True, direct=True)
    C.to_ase().write("0.vasp", format='vasp', vasp5=True, direct=True)
    pmg_s1 = AseAtomsAdaptor.get_structure(C1)
    spg2 = get_symmetry_dataset(C1, symprec=1e-4)['international']
    print(spg1, spg2, sm.StructureMatcher().fit(pmg_s1, C.to_pymatgen()))
