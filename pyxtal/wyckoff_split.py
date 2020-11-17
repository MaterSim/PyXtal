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


    def __init__(self, G=197, idx=None, wp1=[0, 1], group_type='t'):
        
        self.G = sym.Group(G)  # Group object
        if group_type == 't':
            self.wyc = self.G.get_max_t_subgroup()
        else:
            self.wyc = self.G.get_max_k_subgroup()
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
        #print(G, idx, len(self.wyc['subgroup']))
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
            if group_type == 't':
                G1_orbits, G2_orbits = self.split_t(wp1, self.wp2_lists[i])
            else:
                G1_orbits, G2_orbits = self.split_k(wp1, self.wp2_lists[i])
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
    
    def split_t(self, wp1, wp2_lists):
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
                        #print("wp1", SymmOp(gen).as_xyz_string(), in_lists(tmp, wp1_generators_visited), in_lists(tmp, wp1_generators))
                        #print("sgb", SymmOp(new_basis_orbit).as_xyz_string())
                        #print("gb", SymmOp(tmp).as_xyz_string())
                        #for w in wp1_generators:
                        #    print(SymmOp(w).as_xyz_string())
                        if not in_lists(tmp, wp1_generators_visited) and in_lists(tmp, wp1_generators):
                        #if in_lists(tmp, wp1_generators):
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
                    if int(len(temp)*factor) >= len(wp2):           
                        wp1_generators_visited.extend(temp)
                        g1_orbits = [SymmOp(orbit) for orbit in g1_orbits]
                        g2_orbits = [SymmOp(orbit) for orbit in g2_orbits]
                        G1_orbits.append(g1_orbits)
                        G2_orbits.append(g2_orbits)
                        #print("adding unique generators", len(g1_orbits), len(wp2), int(len(temp)*factor))
                        break
            #print("EEEEEE", len(g1_orbits), len(wp2))
            self.check_orbits(g1_orbits, wp1, wp2, wp2_lists)

        return G1_orbits, G2_orbits
        
    def split_k(self, wp1, wp2_lists):
        """
        split the generators in w1 to different w2s for k_subgroups
        """
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
                trans_generator = np.matmul(self.inv_R, gen)
                
                g1_orbits = []
                g2_orbits = []
                strs = []
                for i, wp in enumerate(wp2):
                    new_basis_orbit = np.matmul(wp.as_dict()['matrix'], trans_generator)
                    old_basis_orbit = np.matmul(self.R, new_basis_orbit).round(3)
                    tmp = deepcopy(old_basis_orbit)
                    tmp[3,:] = [0, 0, 0, 1]
                    if factor <= 4:
                        if i==0: 
                            if not in_lists(tmp, wp1_generators_visited, PBC=False):
                                good_generator = True
                            else:
                                break
                    else:
                        good_generator = True
                    g1_orbits.append(old_basis_orbit)        
                    g2_orbits.append(new_basis_orbit)        
                if good_generator:
                    if len(g1_orbits) >= len(wp2):           
                        wp1_generators_visited.extend(g1_orbits)
                        g1_orbits = [SymmOp(orbit) for orbit in g1_orbits]
                        g2_orbits = [SymmOp(orbit) for orbit in g2_orbits]
                        G1_orbits.append(g1_orbits)
                        G2_orbits.append(g2_orbits)
                        break
            self.check_orbits(g1_orbits, wp1, wp2, wp2_lists)
        return G1_orbits, G2_orbits
 
    def check_orbits(self, g1_orbits, wp1, wp2, wp2_lists):
        if len(g1_orbits) < len(wp2):
            s1 = str(wp1.multiplicity)+wp1.letter
            s2 = ""
            for wp2 in wp2_lists:
                s2 += str(wp2.multiplicity)+wp2.letter
                s2 += ', '
            g, h = self.G.number, self.H.number
            print("Error between {:d}[{:s}] -> {:d}[{:s}]".format(g, s1, h, s2))
            print(self.R)
            print(g1_orbits)
            raise ValueError("Cannot find the generator for wp2")

    def __str__(self):
        s = "Wycokff split from {:d} to {:d}\n".format(self.G.number, self.H.number)
        for i, wp1 in enumerate(self.wp1_lists):
            s += "Origin: {:d}{:s}\n".format(wp1.multiplicity, wp1.letter)
            s += "After spliting\n"
        
            for j, wp2 in enumerate(self.wp2_lists[i]):
                s += "{:d}{:s}\n".format(wp2.multiplicity, wp2.letter)
                g1s, g2s, Hs = self.G1_orbits[i][j], self.G2_orbits[i][j], self.H_orbits[i][j]
                for g1_orbit, g2_orbit, h_orbit in zip(g1s, g2s, Hs):
                    s += "{:30s} -> {:30s} -> {:30s}\n".format(g1_orbit.as_xyz_string(), \
                                                               g2_orbit.as_xyz_string(), \
                                                               h_orbit.as_xyz_string())
        return s
        
    def __repr__(self):
        return str(self)


def in_lists(mat1, mat2, eps=1e-4, PBC=True):
    if len(mat2) == 0:
        return False
    else:
        for mat in mat2:
            if np.array_equal(mat[:3,:3], mat1[:3,:3]):
                diffs = np.abs(mat[:3,3] - mat1[:3,3])
                if PBC:
                    diffs -= np.floor(diffs)
                #print("diffs", diffs)
                if (diffs*diffs).sum() < 1e-2:
                    return True
        return False

        
if __name__ == "__main__":
    import pymatgen.analysis.structure_matcher as sm
    from random import choice
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    from pyxtal.crystal import random_crystal
    
    while True:
        sg = 21 #choice(range(141,231))
        s1 = random_crystal(sg, ['B'], [3], sites=[['4f','2d']])
        if s1.valid:
            break
    print(s1)
    pmg_s1 = s1.to_pymatgen()
    sga1 = SpacegroupAnalyzer(pmg_s1, symprec=1e-4).get_space_group_symbol()
    #s2s = s1.subgroup(group_type='k', eps=0)
    s2s = s1.subgroup(group_type='k', eps=0)#, idx=[5])
    for i, s2 in enumerate(s2s):
        pmg_s2 = s2.to_pymatgen()
        try:
            sga2 = SpacegroupAnalyzer(pmg_s2, symprec=1e-4).get_space_group_symbol()
            print(i, sga1, sga2, sm.StructureMatcher().fit(pmg_s1, pmg_s2))
        except:
            print("something is wrong here")
            print(s2)
            import sys; sys.exit()
#s1 = s2

