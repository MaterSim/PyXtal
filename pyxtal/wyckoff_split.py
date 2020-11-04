import numpy as np
import pyxtal.symmetry as sym
from copy import copy
from pymatgen.core.operations import SymmOp
from random import choice

Wyc = {}
Wyc[197] = {"subgroup":[146, 146, 146, 146, 23], 
            "type": ['t', 't', 't', 't', 't'],
            "transformation": [
                               [[-1,0,1/2,0],[1,-1,1/2,0],[0,1,1/2,0]],
                               [[1,0,-1/2,0],[1,-1,1/2,0],[0,-1,-1/2,0]],
                               [[-1,0,1/2,0],[-1,1,-1/2,0],[0,-1,-1/2,0]],
                               [[1,0,-1/2,0],[-1,1,-1/2,0],[0,1,1/2,0]],
                               [[1,0,0,0],[0,1,0,0],[0,0,1,0]],
                              ],
            "relations": [
                        [['3a'], ['9b'], ['3a','9b'], ['9b','9b'], ['9b','9b'], ['9b','9b','9b','9b']],
                        [['3a'], ['9b'], ['3a','9b'], ['9b','9b'], ['9b','9b'], ['9b','9b','9b','9b']],
                        [['3a'], ['9b'], ['3a','9b'], ['9b','9b'], ['9b','9b'], ['9b','9b','9b','9b']],
                        [['3a'], ['9b'], ['3a','9b'], ['9b','9b'], ['9b','9b'], ['9b','9b','9b','9b']],
                        [['2a'], ['2b','2c','2d'], ['8k'], ['4e','4g','4i'], ['4f','4h','4j'], ['8k','8k','8k']]
                        ],
            }

Wyc[227] = {"subgroup":[216, 210, 203, 166, 166, 166, 166, 141, 141, 141, 227], 
            "type": ['t', 't', 't', 't', 't', 't', 't', 't', 't', 't', 'k'],
            "transformation": [[[1,0,0,1/8],[0,1,0,1/8],[0,0,1,1/8]], 
                               [[1,0,0,1/8],[0,1,0,1/8],[0,0,1,1/8]],
                               [[1,0,0,0],[0,1,0,0],[0,0,1,0]],
                               [[-1/2,0,1,0],[1/2,-1/2,1,0],[0,1/2,1,0]],
                               [[1/2,0,-1,1/4],[1/2,-1/2,1,0],[0,-1/2,-1,1/4]],
                               [[-1/2,0,1,0],[-1/2,1/2,-1,1/4],[0,-1/2,-1,1/4]],
                               [[1/2,0,-1,1/4],[-1/2,1/2,-1,1/4],[0,-1/2,-1,0]],
                               [[1/2,1/2,0,1/4],[-1/2,1/2,0,1/4],[0,0,1,0]],
                               [[0,0,1,0],[1/2,1/2,0,1/4],[-1/2,1/2,0,1/4]],
                               [[-1/2,1/2,0,1/4],[0,0,1,0],[1/2,1/2,0,1/4]],
                              ],
            "relations": [[['4a','4d'], ['4b','4c'], ['16e'], ['16e'], ['16e','16e'], ['24f','24g'], ['48h','48h'], ['96i'], ['96i','96i']],
                         [['8a'], ['8b'], ['16c'], ['16d'], ['32e'], ['48f'], ['96h'], ['48g', '48g'], ['96h', '96h']],
                         [['8a'], ['8b'], ['16c'], ['16d'], ['32e'], ['48f'], ['96g'], ['96g'], ['96g', '96g']],
                         [['6c'], ['6c'], ['3a','9d'], ['3b','9e'], ['6c','18h'], ['18h','18h'], ['18h','18h','36i'], ['18f','18g','36i'], ['36i','36i','36i','36i']],
                         [['6c'], ['6c'], ['3a','9d'], ['3b','9e'], ['6c','18h'], ['18h','18h'], ['18h','18h','36i'], ['18f','18g','36i'], ['36i','36i','36i','36i']],
                         [['6c'], ['6c'], ['3a','9d'], ['3b','9e'], ['6c','18h'], ['18h','18h'], ['18h','18h','36i'], ['18f','18g','36i'], ['36i','36i','36i','36i']],
                         [['6c'], ['6c'], ['3a','9d'], ['3b','9e'], ['6c','18h'], ['18h','18h'], ['18h','18h','36i'], ['18f','18g','36i'], ['36i','36i','36i','36i']],
                         [['4a'], ['4b'], ['8c'], ['8d'], ['16h'], ['8e','16g'], ['16h','32i'], ['16f','32i'], ['32i','32i','32i']],
                         [['4a'], ['4b'], ['8c'], ['8d'], ['16h'], ['8e','16g'], ['16h','32i'], ['16f','32i'], ['32i','32i','32i']],
                         [['4a'], ['4b'], ['8c'], ['8d'], ['16h'], ['8e','16g'], ['16h','32i'], ['16f','32i'], ['32i','32i','32i']]],
            }


class wyckoff_split:
    """
    Class for performing wyckoff split between two space groups.
    Essentially, this code is to look for the database from the
    internationally crystallographic table and find the group-subgroup
    relations

    Args:
        G: int (1-230), number of super space group
        H: int (1-230), number of super space group
        w1: string ("4a") or integer (1)
    """


    def __init__(self, G=197, H=146, wp1=[0, 1]):
        
        self.G = sym.Group(G)  # Group object
        self.H = sym.Group(H)  # Group object
           
        id_lists = []
        for wp in wp1:
            if type(wp) == int:
                id_lists.append(wp)
            else:
                id_lists.append(sym.index_from_letter(wp[-1], self.G))
        self.wp1_indices = id_lists
        self.wp1_lists = [self.G[id] for id in id_lists] # a WP object
        self.query_wp2()

        self.G_orbits = []
        self.H_orbits = []
        for i, wp1 in enumerate(self.wp1_lists):
            G_orbits, H_orbits = self.split(wp1, self.wp2_lists[i])
            self.G_orbits.append(G_orbits)
            self.H_orbits.append(H_orbits)

    def query_wp2(self):
        """
        query the wp2 and transformation matrix from the given {G, H, wp1}
        """
        ids = []
        wyc = Wyc[self.G.number]
        ids = [id for id in range(len(wyc['subgroup'])) if wyc['subgroup'][id]==self.H.number]
        id = choice(ids)
        trans = np.array(wyc['transformation'][id])  #np.array([[-1.,0.,.5,0.],[1.,-1.,.5,0.],[0.,1.,.5,0.],[0.,0.,0.,1.]])
        self.R = np.zeros([4,4])
        self.R[:3,:3] += trans[:3,:3]
        self.R[3,3] = 1
        self.inv_R = np.linalg.inv(self.R)
        self.t = trans[:,3]
        subgroup_relations = wyc['relations'][id]
        subgroup_relations.reverse()
        #subgroup_relations = [['9b','9b','9b','9b'],
        #                      ['9b','9b'],
        #                      ['9b','9b'],
        #                      ['3a','9b'],
        #                      ['9b'],
        #                      ['3a']]
        print(self.wp1_indices)
        wp2_lists = []
        for wp1_index in self.wp1_indices:
            wp2_list = []
            for letter in subgroup_relations[wp1_index]:
                id = sym.index_from_letter(letter[-1], self.H)
                wp2_list.append(self.H[id])
            wp2_lists.append(wp2_list)
        self.wp2_lists = wp2_lists
    
    def split(self, wp1, wp2_lists):
        """
        split the generators in w1 to different w2s
        """
                
        # wyckoff objects
        wp1_generators_visited = []
        wp1_generators = [wp.as_dict()['matrix'] for wp in wp1]
        
        # convert them to numpy array
        for generator in wp1_generators:
            generator = np.array(generator)
            generator[:,3] -= np.floor(generator[:,3])

        G_orbits = []
        H_orbits = []
        factor = np.linalg.det(self.R)

        for wp2 in wp2_lists:
            # try all generators here
            for gen in wp1_generators:
                good_generator = True
                trans_generator = np.matmul(self.inv_R, gen)
                
                G1_orbits = []
                H1_orbits = []
                strs = []
                for i, wp in enumerate(wp2):
                    new_basis_orbit = np.matmul(wp.as_dict()['matrix'], trans_generator)
                    old_basis_orbit = np.matmul(self.R, new_basis_orbit).round(3)
                    
                    tmp = copy(old_basis_orbit)
                    tmp[:,3] -= np.floor(tmp[:,3])
                    if i==0 and (tmp.tolist() in wp1_generators_visited):
                        good_generator = False
                        break
                    #str1 = SymmOp(new_basis_orbit).as_xyz_string()
                    #str2 = SymmOp(old_basis_orbit).as_xyz_string()
                    #print("{:32s} --> {:32s}".format(str2, str1))
                    new_basis_orbit[3,:3] += self.t
                    G1_orbits.append(old_basis_orbit)        
                    H1_orbits.append(new_basis_orbit)        
                
                # remove duplicates due to peridoic boundary conditions
                if good_generator:
                    temp=[]
                    for gen in G1_orbits:
                        gen[:,3] -= np.floor(gen[:,3])
                        gen_list = gen.tolist()
                        if gen.tolist() not in temp:
                            temp.append(gen_list)
                    if len(temp)*factor == len(wp2):           
                        wp1_generators_visited.extend(temp)
                        G1_orbits = [SymmOp(orbit) for orbit in G1_orbits]
                        H1_orbits = [SymmOp(orbit) for orbit in H1_orbits]
                        G_orbits.append(G1_orbits)
                        H_orbits.append(H1_orbits)
                        break
        return G_orbits, H_orbits
        
    def __str__(self):
        s = "Wycokff split from {:d} to {:d}\n".format(self.G.number, self.H.number)
        for i, wp1 in enumerate(self.wp1_lists):
            s += "Origin: {:d}{:s}\n".format(wp1.multiplicity, wp1.letter)
            s += "After spliting\n"
        
            for j, wp2 in enumerate(self.wp2_lists[i]):
                s += "{:d}{:s}\n".format(wp2.multiplicity, wp2.letter)
                for g_orbit, h_orbit in zip(self.G_orbits[i][j], self.H_orbits[i][j]):
                    s += "{:30s} -> {:30s}\n".format(g_orbit.as_xyz_string(), h_orbit.as_xyz_string())
        return s
        
    def __repr__(self):
        return str(self)

        
if __name__ == "__main__":

    from pyxtal.crystal import random_crystal
    from pyxtal.operations import apply_ops
    from ase import Atoms
    import pymatgen.analysis.structure_matcher as sm
    from spglib import get_symmetry_dataset
    from pymatgen.io.ase import AseAtomsAdaptor
    #sites = ['24f','6b']
    #G, H, fac = 197, 23, 2
    sites = ['8a', '32e']
    G, H, fac = 227, 216, 4
    numIons = int(sum([int(i[:-1]) for i in sites])/fac)
    print(numIons)
    C = random_crystal(G, ['C'], [numIons], sites=[sites])
    spg1 = get_symmetry_dataset(C.to_ase(), symprec=1e-4)['international']

    splitter = wyckoff_split(G=G,H=H,wp1=sites)
    print(splitter)
    lat1 = np.dot(C.lattice.matrix, splitter.R[:3,:3].T)
    pos1 = None
    for i, site in enumerate(C.atom_sites):
        pos = site.position
        for ops in splitter.H_orbits[i]:
            pos0 = pos + 0.05*(np.random.sample(3) - 0.5)
            #print(pos0)
            pos_tmp = apply_ops(pos0, ops)
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
