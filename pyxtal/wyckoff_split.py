import numpy as np
import pyxtal.symmetry as sym
from copy import copy
from pymatgen.core.operations import SymmOp

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


    def __init__(self, G=197, H=146, wp1=1):
        
        self.G = sym.Group(G)  # Group object
        self.H = sym.Group(H)  # Group object
        if type(wp1) == int:
            id = wp1
        else:
            id = sym.index_from_letter(wp1[-1], self.H)
        self.wp1 = self.G[id] # a WP object
        self.query()
        #print(str(self.wp1)) 
        #print(self.wp2_lists) 
        self.split()
    
    def query(self):
        """
        query the wp2 and transformation matrix from the given {G, H, wp1}
        """
        self.T = np.array([[-1.,0.,.5,0.],[1.,-1.,.5,0.],[0.,1.,.5,0.],[0.,0.,0.,1.]])
        self.inv_T = np.linalg.inv(self.T)
        subgroup_relations = [['9b','9b','9b','9b'],
                              ['9b','9b'],
                              ['9b','9b'],
                              ['3a','9b'],
                              ['9b'],
                              ['3a']]
        
        wp1_index = self.wp1.index
        wp2_lists = []
        for letter in subgroup_relations[wp1_index]:
            id = sym.index_from_letter(letter[-1], self.H)
            wp2_lists.append(self.H[id])
        self.wp2_lists = wp2_lists

    
    def split(self):
                
        # wyckoff objects
        wp1_generators_visited = []
        wp1_generators = [wp.as_dict()['matrix'] for wp in self.wp1]
        
        # convert them to numpy array
        for generator in wp1_generators:
            generator = np.array(generator)
            generator[:,3] -= np.floor(generator[:,3])

        G_orbits = []
        # split the generators in w1 to different w2s
        
        
        for wp2 in self.wp2_lists:
            # try all generators here
            for gen in wp1_generators:
                
                good_generator = True
                trans_generator = np.matmul(self.inv_T, gen)
                
                G1_orbits = []
                strs = []
                for i, wp in enumerate(wp2):
                    new_basis_orbit = np.matmul(wp.as_dict()['matrix'], trans_generator)
                    old_basis_orbit = np.matmul(self.T, new_basis_orbit).round(3)
                    
                    tmp = copy(old_basis_orbit)
                    tmp[:,3] -= np.floor(tmp[:,3])
                    
                    if i==0 and tmp.tolist() in wp1_generators_visited:
                        good_generator = False
                        break
                        
                    #str1 = SymmOp(new_basis_orbit).as_xyz_string()
                    #str2 = SymmOp(old_basis_orbit).as_xyz_string()
                    #print("{:32s} --> {:32s}".format(str2, str1))
                    G1_orbits.append(old_basis_orbit)        
                
                # remove duplicates due to peridoic boundary conditions
                if good_generator:
                    temp=[]
                    for gen in G1_orbits:
                        gen[:,3] -= np.floor(gen[:,3])
                        gen_list = gen.tolist()
                        if gen.tolist() not in temp:
                            temp.append(gen_list)
            
                    wp1_generators_visited.extend(temp)
                    G1_orbits = [SymmOp(orbit) for orbit in G1_orbits]
                    G_orbits.append(G1_orbits)
                    break
                    
        self.G_orbits = G_orbits
        
    def __str__(self):
        s = "Wycokff split from {:d} to {:d}\n".format(self.G.number, self.H.number)
        s += "Origin: {:d}{:s}\n".format(self.wp1.multiplicity, self.wp1.letter)
        s += "After spliting\n"
        
        for i, wp2 in enumerate(self.wp2_lists):
            s += "{:d}{:s}\n".format(wp2.multiplicity, wp2.letter)
        
            for orbit in self.G_orbits[i]:
                s += orbit.as_xyz_string() + '\n'
        return s
        
    def __repr__(self):
        return str(self)

        
if __name__ == "__main__":
    
    splitter = wyckoff_split()
    print(splitter)
