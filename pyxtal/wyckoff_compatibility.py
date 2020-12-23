import pyxtal.symmetry as sym
import numpy as np
from copy import deepcopy


class wyckoff_compatibility():
    #see the example at the bottom for how the the atoms and positions are supposed to be organized.

    def __init__(self, H, group_type='t', idxs=None, atoms=None, positions=None):
        wyc_supergroups=sym.Group(H).get_min_supergroup(group_type)
        self.results={}
        self.supergroups=[]
        if idxs==None:
            for i,G in enumerate(wyc_supergroups['supergroup']):
                self.results[G]=self.compatibility_calculator(G,H,wyc_supergroups['relations'][i],atoms,positions)
                self.supergroups.append(G)
        else:
            for idx in idxs:
                G=wyc_supergroups['supergroup'][idx]
                self.results[G]=self.compatibility_calculator(G,H,wyc_supergroups['relations'][idx],atoms,positions)
                self.supergroups.append(G)




    #The whole functionality of the compatiblity splitter is to find out if the G supergroups
    #have a flexible enough splitting  scheme to be able to let the atoms of
    #subgroup fit without losing or gaining atoms, or leaving holes in the material.


    #In my example material from the literature, I have 4 flourine, 1 Na, and 1 Sb.
    #They all have their wyckoff positions that they occupy in the IT:173 group, so
    #so this splitter is trying to find if there are openings in the supergroup wyckoff
    #positions that will be able to fit this structure.
    def compatibility_calculator(self,G,H,relations,atoms,positions):

        H=sym.Group(H)
        G=sym.Group(G)
        results={}
        wyc_position_list=[(str(x.multiplicity)+x.letter) for x in G]
        wyc_position_list.reverse()
        atom_basis=np.unique(atoms) #['Na','Sb','F']
        atom_positions=[]
        for atom in atom_basis:
            atom_positions.append([positions[i] for i,x in enumerate(atoms) if x==atom])
            #[ ['2b'] , ['6c']  , ['6c','6c','6c','2b']  ]


        species_list=[]
        atom_number_list=[]
        positions_list=[]
        good_splittings_list=[]
        compatible_split=True

        #The rest of this is a lot of integer math. The goal is to find all the
        #integer combinations of supergroup wyckoff positions such that the union of
        #all their splitting schemes will be exactly what's in the atom_positions list.


        for i,species in enumerate(atom_basis):#each element is solved one at a time

            position_basis=np.unique(atom_positions[i])
            position_counts=[atom_positions[i].count(x) for x in position_basis]
            possible_wyc_indexes=[]


            total_units=0
            for j,x in enumerate(position_basis):
                total_units+=int(x[:-1])*position_counts[j]
            #a criteria that is used to eliminate a lot of possible
            #combinations later on is that the total sum multiplicity
            #of all the positions occupied by a species is fixed. for
            #example, when solving fluorine, this retricts the possible position combinations to
            #only the union of supergroup positions whose splittings add up to 20.



            for j,split in enumerate(relations):
                if np.all([x in position_basis for x in split]):
                    possible_wyc_indexes.append(j)
                    #collects all inidices of the supergroup positions that splits into
                    #only the positions that the species is allowed to occupy


            #Here, i try to employ some matrix organization. I am using n sized lists as essentially vectors to represent
            #the possible combinations of split positions that each of the possible_wyc_positions.
            #For example, when solving fluorine, the "position_basis" vectors would be [6c,2b], and for something like
            #The 12i position in the 182->173 splitting that goes 12i->[6c,6c] would be represented as a [2,0]. All the
            #position_blocks holds the vectors. All the possible supergroup positions are turned into 2 item
            #arrays and stored into the position_blocks variable. Then the task is just to find all the integer combinations
            #of those vectors that will add up to [3,1] which is stored in the position_counts variables.
            possible_wyc_positions=[wyc_position_list[x] for x in possible_wyc_indexes]
            position_blocks=[np.array([relations[j].count(pos) for pos in position_basis]) for j in possible_wyc_indexes]
            position_block_units=[sum([int(x[:-1])*block[j] for j,x in enumerate(position_basis)]) for block in position_blocks]
            # the position_block_units stores the total sum multiplicty offered by the splitting from the available
            #supergroup wyckoff positions. The the integer combinations of these numbers must add up to the total_units number above
            #as well for a combination to be valid




            #This is just a brute force search and test method for finding the combinations of splittings that
            #satisifies both the unit sum and vector addition requirement above.

            #It systematically creates all the ways you can distribute integers into lists of length n, where n is the number
            #of possible wyckoff positions we have to work with.

            #at the start, combo_storage holds all the possible ways to distribute one integer into an n sized list.
            #Then, each list is dotted with the stored sum multiplicites. If a combo gives a number that is too low, then it's sent back
            #combo_storage to have more integers distributed. If it's too high, then the combo is deleted for going too far.
            #When a list gives the exact required sum multiplicity, then it's tested for the linear combination of position independent_vectors
            #requirement. The n sized vectors get turned back into combos of supergroup wyckoff letters once it passes.

            #I believe that this scheme should be able to grow and traverse all the integer combinations of wyckoff positions possible until
            #they eventually either land on the right number, or gets discarded for overshooting.
            bins=len(position_block_units)
            combo_storage=[np.zeros(bins)]
            good_list=[]
            while len(combo_storage)!=0:
                holder=[]
                for j,x in enumerate(combo_storage):
                    for k in range(bins):
                        trial=np.array(deepcopy(x))
                        trial[k]+=1
                        sum_units=np.dot(trial,position_block_units)
                        if sum_units>total_units:
                            continue
                        elif sum_units<total_units:
                            holder.append(trial)
                        else:
                            tester=np.zeros(len(position_counts))
                            for l,z in enumerate(trial):
                                tester+=z*position_blocks[l]
                            if np.all(tester==position_counts):
                                supergroup_positions=[]
                                for l,number in enumerate(trial):
                                    if number==0:
                                        continue
                                    elif number==1:
                                        supergroup_positions.append(possible_wyc_positions[l][-1])
                                    else:
                                        supergroup_positions.append((str(int(number))+possible_wyc_positions[l][-1]))
                                if supergroup_positions not in good_list:
                                    good_list.append(supergroup_positions)

                combo_storage=holder

            if len(good_list)==0:
                good_list=[None]
                compatible_split=False
            species_list.append(species)
            atom_number_list.append(atoms.count(species))
            positions_list.append(atom_positions[i])
            good_splittings_list.append(good_list)

        results={'minimum_supergroup':G,
                 'species':species_list,
                 'species_number':atom_number_list,
                 'species_positions':positions_list,
                 'compatible_supergroup_positions':good_splittings_list,
                 'compatible_split':compatible_split
                }
        return results










if __name__ == "__main__":
    H=173
    atoms=['Sb','Na','F','F','F','F']
    positions=['6c','2b','6c','6c','6c','2b']
    wc=wyckoff_compatibility(H,group_type='t',atoms=atoms,positions=positions)



    s=wc.supergroups[0]
    print('MATERIAL: Na-Sb3-F10', "Atomic Species:",wc.results[s]['species'], 'Group Number:',H)
    print("Atomic multiplicities of", wc.results[s]['species_number']," and with wyckoff positions ", wc.results[s]['species_positions'])




    for G in wc.supergroups:
        print("Splitting from ", G, "->", H)
        for i,positions in enumerate(wc.results[G]['compatible_supergroup_positions']):
            print("SPECIES=",wc.results[G]["species"][i])
            for position in positions:
                print(position)
        print('COMPATIBLE=',wc.results[G]['compatible_split'])
        print('\n')
