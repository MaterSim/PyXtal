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

    def compatibility_calculator(self,G,H,relations,atoms,positions):
        H=sym.Group(H)
        G=sym.Group(G)
        results={}

        wyc_position_list=[(str(x.multiplicity)+x.letter) for x in G]
        wyc_position_list.reverse()
        atom_basis=np.unique(atoms)
        atom_positions=[]
        for atom in atom_basis:
            atom_positions.append([positions[i] for i,x in enumerate(atoms) if x==atom])



        species_list=[]
        atom_number_list=[]
        positions_list=[]
        good_splittings_list=[]
        compatible_split=True
        for i,species in enumerate(atom_basis):
            position_basis=np.unique(atom_positions[i])
            position_counts=[atom_positions[i].count(x) for x in position_basis]
            possible_wyc_indexes=[]
            total_units=0
            for j,x in enumerate(position_basis):
                total_units+=int(x[:-1])*position_counts[j]
            for j,split in enumerate(relations):
                if np.all([x in position_basis for x in split]):
                    possible_wyc_indexes.append(j)
            possible_wyc_positions=[wyc_position_list[x] for x in possible_wyc_indexes]
            position_blocks=[np.array([relations[j].count(pos) for pos in position_basis]) for j in possible_wyc_indexes]

            position_block_units=[sum([int(x[:-1])*block[j] for j,x in enumerate(position_basis)]) for block in position_blocks]
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
    wc=wyckoff_compatibility(H,atoms=atoms,positions=positions)
    print('MATERIAL: Na-Sb3-F10', "Atomic Species:",wc.results[176]['species'], 'Group Number:',H)
    print("Atomic multiplicities of", wc.results[176]['species_number']," and with wyckoff positions ", wc.results[176]['species_positions'])




    for G in wc.supergroups:
        print("Splitting from ", G, "->", H)
        for i,positions in enumerate(wc.results[G]['compatible_supergroup_positions']):
            print("SPECIES=",wc.results[G]["species"][i])
            for position in positions:
                print(position)
        print('COMPATIBLE=',wc.results[G]['compatible_split'])
        print('\n')
