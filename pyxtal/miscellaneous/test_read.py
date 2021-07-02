from pyxtal.molecular_crystal import molecular_crystal
import numpy as np
from pymatgen.core.structure import Structure
from rdkit import Chem
from rdkit.Chem import AllChem
from pyxtal.symmetry import Group
import pymatgen.analysis.structure_matcher as sm
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pkg_resources import resource_filename

lists = [
("s1c2ccc3scc4ccc(c1)c2c34", "HAHCOI.cif"), #P21 -> P21
("s1ccc2cc3sc4cc5ccsc5cc4c3cc12", "JAPWIH.cif"), #Pmn21 -> P21
("s1c(c2ccccc2)c(c2ccccc2)c2c1c(c(s2)c1ccccc1)c1ccccc1", "WEXBOS.cif"), #P21/c->Pc
("s1c(/C=N/[C@H](C)c2ccc(F)cc2)ccc1/C=N/[C@H](C)c1ccc(F)cc1","LAGNAL.cif"), #
("CC1=CC2=C(S1)C3=CC4=C(C=C3C=C2)C5=C(C=C4)C=C(S5)C", "LUFHAW.cif"), #P21/n->Pn
("s1c2c(c3c1SCCC3)cc1sc3SCCCc3c1c2", "MERQIM.cif"), #P21/a - > Pc
("s1cccc1c1c(F)c(OC)c(c2sccc2)c(F)c1OC", "YICMOP.cif"), #P21/n -> Pn
("s1c(ccc1C1=N[C@H](CO1)CC)C1=N[C@H](CO1)CC", "DAFJOK.cif"), #C2 -> P21
("CN1C2=C(N=C1C3=CC=C(S3)C=O)N(C(=O)N(C2=O)C)C", "LUTXII.cif"), #P21/m -> P21
]

for mol in lists:
    (smi, cif) = mol
    m = Chem.MolFromSmiles(smi)
    cif = resource_filename("pyxtal", "database/cifs/") + cif
    m2 = Chem.AddHs(m)
    m3 = Chem.AddHs(m2)
    AllChem.EmbedMultipleConfs(m3)
    Chem.rdmolfiles.MolToXYZFile(m3, '1.xyz') 

    # Original Structure
    struc = molecular_crystal(2, ['1.xyz'], [2], seed=cif)    
    #print(struc.mol_sites[0])
    s = struc.copy()
    site = s.mol_sites[0].make_gen_wyckoff_site()
    print("Site distance", site.check_distances())
    s.mol_sites = [site]
    s.lattice = site.lattice
    s.group = Group(site.wp.number)
    #print([op.as_xyz_string() for op in s.mol_sites[0].wp.ops])
    pmg1 = s.to_pymatgen()
    pmg2 = struc.to_pymatgen()
    print(SpacegroupAnalyzer(pmg1).get_space_group_symbol())
    print(cif, s.group.symbol, struc.group.symbol, "Match: ", sm.StructureMatcher().fit(pmg1, pmg2))
    #print(struc.to_file())
    s3 = Structure.from_str(s.to_file(), fmt="cif")
    s0 = Structure.from_file(cif)
    print("read-source Match: ", sm.StructureMatcher().fit(s3, s0))
    print(np.sort(s3.lattice.abc))
    print(np.sort(s0.lattice.abc))
    #print(pmg1.to(fmt='cif'))
    #print(pmg2.to(fmt='cif'))
    #print(s0.to(fmt='cif'))



#from pyxtal.symmetry import Wyckoff_position
#for cif in ["HAHCOI.cif", "JAPWIH.cif", "LAGNAL.cif", 
#            "MERQIM.cif", "WEXBOS.cif", "YICMOP.cif",
#            "LUFHAW.cif",
#           ]:
#    pmg_struc = Structure.from_file(cif)
#    sga = SpacegroupAnalyzer(pmg_struc)
#    ops = sga.get_symmetry_operations()
#    strings = [op.as_xyz_string() for op in ops]
#    print(sga.get_space_group_symbol(), strings)
#    wyc, perm = Wyckoff_position.from_symops(ops, sga.get_space_group_number())
#    #print(perm)
