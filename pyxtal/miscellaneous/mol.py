from rdkit import Chem
from rdkit.Chem import AllChem
import pymatgen as mg
from glob import glob
data = {}
data["HAHCOI"] = "s1c2ccc3scc4ccc(c1)c2c34"
data["JAPWIH"] = "s1ccc2cc3sc4cc5ccsc5cc4c3cc12"
data["WEXBOS"] = "s1c(c2ccccc2)c(c2ccccc2)c2c1c(c(s2)c1ccccc1)c1ccccc1"
data["LAGNAL"] = "s1c(/C=N/[C@H](C)c2ccc(F)cc2)ccc1/C=N/[C@H](C)c1ccc(F)cc1"
data["YICMOP"] = "s1cccc1c1c(F)c(OC)c(c2sccc2)c(F)c1OC"
data["MERQIM"] = "s1c2c(c3c1SCCC3)cc1sc3SCCCc3c1c2"
data["LUFHAW"] = "CC1=CC2=C(S1)C3=CC4=C(C=C3C=C2)C5=C(C=C4)C=C(S5)C"
#smi, cif = "s1c2ccc3scc4ccc(c1)c2c34", "HAHCOI.cif"
#smi, cif = "s1ccc2cc3sc4cc5ccsc5cc4c3cc12", "JAPWIH.cif"
#smi, cif = "s1c(c2ccccc2)c(c2ccccc2)c2c1c(c(s2)c1ccccc1)c1ccccc1", "WEXBOS.cif"
#smi, cif = "s1c(/C=N/[C@H](C)c2ccc(F)cc2)ccc1/C=N/[C@H](C)c1ccc(F)cc1","LAGNAL.cif"
#smi, cif = "s1cccc1c1c(F)c(OC)c(c2sccc2)c(F)c1OC", "YICMOP.cif"
#smi, cif = "s1c2c(c3c1SCCC3)cc1sc3SCCCc3c1c2", "MERQIM.cif"
#smi, cif = "CC1=CC2=C(S1)C3=CC4=C(C=C3C=C2)C5=C(C=C4)C=C(S5)C", "LUFHAW.cif"

for file in glob("*.cif"):
    name = file[:-4]
    if name in data.keys():
        smi = data[name]
        m = Chem.MolFromSmiles(smi)
        m2 = Chem.AddHs(m)
        AllChem.EmbedMolecule(m2)
        cids = AllChem.EmbedMultipleConfs(m2, numConfs=1)
        xyz = Chem.rdmolfiles.MolToXYZBlock(m2, 0)
        mol = mg.Molecule.from_str(xyz, fmt="xyz")
        mol.to(filename=name+".xyz")
