from pyxtal.molecule import find_id_from_smile

smiles = [
          'o1cc2CCc2c1',
          'C1=C(SC=C1O)C#N',
          'C1CC2CC(C1)C(=O)NC2=O',
          'Nc2cccc(=NS(=O)(=O)c1ccccc1)[nH]2',
          'CCC',
          'N1C(=O)NC(=O)C1',
          'c1(c(cc(c(c1)NC(=O)C)C)N(=O)=O)N(=O)=O',
          'N1CCC1',
          'C=CC=O',
          's1c(=S)n(c2ccccc12)C(=S)N(C)C',
          'CC1=NC(=NC=C1)N.CC1=CC=CC=C1C(=O)O',
          'C1=CC=C(C(=C1)[N+]#N)[O-]',
          'C1=C(C(=CC(=C1Cl)Cl)[N+](=O)[O-])[N+](=O)[O-]',
          'CC(=O)C(=[N+]=[N-])S(=O)(=O)C1=CC=C(C=C1)Cl',
          'c2cnc1[nH+]cccc1c2.C(=CC(=O)[O-])C(=O)O',
          'CC1=CC=C(C=C1)S(=O)(=O)C2=C(N=C(S2)C3=CC=C(C=C3)NC(=O)OCC4=CC=CC=C4)C',
          'C1=C(C=C(C(=C1O)O)O)C(=O)O',
          'C(#N)C1=C(SC2=NSC(=C2S1)C#N)C#N',
          'C1=CC=C(C(=C1)C(=O)O)NC2=CC=C(C=C2)CCC3=CC(=C(C=C3)Cl)Cl',
          'NC(=[NH2+])S/C=C/C(=O)O', #XXIV
          'C1=C(C=C(C=C1[N+](=O)[O-])[N+](=O)[O-])C(=O)O', #XXV
          'CC1=CC2=C(C=C1)N3CC4=C(C=CC(=C4)C)N(C2)C3',     #XXV
          'C1=CC=C2C(=C1)C=CC(=C2C3=C(C=CC4=CC=CC=C43)NC(=O)C5=CC=CC=C5Cl)NC(=O)',
         ]
for smile in smiles:
    print(smile)
    print(find_id_from_smile(smile))
