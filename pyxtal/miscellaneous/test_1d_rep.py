from pyxtal import pyxtal
import pymatgen.analysis.structure_matcher as sm
from pyxtal.representation import representation
from pyxtal.db import database
from warnings import simplefilter
simplefilter(action='ignore', category=DeprecationWarning)

from ccdc.crystal import PackingSimilarity
from ccdc.io import CrystalReader
engine = PackingSimilarity()

import numpy as np
np.seterr(all="ignore")

matcher = sm.StructureMatcher(ltol=0.3, stol=0.3, angle_tol=10)
codes = ['XULDUD', 'ABACOW', 'DAZJOG', 'BINREY', 'CUVLIO', 'FUPKIK', 'FURWOG', 'HUJFAV', 
         'IJEHOU', 'IFAQUE', 'JUMTAP', 'LELVIK', 'LOPQAL']

struc = pyxtal(molecular=True)
for i, code in enumerate(codes):
        struc.from_CSD(code)
        smiles = struc.tag['smiles'].split('.')
        rep = struc.get_1D_representation()
        pmg1 = struc.to_pymatgen()
        struc1 = representation.from_string(rep.to_string(), smiles).to_pyxtal()
        pmg2 = struc1.to_pymatgen()
        pmg1.remove_species('H')
        pmg2.remove_species('H')
        print(i, code, struc.group.symbol)
        if matcher.fit(pmg1, pmg2):
            #print(code, smiles, matcher.fit(pmg1, pmg2))
            print(rep.to_string())
            rep.to_standard_setting()
            print(rep.to_string())
        else:
            struc.to_file('1.cif')
            struc1.to_file('2.cif')
            wp = struc.mol_sites[0]
            xyz, _ = wp._get_coords_and_species(absolute=True, first=True)
            print(code)
            print(wp.molecule.torsionlist)
            wp.molecule.get_rmsd(xyz, debug=True)

            cryst0 = CrystalReader('1.cif')[0]
            cryst1 = CrystalReader('2.cif')[0]
            h = engine.compare(cryst0, cryst1) 
            #if h is None or h.nmatched_molecules<=10:
            #    print(h.nmatched_molecules)
            #    import sys; sys.exit()

smiles = ['CC1=CC=C(C=C1)S(=O)(=O)C2=C(N=C(S2)C3=CC=C(C=C3)NC(=O)OCC4=CC=CC=C4)C']
strings = [
         "83 14.08  6.36 25.31  83.9 1 0.72 0.40 0.27  131.6  -17.0 -120.0  -83.8 -134.1 -174.5 -175.7 -168.8  173.9  178.0 -157.4 0",
         "81 14.08  6.36 25.31  83.9 1 0.83 0.40 0.63  136.6  -21.6 -151.1 -101.1 -131.2  154.7 -176.4 -147.8  178.2 -179.1  -53.3 0",
         "81 14.08  6.36 25.31  83.9 1 0.59 0.81 0.39 -117.8  -50.1  -95.3  -25.8  -80.6  164.7  155.9 -124.9 -159.2  178.6 -154.7 0",
         "81 14.08  6.36 25.31  83.9 1 0.03 0.84 0.89  149.1   -8.0  -37.8  -39.9 -104.2  176.2 -179.6  137.8 -178.5 -173.3 -103.6 0",
         "81 14.08  6.36 25.31  83.9 1 0.75 0.09 0.01  133.8  -19.5  -55.1  -86.7  -91.7 -175.0 -170.4 -176.8  173.3 -164.8  -58.4 0",
         "81 14.08  6.36 25.31  83.9 1 0.72 0.44 0.01  135.2   27.5   97.2 -101.1 -105.1  -29.7 -169.7  -50.1  172.2 -173.1  131.6 0",
         "82 14.00  6.34 25.26  83.6 1 0.21 0.08 0.54  146.0  -12.0   50.2  108.0  112.3 -166.3 -158.7  -35.5  172.3 -168.7  133.0 0",
         #"83 14.08  6.36 25.31  83.9 1 0.25 0.91 0.06  136.4   20.0   52.2  -55.5  -79.7 -176.4 -160.9  109.2  168.9 -179.7  -25.7 0",
         #"82 14.08  6.36 25.31  83.9 1 0.55 0.04 0.64 -130.1   15.3  -10.7   74.9   64.0  177.4  180.0    8.7  177.8  174.6 -174.2 0",
         #"82 14.08  6.36 25.31  83.9 1 0.99 0.02 0.62   38.4   11.6  -16.2  -86.3   71.3  -19.5  164.9 -176.4 -168.2  174.4   88.1 0",
         #"83 14.08  6.36 25.31  83.9 1 0.67 0.40 0.87  136.6  -21.6   28.9 -101.1 -131.2  154.7 -176.4 -147.8  178.2 -179.1  -53.3 0",
         #"82 14.08  6.36 25.31  83.9 1 0.18 0.53 0.12  -84.9   36.3   89.1   15.6   54.8  140.6 -179.4 -138.3  178.5  169.7    4.7 0",
         #"82 14.08  6.36 25.31  83.9 1 0.47 0.75 0.11  -44.4  -26.0 -138.6  -50.5  -64.5 -174.1  161.8   20.8 -167.9  171.7   68.2 0",
         #"82 14.08  6.36 25.31  83.9 1 0.75 0.43 0.55  -44.9   16.9   58.5  -85.6  -87.7 -175.5 -170.6 -174.4  175.3 -170.0  -82.7 0",
         #"83 14.08  6.36 25.31  83.9 1 0.52 0.50 0.62  -39.5   -0.4  -32.5  -59.9  -76.0   11.3  169.5  166.5 -148.6  -87.3   82.6 0",
         #"83 14.08  6.36 25.31  83.9 1 0.42 0.98 0.67  141.3  -16.3 -123.8   62.9  -73.5  -39.4  178.8 -151.7 -175.1 -171.6   35.4 0",
         #"82 14.08  6.36 25.31  83.9 1 0.96 0.16 0.37  -41.7   45.2   -6.2  -93.1  -81.2    8.3 -146.8  147.2  147.2  -76.0   84.9 0",
         #"83 14.08  6.36 25.31  83.9 1 0.92 0.68 0.13  110.3  -46.9 -152.1  -80.7  -79.4 -172.6 -173.1 -140.8  169.0  165.6  155.7 0",
         ]
for i, string in enumerate(strings):
    print(str(i) + '  ' +string)
    rep4 = representation.from_string(string, smiles)
    pmg1 = rep4.to_pyxtal().to_pymatgen(); pmg1.remove_species('H')
    rep4.to_standard_setting()
    pmg2 = rep4.to_pyxtal().to_pymatgen(); pmg2.remove_species('H')
    print(i, rep4.to_string(), matcher.fit(pmg1, pmg2))
