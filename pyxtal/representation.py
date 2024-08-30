import numpy as np

from pyxtal.lattice import Lattice
from pyxtal.molecule import find_rotor_from_smile
from pyxtal.symmetry import Group
from pyxtal.wyckoff_site import atom_site, mol_site


class representation_atom:
    """
    A class to handle the 1D representation of atomic crystal
    Works for Zprime > 1

    Args:
        x: a list of [cell, site_1, site_2, ...]
    """

    def __init__(self, x):
        self.x = x

    def __str__(self):
        return self.to_string()

    @classmethod
    def from_pyxtal(cls, struc, standard=True):
        """
        Initialize 1D rep. from the pyxtal object

        Args:
            struc: pyxtal object
        """
        if standard and not struc.standard_setting:
            pmg = struc.to_pymatgen()
            struc.from_seed(pmg, standard=True)
        symmetry = [struc.atom_sites[0].wp.hall_number]
        lat = struc.lattice.encode()
        vector = [symmetry + lat]
        vector.extend([site.encode() for site in struc.atom_sites])
        x = vector
        return cls(x)

    def to_standard_setting(self):
        xtal = self.to_pyxtal()
        self.x = representation.from_pyxtal(xtal, standard=True).x

    def to_pyxtal(self):
        """
        Export the pyxtal structure

        Args:
            smiles: list of smiles
            compoisition: list of composition
        """
        from pyxtal import pyxtal

        # symmetry
        v = self.x[0]
        struc = pyxtal()
        struc.group, _number = Group(v[0], use_hall=True), v[0]

        # lattice
        ltype = struc.group.lattice_type
        if ltype == "triclinic":
            a, b, c, alpha, beta, gamma = v[1], v[2], v[3], v[4], v[5], v[6]
        elif ltype == "monoclinic":
            a, b, c, alpha, beta, gamma = v[1], v[2], v[3], 90, v[4], 90
        elif ltype == "orthorhombic":
            a, b, c, alpha, beta, gamma = v[1], v[2], v[3], 90, 90, 90
        elif ltype == "tetragonal":
            a, b, c, alpha, beta, gamma = v[1], v[1], v[2], 90, 90, 90
        elif ltype == "hexagonal":
            a, b, c, alpha, beta, gamma = v[1], v[1], v[2], 90, 90, 120
        else:
            a, b, c, alpha, beta, gamma = v[1], v[1], v[1], 90, 90, 90
        try:
            struc.lattice = Lattice.from_para(a, b, c, alpha, beta, gamma, ltype=ltype)
        except:
            print(a, b, c, alpha, beta, gamma, ltype)
            raise ValueError("Problem in Lattice")

        # sites
        struc.numIons = []
        struc.atom_sites = []
        species = []

        for _x in self.x[1:]:
            dicts = {}
            specie, index, pos = _x[0], _x[1], _x[2:]
            dicts["specie"] = specie
            dicts["index"] = index
            dicts["dim"] = 3
            dicts["PBC"] = [1, 1, 1]
            dicts["hn"] = struc.group.hall_number
            wp = struc.group[index]
            dicts["position"] = wp.get_position_from_free_xyzs(pos)
            site = atom_site.load_dict(dicts)
            struc.atom_sites.append(site)

            if specie not in species:
                species.append(specie)
                struc.numIons.append(wp.multiplicity)
            else:
                for i, _specie in enumerate(species):
                    if _specie == specie:
                        struc.numIons[i] += site.wp.multiplicity

        struc.species = species
        struc._get_formula()
        struc.source = "1D rep."
        struc.valid = True
        struc.standard_setting = site.wp.is_standard_setting()

        return struc

    def to_array(self):
        """
        Export only varibles to a 1d numpy array
        """
        cells, xyzs = self.x[0][1:], self.x[1:]
        x = cells
        for xyz in xyzs:
            x = np.hstack((x, xyz[2:]))
        return x

    def to_string(self, time=None, eng=None, tag=None):
        """
        Export string representation

        Args:
            time: float
            eng: float
            tag: string
        """
        x = self.x
        strs = f"{int(x[0][0]):3d} "

        # data for cell
        if x[0][0] <= 348:
            num = 4
        elif x[0][0] <= 488:
            num = 3
        else:  # cubic
            num = 2

        for c in x[0][1:num]:
            strs += f"{c:5.2f} "
        for c in x[0][num:]:
            strs += f"{c:5.1f} "

        # data for atoms
        strs += f"{len(x) - 1:d} "  # Number of sites
        for i in range(1, len(x)):
            strs += f"{x[i][0]:s} "
            strs += f"{x[i][1]:d} "
            for v in x[i][2:]:
                strs += f"{v:6.4f} "

        if time is not None:
            strs += f"{time:5.2f}"

        if eng is not None:
            strs += f"{eng:11.3f}"

        if tag is not None:
            strs += f" {tag:s}"

        return strs


class representation:
    """
    A class to handle the 1D representation of molecular crystal
    Works for Zprime > 1

    Args:
        x: a list of [cell, site_1, site_2, ...]
        smiles: a list of [smile_1, smile_2, ...]
    """

    def __init__(self, x, smiles=None):
        if smiles is not None:
            self.smiles = []
            for _i, smile in enumerate(smiles):
                if smile.endswith(".smi"):
                    smile = smile[:-4]
                self.smiles.append(smile)
        else:
            self.smiles = None
        self.x = x

    def __str__(self):
        return self.to_string()

    @classmethod
    def from_pyxtal(cls, struc, standard=False):
        """
        Initialize 1D rep. from the pyxtal object

        Args:
            struc: pyxtal object
        """
        if standard and not struc.standard_setting:
            # struc.optimize_lattice(standard=True)
            pmg = struc.to_pymatgen()
            struc.from_seed(pmg, molecules=struc.molecules, standard=True)
        symmetry = [struc.mol_sites[0].wp.hall_number]
        lat = struc.lattice.encode()
        vector = [symmetry + lat]
        smiles = []
        for site in struc.mol_sites:
            vector.append(site.encode())
            smiles.append(site.molecule.smile)
        x = vector
        if smiles[0] is None: smiles = None
        return cls(x, smiles)

    @classmethod
    def from_string(cls, inputs, smiles, composition=None):
        """
        Initialize 1D rep. from the string

        Args:
            inputs: input string
            smiles: list of smiles
        """
        # parse the cell
        if composition is None:
            composition = [1] * len(smiles)

        inputs = [float(tmp) for tmp in inputs.split()]
        hn = int(inputs[0])
        if hn <= 2:
            n_cell = 8
        elif hn <= 107:
            n_cell = 6
        elif hn <= 348:
            n_cell = 5
        elif hn <= 488:
            n_cell = 4
        else:
            n_cell = 3  # cubic
        cell = [hn] + inputs[1 : n_cell - 1]

        x = [cell]
        n_site = int(inputs[n_cell - 1])
        if n_site != sum(composition):
            msg = f"Composition is inconsistent: {sum(composition):d}/{n_site:d}\n"
            msg += str(inputs)
            raise ValueError(msg)
        # n_cell += 1

        for i, smile in enumerate(smiles):
            if smile.endswith(".smi"):
                smile = smile[:-4]
            for _c in range(composition[i]):
                if smile in ["Cl-"]:
                    n_mol = 4
                else:
                    n_torsion = len(find_rotor_from_smile(smile))
                    n_mol = 8 + n_torsion  # (wp_id, x, y, z, ori_x, ori_y, ori_z, inv) + torsion
                # inversion
                # print(n_mol, n_cell, len(inputs))
                inputs[n_cell] = int(inputs[n_cell])
                inputs[n_cell + n_mol - 1] = int(inputs[n_cell + n_mol - 1])
                x.append(inputs[n_cell : n_cell + n_mol])  # ; print('string', x[-1])
                n_cell += n_mol
        assert n_cell == len(inputs)
        return cls(x, smiles)

    def to_standard_setting(self):
        xtal = self.to_pyxtal()
        xtal.optimize_lattice(standard=True)
        rep0 = representation.from_pyxtal(xtal)
        self.x = rep0.x

    def to_pyxtal(self, smiles=None, composition=None):
        """
        Export the pyxtal structure

        Args:
            smiles: list of smiles
            compoisition: list of composition
        """
        from pyxtal import pyxtal

        if smiles is None:
            smiles = self.smiles

        if composition is None:
            composition = [1] * len(smiles)

        if sum(composition) + 1 != len(self.x):
            msg = "Composition is inconsistent:\n"
            msg += str(composition) + "\n"
            msg += self.to_string()
            raise ValueError(msg)

        # symmetry
        v = self.x[0]
        struc = pyxtal(molecular=True)
        struc.group, _number = Group(v[0], use_hall=True), v[0]

        # lattice
        ltype = struc.group.lattice_type
        if ltype == "triclinic":
            a, b, c, alpha, beta, gamma = v[1], v[2], v[3], v[4], v[5], v[6]
        elif ltype == "monoclinic":
            a, b, c, alpha, beta, gamma = v[1], v[2], v[3], 90, v[4], 90
        elif ltype == "orthorhombic":
            a, b, c, alpha, beta, gamma = v[1], v[2], v[3], 90, 90, 90
        elif ltype == "tetragonal":
            a, b, c, alpha, beta, gamma = v[1], v[1], v[2], 90, 90, 90
        elif ltype == "hexagonal":
            a, b, c, alpha, beta, gamma = v[1], v[1], v[2], 90, 90, 120
        else:
            a, b, c, alpha, beta, gamma = v[1], v[1], v[1], 90, 90, 90
        try:
            struc.lattice = Lattice.from_para(a, b, c, alpha, beta, gamma, ltype=ltype)
        except:
            print(a, b, c, alpha, beta, gamma, ltype)
            raise ValueError("Problem in Lattice")

        # sites
        struc.numMols = [0] * len(smiles)
        struc.molecules = []
        struc.mol_sites = []

        count = 1
        for i, comp in enumerate(composition):
            smile = smiles[i]
            if smile.endswith(".smi"):
                smile = smile[:-4]
            for _j in range(comp):
                v = self.x[count]
                dicts = {}
                dicts["smile"] = smile
                dicts["type"] = i
                dicts["dim"] = 3
                dicts["PBC"] = [1, 1, 1]
                # dicts['number'] = number
                dicts["hn"] = struc.group.hall_number
                dicts["index"] = v[0]
                dicts["lattice"] = struc.lattice.matrix
                dicts["lattice_type"] = ltype
                dicts["center"] = v[1:4]
                if smile not in ["Cl-"]:
                    dicts["orientation"] = np.array(v[4:7])
                    dicts["rotor"] = v[7:-1]  # ; print('ro', dicts['rotor'])
                    dicts["reflect"] = int(v[-1])
                site = mol_site.from_1D_dicts(dicts)

                bypass = False
                for mol_id, molecule in enumerate(struc.molecules):
                    if str(site.molecule) == str(molecule):
                        site.type = mol_id
                        struc.numMols[mol_id] += site.wp.multiplicity
                        bypass = True
                        break
                if not bypass:
                    struc.molecules.append(site.molecule)
                    site.type = len(struc.molecules) - 1
                    struc.numMols[site.type] += site.wp.multiplicity

                # site.type = i
                struc.mol_sites.append(site)
                # move to next rep
                count += 1

        struc._get_formula()
        struc.source = "1D rep."
        struc.valid = True
        struc.standard_setting = site.wp.is_standard_setting()

        return struc

    def to_string(self, time=None, eng=None, tag=None):
        """
        Export string representation

        Args:
            time: float
            eng: float
            tag: string
        """
        x = self.x
        strs = f"{int(x[0][0]):3d} "

        # data for cell
        if x[0][0] <= 348:
            num = 4
        elif x[0][0] <= 488:
            num = 3
        else:  # cubic
            num = 2

        for c in x[0][1:num]:
            strs += f"{c:5.2f} "
        for c in x[0][num:]:
            strs += f"{c:5.1f} "

        # data for molecule
        strs += f"{len(x) - 1:d} "  # ; print(x[1])
        for i in range(1, len(x)):
            strs += f"{x[i][0]:d} "
            for v in x[i][1:4]:
                strs += f"{v:4.2f} "
            for v in x[i][4:-1]:
                strs += f"{v:6.1f} "
            strs += f"{int(x[i][-1]):d} "

        if time is not None:
            strs += f"{time:5.2f}"

        if eng is not None:
            strs += f"{eng:11.3f}"

        if tag is not None:
            strs += f" {tag:s}"

        return strs

    def same_smiles(self, smiles):
        if len(self.smiles) == smiles:
            return all(s2 == s2 for s1, s2 in zip(self.smiles, smiles))
        else:
            return False

    def get_dist(self, rep):
        """
        get distance with the other rep1
        Now only supports Z'=1
        """
        from pyxtal.symmetry import Wyckoff_position as WP

        if self.same_smiles(rep.smiles):
            msg = "different smiles"
            print(msg)
            return None
        elif len(self.x) != len(rep.x):
            msg = "different number of sites"
            print(msg)
            return None
        elif self.x[0][0] != rep.x[0][0]:
            msg = "different space group numbers"
            print(msg)
            return None
        else:
            diffs = []
            wp = WP.from_group_and_index(self.x[0][0], 0, use_hall=True)
            for i in range(len(self.x)):
                np.zeros(len(self.x[i]))
                tmp1 = np.array(self.x[i])
                tmp2 = np.array(rep.x[i])
                # cell difference
                if i == 0:
                    diff_cell = tmp2 - tmp1
                    diffs.extend(diff_cell)
                # site difference
                else:
                    # symmmetry variation
                    xyzs = wp.apply_ops(tmp2[:3])
                    diff_xyzs = xyzs - tmp1[:3]
                    diff_xyzs -= np.rint(diff_xyzs)
                    id = np.argmin(np.linalg.norm(diff_xyzs, axis=1))
                    diff_xyz = diff_xyzs[id]
                    diff_ori = tmp2[3:6] - tmp1[3:6]
                    diff_ori /= [360.0, 180.0, 360.0]
                    diff_ori -= np.rint(diff_ori)
                    diff_ori *= [360.0, 180.0, 360.0]
                    diff_tor = tmp2[6:] - tmp1[6:]
                    diff_tor /= 360.0
                    diff_tor -= np.rint(diff_tor)
                    diff_tor *= 360.0
                    diffs.extend(diff_xyz)
                    diffs.extend(diff_ori)
                    diffs.extend(diff_tor)
            return np.array(diffs)


if __name__ == "__main__":
    # aspirin
    smiles = ["CC(=O)OC1=CC=CC=C1C(=O)O"]
    x = [
        [81, 11.43, 6.49, 11.19, 83.31],
        [0, 0.77, 0.57, 0.53, 48.55, 24.31, 145.94, -77.85, -4.40, 170.86, False],
    ]
    # rep0 = representation(x, smiles)
    # print(rep0.to_string())
    rep1 = representation(x, smiles)
    xtal = rep1.to_pyxtal()
    print(xtal)
    rep2 = representation.from_pyxtal(xtal)
    print(rep2.to_pyxtal())
    print(rep2.to_string())

    print("Test read from string")
    string = "82 11.43  6.49 11.19 83.31 1  0 0.77  0.57  0.53 48.55 24.31 145.9 -77.85 -4.40 170.9 0"
    rep3 = representation.from_string(string, smiles)
    print(rep3.to_string())
    print(rep3.to_pyxtal())
    # x = rep3.to_pyxtal(); x.optimize_lattice(standard=True); print(x)
    rep3.to_standard_setting()
    print(rep3.to_pyxtal())
    print(rep3.to_string())

    print("Test other cases")
    string1 = "81 14.08  6.36 25.31  83.9 1 0 0.83 0.40 0.63  136.6  -21.6 -151.1 -101.1 -131.2  154.7 -176.4 -147.8  178.2 -179.1  -53.3 0"  # noqa: E501
    string2 = "81 14.08  6.36 25.31  83.9 1 0 0.03 0.84 0.89  149.1   -8.0  -37.8  -39.9 -104.2  176.2 -179.6  137.8 -178.5 -173.3 -103.6 0"  # noqa: E501
    smiles = ["CC1=CC=C(C=C1)S(=O)(=O)C2=C(N=C(S2)C3=CC=C(C=C3)NC(=O)OCC4=CC=CC=C4)C"]
    rep4 = representation.from_string(string1, smiles)
    rep5 = representation.from_string(string2, smiles)
    print(string1)
    print(string2)
    print(rep4.get_dist(rep5))

    from pyxtal import pyxtal

    xtal = pyxtal()
    xtal.from_seed("pyxtal/database/cifs/Fd3.cif")
    xtal.from_seed("pyxtal/database/cifs/NaSb3F10.cif")
    rep = representation_atom.from_pyxtal(xtal)
    print(rep)
    print(xtal)
    print(rep.to_pyxtal())
    # strings = [
    # "83 14.08  6.36 25.31  83.9 1 0.72 0.40 0.27  131.6  -17.0 -120.0  -83.8 -134.1 -174.5 -175.7 -168.8  173.9  178.0 -157.4 0",  # noqa: E501
    # "81 14.08  6.36 25.31  83.9 1 0.59 0.81 0.39 -117.8  -50.1  -95.3  -25.8  -80.6  164.7  155.9 -124.9 -159.2  178.6 -154.7 0",  # noqa: E501
    # "81 14.08  6.36 25.31  83.9 1 0.75 0.09 0.01  133.8  -19.5  -55.1  -86.7  -91.7 -175.0 -170.4 -176.8  173.3 -164.8  -58.4 0",  # noqa: E501
    # "81 14.08  6.36 25.31  83.9 1 0.72 0.44 0.01  135.2   27.5   97.2 -101.1 -105.1  -29.7 -169.7  -50.1  172.2 -173.1  131.6 0",  # noqa: E501
    # "82 14.00  6.34 25.26  83.6 1 0.21 0.08 0.54  146.0  -12.0   50.2  108.0  112.3 -166.3 -158.7  -35.5  172.3 -168.7  133.0 0",  # noqa: E501
    # "81 14.08  6.36 25.31  83.9 1 0.05 0.30 0.89  -68.2   41.2  148.8  -66.9  -85.0 -167.4  172.3 -166.2 -178.3  166.4  -45.9 0",  # noqa: E501
    # ]

    # import pymatgen.analysis.structure_matcher as sm
    # matcher = sm.StructureMatcher(ltol=0.3, stol=0.3, angle_tol=10)
    # for i, string in enumerate(strings):
    #    print(str(i) + '  ' +string)
    #    rep4 = representation.from_string(string, smiles)
    #    pmg1 = rep4.to_pyxtal().to_pymatgen(); pmg1.remove_species('H')
    #    rep4.to_standard_setting()
    #    pmg2 = rep4.to_pyxtal().to_pymatgen(); pmg2.remove_species('H')
    #    print(i, rep4.to_string(), matcher.fit(pmg1, pmg2))
