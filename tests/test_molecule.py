# python -m unittest pyxtal/test_all.py
import importlib.util
import os
import unittest

import numpy as np
import pymatgen.analysis.structure_matcher as sm
from pymatgen.core.structure import Molecule
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from pyxtal import pyxtal
from pyxtal.lattice import Lattice
from pyxtal.molecule import Orientation, pyxtal_molecule
from pyxtal.symmetry import Group, Wyckoff_position


def resource_filename(package_name, resource_path):
    package_path = importlib.util.find_spec(package_name).submodule_search_locations[0]
    return os.path.join(package_path, resource_path)


cif_path = resource_filename("pyxtal", "database/cifs/")
l01 = Lattice.from_matrix([[4.08, 0, 0], [0, 9.13, 0], [0, 0, 5.50]])
l02 = Lattice.from_para(4.08, 9.13, 5.50, 90, 90, 90)
wp1 = Wyckoff_position.from_group_and_index(36, 0)
wp2 = Wyckoff_position.from_group_and_letter(36, "4a")

REPRODUCIBILITY_SEEDS = [0, 1, 7, 42, 123, 456]
REPRODUCIBILITY_REPEATS = 5
REPRODUCIBILITY_CASES = [
    {"group": 14, "species": ["aspirin"], "numIons": [4]},
    {"group": 19, "species": ["H2O"], "numIons": [4]},
    {"group": 36, "species": ["Benzene"], "numIons": [4]},
]


def fingerprint_molecular_crystal(crystal):
    """Return arrays that should match across repeated runs with the same seed."""
    lattice = np.array(crystal.lattice.get_para(), dtype=float)
    positions = np.array([site.position for site in crystal.mol_sites], dtype=float)
    orientations = np.array([site.orientation.matrix for site in crystal.mol_sites], dtype=float)
    return lattice, positions, orientations


def generate_molecular_crystal(seed, group, species, numIons):
    crystal = pyxtal(molecular=True)
    crystal.from_random(3, group, species, numIons, random_state=seed)
    assert crystal.valid
    return crystal


class TestOrientationRandomState(unittest.TestCase):
    def test_copy_spawns_independent_rng(self):
        base = Orientation(np.eye(3), degrees=2, random_state=np.random.default_rng(42))
        o1 = base.copy()
        o2 = base.copy()
        m1 = o1.change_orientation(flip=True)
        m2 = o2.change_orientation(flip=True)
        assert not np.allclose(m1, m2)

    def test_change_orientation_reproducible(self):
        def orient_matrix(seed):
            ori = Orientation(np.eye(3), degrees=2, random_state=np.random.default_rng(seed))
            return ori.copy().change_orientation(flip=True)

        assert np.allclose(orient_matrix(7), orient_matrix(7))
        assert not np.allclose(orient_matrix(7), orient_matrix(8))

    def test_set_axis_and_get_matrix_use_local_rng(self):
        ori = Orientation(np.eye(3), degrees=2, random_state=np.random.default_rng(0))
        ori.set_axis()
        assert ori.axis.shape == (3,)
        assert np.linalg.norm(ori.axis) > 0
        assert ori.get_matrix().shape == (3, 3)

    def test_seeded_molecular_generation_reproducible(self):
        lat1, pos1, _ = fingerprint_molecular_crystal(generate_molecular_crystal(123, 14, ["aspirin"], [4]))
        lat2, pos2, _ = fingerprint_molecular_crystal(generate_molecular_crystal(123, 14, ["aspirin"], [4]))
        lat3, pos3, _ = fingerprint_molecular_crystal(generate_molecular_crystal(456, 14, ["aspirin"], [4]))
        assert np.allclose(lat1, lat2, lat3)
        assert np.allclose(pos1, pos2, pos3)
        assert not (np.allclose(lat1, lat3) and np.allclose(pos1, pos3))

    def test_repeated_generation_with_random_states(self):
        """Run several seeds multiple times and require identical fingerprints."""
        fingerprints_by_seed = {}

        for case in REPRODUCIBILITY_CASES:
            for seed in REPRODUCIBILITY_SEEDS:
                runs = []
                for _ in range(REPRODUCIBILITY_REPEATS):
                    crystal = generate_molecular_crystal(
                        seed,
                        case["group"],
                        case["species"],
                        case["numIons"],
                    )
                    runs.append(fingerprint_molecular_crystal(crystal))

                reference = runs[0]
                for i, result in enumerate(runs[1:], start=2):
                    assert np.allclose(reference[0], result[0]), (
                        f"lattice mismatch for seed={seed}, case={case}, run={i}"
                    )
                    assert np.allclose(reference[1], result[1]), (
                        f"position mismatch for seed={seed}, case={case}, run={i}"
                    )
                    assert np.allclose(reference[2], result[2]), (
                        f"orientation mismatch for seed={seed}, case={case}, run={i}"
                    )

                key = (case["group"], tuple(case["species"]), seed)
                fingerprints_by_seed[key] = reference

        # Different seeds should not all collapse to the same structure.
        unique_lattices = {
            tuple(np.round(fp[0], 6)) for fp in fingerprints_by_seed.values()
        }
        assert len(unique_lattices) > 1

    def test_repeated_orientation_changes_with_random_states(self):
        """Orientation copy/rotate path should be reproducible across repeated runs."""
        seeds = REPRODUCIBILITY_SEEDS
        repeats = REPRODUCIBILITY_REPEATS

        for seed in seeds:
            matrices = []
            for _ in range(repeats):
                ori = Orientation(np.eye(3), degrees=2, random_state=np.random.default_rng(seed))
                matrices.append(ori.copy().change_orientation(flip=True))

            reference = matrices[0]
            for i, matrix in enumerate(matrices[1:], start=2):
                assert np.allclose(reference, matrix), (
                    f"orientation mismatch for seed={seed}, run={i}"
                )

            other_seed = seed + 1000
            other = Orientation(
                np.eye(3),
                degrees=2,
                random_state=np.random.default_rng(other_seed),
            ).copy().change_orientation(flip=True)
            if seed not in (0, other_seed):
                assert not np.allclose(reference, other)


class TestMolecule(unittest.TestCase):
    def test_get_orientations_in_wp(self):
        m = pyxtal_molecule("Benzene")
        g = Group(61)
        assert len(m.get_orientations_in_wp(g[0])) == 1
        assert len(m.get_orientations_in_wp(g[1])) == 1
        assert len(m.get_orientations_in_wp(g[2])) == 1


class TestMolecular(unittest.TestCase):
    def test_single_specie(self):
        struc = pyxtal(molecular=True)
        struc.from_random(3, 36, ["H2O"], sites=[["8b"]])
        struc.to_file("tmp-molecular.cif")
        os.remove("tmp-molecular.cif")
        assert struc.valid

        # test space group
        pmg_struc = struc.to_pymatgen()
        sga = SpacegroupAnalyzer(pmg_struc)
        # print(sga.get_space_group_symbol())
        assert sga.get_space_group_number() >= 36

        # test rotation
        struc.mol_sites[0].orientation.axis
        struc.mol_sites[0].rotate(ax_vector=[1, 0, 0], angle=90)
        pmg_struc = struc.to_pymatgen()
        sga = SpacegroupAnalyzer(pmg_struc)
        # pmg_struc.to("cif", "1.cif")
        assert sga.get_space_group_symbol() == "Cmc2_1"
        # print(pmg_struc.frac_coords[:3])

    def test_sites(self):
        struc = pyxtal(molecular=True)
        struc.from_random(3, 19, ["H2O"])
        pmg_struc = struc.to_pymatgen()
        sga = SpacegroupAnalyzer(pmg_struc)
        assert sga.get_space_group_symbol() == "P2_12_12_1"

        struc = pyxtal(molecular=True)
        struc.from_random(3, 36, ["H2O"], [8], sites=[["4a", "4a"]])
        pmg_struc = struc.to_pymatgen()
        sga = SpacegroupAnalyzer(pmg_struc)
        assert sga.get_space_group_symbol() == "Cmc2_1"

    def test_sites_xyz(self):
        struc = pyxtal(molecular=True)
        sites = [{"4e": [0.77, 0.57, 0.53]}]
        lat = Lattice.from_para(11.43, 6.49, 11.19, 90, 83.31, 90, ltype="monoclinic")
        struc.from_random(3, 14, ["aspirin"], [4], lattice=lat, sites=sites)
        assert struc.valid

    def test_read(self):
        # test reading structure from external
        struc = pyxtal(molecular=True)
        struc.from_seed(seed=cif_path + "aspirin.cif", molecules=["aspirin"])
        assert struc.lattice.ltype == "monoclinic"
        pmg_struc = struc.to_pymatgen()
        sga = SpacegroupAnalyzer(pmg_struc)
        assert sga.get_space_group_symbol() == "P2_1/c"
        C = struc.subgroup_once(eps=0, H=4)
        pmg_s2 = C.to_pymatgen()
        assert sm.StructureMatcher().fit(pmg_struc, pmg_s2)

    def test_big_molecule(self):
        # print("test_big_molecule")
        for mol in ["ROY", "aspirin"]:
            struc = pyxtal(molecular=True)
            struc.from_random(3, 19, [mol], factor=1.4)
            assert struc.valid
            pair = struc.check_short_distances()
            if len(pair) > 0:
                print("short distances were detected")
                print(mol)
                print(pair)
            assert len(pair) == 0

    def test_from_random_site(self):
        spg, wps, elements, numIons = 224, [["24j"]], ["C"], [24]
        for _i in range(10):
            c = pyxtal()
            c.from_random(3, spg, elements, numIons, sites=wps)
            pair = c.check_short_distances(0.5)
            assert len(pair) == 0

    def test_c60(self):
        struc = pyxtal(molecular=True)
        struc.from_random(3, 36, ["C60"], [4], 1.0)
        assert struc.valid

    def test_mutiple_species(self):
        Li = Molecule(["Li"], [[0.0, 0.0, 0.0]])
        coords = [
            [0.000000, 0.000000, 0.000000],
            [1.200000, 1.200000, -1.200000],
            [1.200000, -1.200000, 1.200000],
            [-1.200000, 1.200000, 1.200000],
            [-1.200000, -1.200000, -1.200000],
        ]
        ps4 = Molecule(["P", "S", "S", "S", "S"], coords)

        for _i in range(3):
            struc = pyxtal(molecular=True)
            struc.from_random(3, 10, [Li, ps4], [6, 2], 1.2, conventional=False)
            if struc.valid:
                assert len(struc.to_pymatgen()) == 16

    def test_distance(self):
        Ag_xyz = """1
        AgC2N2H
        Ag         4.30800        8.26300       -0.2200
        """

        C2N2H7_xyz = """12
        AgC2N2H
        H          5.95800        5.80600       -0.9530
        N          5.24100        6.16800       -1.1210
        N          2.23200        6.99000       -0.6820
        C          4.02900        5.47000       -1.5870
        C          2.78500        5.61100       -0.6610
        H          3.69300        5.63500       -2.5830
        H          4.17400        4.42800       -1.6990
        H          3.12700        5.50500        0.4260
        H          2.05000        5.01200       -0.9300
        H          1.96000        7.20500       -1.3860
        H          1.59400        6.99200       -0.0710
        """
        with open("Ag.xyz", "w") as f:
            f.write(Ag_xyz)
        with open("C2N2H7.xyz", "w") as f:
            f.write(C2N2H7_xyz)

        for _i in range(10):
            c = pyxtal(molecular=True)
            c.from_random(3, 9, ["Ag.xyz", "C2N2H7.xyz"], [12, 12])
            short_bonds = c.check_short_distances(r=1.1)
            assert len(short_bonds) == 0
        os.remove("Ag.xyz")
        os.remove("C2N2H7.xyz")

    def test_molecular_2d(self):
        # print("test_molecular_2d")
        struc = pyxtal(molecular=True)
        struc.from_random(2, 20, ["H2O"])
        assert struc.valid

    def test_molecular_1d(self):
        struc = pyxtal(molecular=True)
        struc.from_random(1, 20, ["H2O"])
        assert struc.valid

    def test_preassigned_sites(self):
        sites = [["4a", "4a"]]
        struc = pyxtal(molecular=True)
        struc.from_random(3, 36, ["H2O"], [8], sites=sites)
        assert struc.valid

    def test_special_sites(self):
        struc = pyxtal(molecular=True)
        struc.from_random(3, 61, ["Benzene"], [4])
        assert struc.valid

if __name__ == "__main__":
    unittest.main()
