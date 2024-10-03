# python -m unittest pyxtal/test_all.py
import unittest
import numpy as np
from pyxtal.lego.SO3 import SO3
from pyxtal import pyxtal
from ase import Atoms
from ase.build import bulk, sort

def get_rotated_cluster(struc, angle=0, axis='x'):
    s_new = struc.copy()
    s_new.rotate(angle, axis)
    cell = 17.22*np.eye(3)
    p_struc = Atoms(s_new.symbols.numbers, positions=s_new.positions, cell=cell, pbc=True)
    return p_struc

def get_perturbed_cluster(struc, p0, p1, eps):
    s_new = struc.copy()
    pos = s_new.positions
    pos[p0, p1] += eps
    cell = 17.22*np.eye(3)
    p_struc = Atoms(s_new.symbols.numbers, positions=pos, cell=cell, pbc=True)
    return p_struc

def get_perturbed_xtal(struc, p0, p1, eps):
    pos = struc.positions.copy()
    pos[p0, p1] += eps
    p_struc = struc.copy()
    p_struc.set_positions(pos)
    return p_struc

def get_dPdR_xtal(xtal, nmax, lmax, rc, eps):
    p0 = SO3(nmax=nmax, lmax=lmax, rcut=rc).calculate(xtal, derivative=True)
    shp = p0['x'].shape
    array1 = p0['dxdr']

    for j in range(shp[0]):
        for k in range(3):
            struc = get_perturbed_xtal(xtal, j, k, eps)
            p1 = SO3(nmax=nmax, lmax=lmax, rcut=rc).calculate(struc)
            array2 = (p1['x'] - p0['x'])/eps
            #if np.linalg.norm(array2) > 1e-2: print(j, k, array2)
            if not np.allclose(array1[:, j, :, k], array2, atol=1e-4):
                print('xtal', j, k, '\n', array1[:, j, :, k], '\n', array2)
            assert(np.allclose(array1[:, j, :, k], array2, atol=1e-4))

# Descriptors Parameters
eps = 1e-8
rc = 2.80
nmax, lmax = 2, 2

# NaCl cluster
cluster = bulk('NaCl', crystalstructure='rocksalt', a=5.691694, cubic=True)
cluster = sort(cluster, tags=[0, 4, 1, 5, 2, 6, 3, 7])
cluster.set_pbc((0,0,0))
cluster = get_rotated_cluster(cluster, angle=0.1) # Must rotate

# Diamond
class TestCluster(unittest.TestCase):
    struc = get_rotated_cluster(cluster)
    p0 = SO3(nmax=nmax, lmax=lmax, rcut=rc).calculate(struc, derivative=True)
    struc = get_rotated_cluster(cluster, 10, 'x')
    p1 = SO3(nmax=nmax, lmax=lmax, rcut=rc).calculate(struc)

    def test_SO3_rotation_variance(self):
        array1 = self.p0['x']
        array2 = self.p1['x']
        assert(np.allclose(array1, array2))

    def test_dPdR_vs_numerical(self):
        shp = self.p0['x'].shape
        array1 = self.p0['dxdr']

        for j in range(shp[0]):
            for k in range(3):
                struc = get_perturbed_cluster(cluster, j, k, eps)
                p2 = SO3(nmax=nmax, lmax=lmax, rcut=rc).calculate(struc)
                array2 = (p2['x'] - self.p0['x'])/eps
                assert(np.allclose(array1[:,j,:,k], array2, atol=1e-3))

class TestXtal(unittest.TestCase):

    def test_dPdR_diamond(self):
        c = pyxtal()
        c.from_prototype('diamond')
        get_dPdR_xtal(c.to_ase(), nmax, lmax, rc, eps)

    def test_dPdR_graphite(self):
        c = pyxtal()
        c.from_prototype('graphite')
        get_dPdR_xtal(c.to_ase(), nmax, lmax, rc, eps)

    def test_dPdR_random(self):
        x = [ 7.952, 2.606, 0.592, 0.926, 0.608, 0.307]
        c = pyxtal()
        c.from_spg_wps_rep(179, ['6a', '6a', '6a', '6a'], x)
        get_dPdR_xtal(c.to_ase(), nmax, lmax, rc, eps)

    def test_dPdR_random_P(self):
        x = [ 7.952, 2.606, 0.592, 0.926, 0.608, 0.307]
        c = pyxtal()
        c.from_spg_wps_rep(179, ['6a', '6a', '6a', '6a'], x)
        atoms = c.to_ase()
        f = SO3(nmax=nmax, lmax=lmax, rcut=rc)
        p0 = f.compute_p(atoms)
        _, p1 = f.compute_dpdr(atoms)
        _, p2 = f.compute_dpdr_5d(atoms)
        assert(np.allclose(p0, p1, atol=1e-3))
        assert(np.allclose(p0, p2, atol=1e-3))

if __name__ == "__main__":
    unittest.main()
