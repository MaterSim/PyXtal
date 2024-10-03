# python -m unittest pyxtal/test_all.py
import unittest
import numpy as np
from pyxtal.lego.SO3 import SO3
from pyxtal import pyxtal
from ase import Atoms
from ase.build import bulk, sort

def calculate_S(atoms, P_ref):
    P = calculator.compute_p(atoms)
    S = np.sum((P - P_ref)**2)
    return S


def numerical_dSdx(x, xtal,  P_ref, eps=1e-4):
    if type(x) == list: x = np.array(x)

    xtal.update_from_1d_rep(x)
    atoms = xtal.to_ase() #* 2
    S0 = calculate_S(atoms, P_ref)
    dSdx = np.zeros(len(x))
    for i in range(len(x)):
        x0 = x.copy()
        x0[i] += eps
        xtal.update_from_1d_rep(x0)
        atoms = xtal.to_ase() #* 2
        S1 = calculate_S(atoms, P_ref)

        x0 = x.copy()
        x0[i] -= eps
        xtal.update_from_1d_rep(x0)
        atoms = xtal.to_ase() #* 2
        S2 = calculate_S(atoms, P_ref)

        dSdx[i] = 0.5*(S1-S2)/eps
    return dSdx


def calculate_dSdx_supercell(x, xtal, P_ref, eps=1e-4):

    xtal.update_from_1d_rep(x)
    atoms = xtal.to_ase() #* 2

    dPdr, P = calculator.compute_dpdr_5d(atoms)

    # Compute dSdr [N, M] [N, N, M, 3, 27] => [N, 3, 27]
    dSdr = np.einsum("ik, ijklm -> jlm", 2*(P - P_ref), dPdr)

    # Get supercell positions
    ref_pos = np.repeat(atoms.positions[:, :, np.newaxis], 27, axis=2)
    for cell in range(27):
        x1, y1, z1 = cell // 9 - 1, (cell // 3) % 3 - 1, cell % 3 - 1
        ref_pos[:, :, cell] += np.array([x1, y1, z1]) @ atoms.cell

    # Compute drdx via numerical func
    drdx = np.zeros([len(atoms), 3, 27, len(x)])

    xtal0 = xtal.copy()
    for i in range(len(x)):
        x0 = x.copy()
        x0[i] += eps
        xtal0.update_from_1d_rep(x0)
        atoms = xtal0.to_ase()

        # Get supercell positions
        pos = np.repeat(atoms.positions[:, :, np.newaxis], 27, axis=2)
        for cell in range(27):
            x1, y1, z1 = cell // 9 - 1, (cell // 3) % 3 - 1, cell % 3 - 1
            pos[:, :, cell] += np.array([x1, y1, z1]) @ atoms.cell

        drdx[:, :, :, i] += (pos - ref_pos)/eps

    # [N, 3, 27] [N, 3, 27, H] => H
    dSdx = np.einsum("ijk, ijkl -> l", dSdr, drdx)
    return dSdx


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

def get_dPdR_xtal(xtal, eps):
    p0 = calculator.calculate(xtal, derivative=True)
    shp = p0['x'].shape
    array1 = p0['dxdr']

    for j in range(shp[0]):
        for k in range(3):
            struc = get_perturbed_xtal(xtal, j, k, eps)
            p1 = calculator.calculate(struc)
            array2 = (p1['x'] - p0['x'])/eps
            #if np.linalg.norm(array2) > 1e-2: print(j, k, array2)
            if not np.allclose(array1[:, j, :, k], array2, atol=1e-4):
                print('xtal', j, k, '\n', array1[:, j, :, k], '\n', array2)
            assert(np.allclose(array1[:, j, :, k], array2, atol=1e-4))

# Descriptors Parameters
eps = 1e-8
rc1 = 1.8
rc2 = 3.5
nmax, lmax = 2, 2
calculator = SO3(nmax=nmax, lmax=lmax, rcut=rc1)
calculator0 = SO3(nmax=nmax, lmax=lmax, rcut=rc2)

# NaCl cluster
cluster = bulk('NaCl', crystalstructure='rocksalt', a=5.691694, cubic=True)
cluster = sort(cluster, tags=[0, 4, 1, 5, 2, 6, 3, 7])
cluster.set_pbc((0,0,0))
cluster = get_rotated_cluster(cluster, angle=0.1) # Must rotate

xtal = pyxtal()
xtal.from_prototype('graphite')
atoms = xtal.to_ase()
P_ref = calculator.compute_p(atoms)[0]

# Diamond
class TestCluster(unittest.TestCase):
    struc = get_rotated_cluster(cluster)
    p0 = calculator0.calculate(struc, derivative=True)
    struc = get_rotated_cluster(cluster, 10, 'x')
    p1 = calculator0.calculate(struc)

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
                p2 = calculator0.calculate(struc)
                array2 = (p2['x'] - self.p0['x'])/eps
                assert(np.allclose(array1[:,j,:,k], array2, atol=1e-3))

class TestXtal(unittest.TestCase):

    def test_dPdR_diamond(self):
        c = pyxtal()
        c.from_prototype('diamond')
        get_dPdR_xtal(c.to_ase(), eps)

    def test_dPdR_graphite(self):
        c = pyxtal()
        c.from_prototype('graphite')
        get_dPdR_xtal(c.to_ase(), eps)

    def test_dPdR_random(self):
        x = [ 7.952, 2.606, 0.592, 0.926, 0.608, 0.307]
        c = pyxtal()
        c.from_spg_wps_rep(179, ['6a', '6a', '6a', '6a'], x)
        get_dPdR_xtal(c.to_ase(), eps)

    def test_dPdR_random_P(self):
        x = [ 7.952, 2.606, 0.592, 0.926, 0.608, 0.307]
        c = pyxtal()
        c.from_spg_wps_rep(179, ['6a', '6a', '6a', '6a'], x)
        atoms = c.to_ase()
        p0 = calculator.compute_p(atoms)
        _, p1 = calculator.compute_dpdr(atoms)
        _, p2 = calculator.compute_dpdr_5d(atoms)
        assert(np.allclose(p0, p1, atol=1e-3))
        assert(np.allclose(p0, p2, atol=1e-3))

class TestSimilarity(unittest.TestCase):

    def test_sim_diamond(self):
        x = [3.0]
        c = pyxtal()
        c.from_spg_wps_rep(227, ['8a'], x, ['C'])
        atoms = c.to_ase()
        x = c.get_1d_rep_x()
        dSdx1 = numerical_dSdx(x, c, P_ref)
        dSdx2 = calculate_dSdx_supercell(x, c, P_ref)
        #print(dSdx1, dSdx2)
        assert(np.allclose(dSdx1, dSdx2, rtol=1e-1, atol=1e+1))

    def test_dPdR_random(self):
        #x = [ 7.952, 2.606, 0.592, 0.926, 0.608, 0.307]
        x = [9.55,  2.60,  0.48,  0.88,  0.76,   0.36]
        c = pyxtal()
        c.from_spg_wps_rep(179, ['6a', '6a', '6a', '6a'], x)
        atoms = c.to_ase()
        dSdx1 = numerical_dSdx(x, c, P_ref)
        dSdx2 = calculate_dSdx_supercell(x, c, P_ref)
        #print(dSdx1, dSdx2)
        assert(np.allclose(dSdx1, dSdx2, rtol=1e-1, atol=1e+1))

if __name__ == "__main__":
    unittest.main()
