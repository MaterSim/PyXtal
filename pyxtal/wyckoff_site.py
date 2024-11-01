"""
Module for handling Wyckoff sites for both atom and molecule
"""

# Standard Libraries
from copy import deepcopy
import numpy as np

# External Libraries
from pymatgen.core import Molecule
from scipy.spatial.distance import cdist
from scipy.spatial.transform import Rotation as R

from pyxtal.constants import rad
from pyxtal.database.element import Element
from pyxtal.lattice import Lattice
from pyxtal.operations import (
    SymmOp,
    create_matrix,
    distance_matrix,
    filtered_coords,
)
from pyxtal.symmetry import Group, Wyckoff_position
from pyxtal.tolerance import Tol_matrix


class atom_site:
    """
    Class for storing atomic Wyckoff positions with a single coordinate.

    Args:
        wp: a `WP <pyxtal.symmetry.Wyckoff_position.html> object
        coordinate (float): a fractional (x, y, z) coordinate
        specie (str): element name or symbol, or atomic number
        search (bool): whether or not search generator for special WP
    """

    def __init__(self, wp=None, coordinate=None, specie=1, search=False):
        self.position = np.array(coordinate)
        self.position -= np.floor(self.position)
        self.specie = Element(specie).short_name
        self.wp = wp
        self.coordination = None

        self._get_dof()
        self.PBC = self.wp.PBC
        self.multiplicity = self.wp.multiplicity
        if search:
            self.search_position()

        self.update()

    def __str__(self):
        if not hasattr(self.wp, "site_symm"):
            self.wp.get_site_symmetry()
        s = "{:>2s} @ [{:7.4f} {:7.4f} {:7.4f}], ".format(
            self.specie, *self.position)
        s += f"WP [{self.wp.get_label()}] "
        if self.coordination is not None:
            s += f" CN [{self.coordination:2d}] "
        s += "Site [{:}]".format(self.wp.site_symm.replace(" ", ""))

        return s

    def __repr__(self):
        return str(self)

    def copy(self):
        """
        Simply copy the structure
        """

        return deepcopy(self)

    def save_dict(self):
        return {
            "position": self.position,
            "specie": self.specie,
            "wp": self.wp.save_dict(),
        }

    def _get_dof(self):
        """
        get the number of dof for the given structures:
        """

        self.dof = self.wp.get_dof()

    def get_bounds(self):
        """
        get the number of dof for the given structures:
        """
        self.bounds = []
        for i in range(self.dof):
            self.bounds.append([0.0, 1.0])
        return self.bounds

    @classmethod
    def load_dict(cls, dicts):
        """
        load the sites from a dictionary
        """

        position = dicts["position"]
        specie = dicts["specie"]
        if "wp" in dicts:
            wp = Wyckoff_position.load_dict(dicts["wp"])
        else:
            hn, index = dicts["hn"], dicts["index"]
            wp = Wyckoff_position.from_group_and_index(
                hn, index, use_hall=True)

        return cls(wp, position, specie)

    def perturbate(self, lattice, magnitude=0.1):
        """
        Random perturbation of the site

        Args:
            lattice: lattice vectors
            magnitude: the magnitude of displacement (default: 0.1 A)
        """
        dis = (np.random.sample(3) - 0.5).dot(lattice)
        dis /= np.linalg.norm(dis)
        dis *= magnitude
        pos = self.position + dis.dot(np.linalg.inv(lattice))
        self.update(pos)

    def search_position(self, tol=1e-3):
        """
        Sometimes, the initial position is not the proper generator
        Needs to find the proper generator
        """

        found = False
        if self.wp.index > 0:
            wp0 = Group(self.wp.number, self.wp.dim)[0]
            for coord in wp0.apply_ops(self.position):
                coord -= np.floor(coord)
                ans = self.wp.ops[0].operate(coord)
                diff = coord - ans
                diff -= np.rint(diff)
                if np.sum(diff**2) < tol:
                    found = True
                    # print(found, coord, coord-ans)
                    self.position = coord - np.floor(coord)
                    break

        if not found:
            print("\nInput xyz", self.position)
            print("Target operation", self.wp.ops[0].as_xyz_str())
            raise ValueError("Cannot generate the desried generator")

    def encode(self):
        """
        transform dict to 1D vector
        [specie, wp.index, free x, y, z]
        """
        xyz = self.wp.get_free_xyzs(self.position)
        # print(self.wp.ops[0].rotation_matrix, self.wp.get_frozen_axis(), self.wp.get_dof())
        # print([self.specie, self.wp.index] + list(xyz))
        return [self.specie, self.wp.index, *list(xyz)]


    def shift_by_swap(self, swap_id):
        """
        check if a shift is needed during swap
        May occur for special WP in the I/A/B/C/F cases
        e.g., in space group 71 (Immm), the permutation
        4j(1/2,0,0.2) -> (0.2,0,1/2) -> 4f(0.7,1/2,0)
        it requires a shift of (0.5,0.5,0.5)
        """
        wp, shift = self.wp.swap_axis(swap_id)
        return shift

    def equivalent_set(self, tran, indices):
        """
        Transform the wp to another equivalent set.
        Needs to update both wp and positions

        Args:
            tran: affine matrix
            indices: the list of transformed wps
        """
        self.position = SymmOp(tran).operate(self.position)
        self.position -= np.floor(self.position)
        self.wp = self.wp.equivalent_set(
            indices[self.wp.index])  # update the wp index
        self.site_symm = self.wp.site_symm  # update the site symmetry
        self.update()

    def update(self, pos=None, reset_wp=False):
        """
        Used to generate coords from self.position
        """

        if pos is None:
            pos = self.position
        if reset_wp:
            self.wp.ops = Group(self.wp.number)[self.wp.index].ops
        self.coords = self.wp.apply_ops(pos)
        self.position = self.coords[0]

    def get_translations(self, pos, axis):
        """
        Get the displacement towards the reference positions

        Args:
            pos: reference position (1*3 vector)
            lattice: 3*3 matrix
            translation:
            axis:
        """

        # diffs0 = pos - self.coords
        diffs0 = self.wp.apply_ops(pos) - self.position
        diffs = diffs0.copy()
        diffs -= np.rint(diffs)
        diffs[:, axis] = 0
        return diffs0 - diffs

    def get_disp(self, pos, lattice, translation):
        """
        return the displacement towards the reference positions

        Args:
            pos: reference position (1*3 vector)
            lattice: 3*3 matrix
            translation:
        """
        coords = self.wp.apply_ops(pos)
        diffs = coords - (self.position + translation)
        # coords = self.wp.apply_ops(self.position + translation)
        # diffs = pos - coords

        diffs -= np.rint(diffs)
        dists = np.linalg.norm(diffs.dot(lattice), axis=1)
        id = np.argmin(dists)

        # print("++++", id, dists[id], id, diffs[id], translation) #; import sys; sys.exit()
        return diffs[id], dists[id]

    def check_with_ws2(self, ws2, lattice, tm, same_group=True):
        """
        Given two Wyckoff sites, checks the inter-atomic distances between them.

        Args:
            - ws2: a different WP object (will always return False if
            two identical WS's are provided)
            - lattice: a 3x3 cell matrix
            - same_group: whether or not two WS's are in the same structure.
            Default value True reduces the calculation cost
        Returns:
            True or False
        """

        # Ensure the PBC values are valid
        if self.PBC != ws2.PBC:
            raise ValueError("PBC values do not match between Wyckoff sites")

        # Get tolerance
        tol = tm.get_tol(self.specie, ws2.specie)
        # Symmetry shortcut method: check only some atoms
        if same_group is True:
            # We can either check one atom in WS1 against all WS2, or vice-versa
            # Check which option is faster
            if self.multiplicity > ws2.multiplicity:
                coords1 = [self.coords[0]]
                coords2 = ws2.coords
            else:
                coords1 = [ws2.coords[0]]
                coords2 = self.coords
            # Calculate distances
            dm = distance_matrix(coords1, coords2, lattice, PBC=self.PBC)
            # Check if any distances are less than the tolerance
            return not (dm < tol).any()
        # No symmetry method: check all atomic pairs
        else:
            dm = distance_matrix(
                self.wp.coords, ws2.coords, lattice, PBC=self.PBC)
            # Check if any distances are less than the tolerance
            return not (dm < tol).any()

    def substitute_with_single(self, ele):
        """
        chemical substitution with another element

        Args:
            ele (str): e.g. 'Zn'
        """
        self.specie = ele
        return self

    def substitute_with_linear(self, eles, direction, lattice):
        """
        chemical substitution with another linear building block, e.g. CN

        Args:
            eles (str): e.g., ['C', 'N']
            neighbors: two nearest neiboring atom xyz
        """
        bond_length = 1.4
        # direction = neighbors[1] - neighbors[0]
        # direction = np.dot(direction, lattice)
        # direction /= np.linalg.norm(direction)

        # Get the fraction coordinates
        shift = np.dot(bond_length / 2 * direction, np.linalg.inv(lattice))
        site1 = self.copy()
        site2 = self.copy()
        site1.specie = eles[0]
        site2.specie = eles[1]
        site1.update(site1.position + shift)
        site2.update(site2.position - shift)
        return site1, site2

    def to_mol_site(self, lattice, molecule, ori=None, reflect=False, type_id=0):
        """
        Transform it to the mol_sites, i.e., to build a molecule on
        the current WP
        """

        if ori is None:
            ori = [0, 0, 0]
        dicts = {}
        dicts["smile"] = molecule.smile
        dicts["type"] = type_id
        dicts["dim"] = 3
        dicts["PBC"] = [1, 1, 1]
        dicts["hn"] = self.wp.hall_number
        dicts["index"] = self.wp.index
        dicts["lattice"] = lattice.matrix
        dicts["lattice_type"] = lattice.ltype
        dicts["center"] = self.position
        if molecule.smile not in ["Cl-"]:
            dicts["orientation"] = np.array(ori)
            dicts["rotor"] = molecule.get_torsion_angles()
            dicts["reflect"] = reflect

        return mol_site.from_1D_dicts(dicts)


class mol_site:
    """
    Class for storing molecular Wyckoff positions and orientations within
    the `molecular_crystal` class. Each `mol_site` object represents an
    entire Wyckoff position, not necessarily a single molecule.
    This is the molecular version of the `Wyckoff_site` class.

    Attributes:
        mol (pyxtal_molecule): A pyxtal_molecule object representing the molecule at the site.
        position (list or array): A 3-vector representing the generating molecule's position
                                  in fractional coordinates.
        orientation (Orientation): An orientation object that describes the molecule's orientation.
        wp (Wyckoff_position): A Wyckoff position object that holds symmetry information.
        lattice (Lattice): A lattice object that defines the crystal structure's unit cell.
        stype (int): An integer specifying the type of molecule. Default is 0.
        symbols (list): List of atomic symbols for the atoms in the molecule.
        numbers (list): List of atomic numbers for the atoms in the molecule.
        PBC (list): Periodic boundary conditions inherited from the Wyckoff position.
        radius (float): Radius of the molecule, typically used for collision detection.
        tols_matrix (numpy array): Tolerance matrix for the molecular structure.

    Args:
        mol (pyxtal_molecule): A `pyxtal_molecule` object that describes the molecule.
        position (list or array): The 3D fractional coordinates of mol_center in the unit cell.
        orientation (Orientation): The orientation object describing the molecule's rotation.
        wp (Wyckoff_position): A `Wyckoff_position` object defining the symmetry of the site.
        lattice (Lattice, optional): The lattice of the crystal. Can be either a Lattice object
                                     or a matrix that will be converted into a Lattice object.
        stype (int, optional): Integer specifying the type of molecule. Default is 0.

    Methods:
        _get_dof(): Internal method to calculate the degrees of freedom (DoF) for the molecule.
    """

    def __init__(self, mol, position, orientation, wp, lattice=None, stype=0):
        # describe the molecule
        self.molecule = mol
        self.wp = wp
        self.position = position  # fractional coordinate of molecular center
        self.orientation = orientation  # pyxtal.molecule.orientation object
        self._get_dof()
        if isinstance(lattice, Lattice):
            self.lattice = lattice
        else:
            self.lattice = Lattice.from_matrix(lattice)
        self.PBC = self.wp.PBC
        # self.mol = mol.mol # A Pymatgen molecule object
        # [site.specie.value for site in self.mol.sites]
        self.symbols = mol.symbols
        self.numbers = self.molecule.mol.atomic_numbers
        self.tols_matrix = mol.tols_matrix
        self.radius = mol.radius
        self.type = stype

    def update_molecule(self, mol):
        self.molecule = mol
        self.numbers = mol.mol.atomic_numbers
        self.symbols = mol.symbols
        self.tols_matrix = mol.tols_matrix
        self.radius = mol.radius

    def update_orientation(self, angles):
        # QZ: Symmetrize the angle to the compatible orientation first
        self.orientation.r = R.from_euler("zxy", angles, degrees=True)
        self.orientation.matrix = self.orientation.r.as_matrix()

    def get_min_dist(self, angle=None):
        """
        Compute the minimum interatomic distance within the WP.

        Returns:
            minimum distance
        """
        if angle is not None:
            matrix = self.orientation.change_orientation(angle, update=False)
        else:
            matrix = self.orientation.matrix

        # Get total distances
        d0, _ = self.get_dists_auto(matrix=matrix)
        d1, _ = self.get_dists_WP(matrix=matrix, ignore=True)
        d_total = np.append(d0, d1, axis=0)

        return d_total.min()

    def optimize_orientation_by_dist(self, ori_attempts=10, verbose=False):
        """
        Optimize the orientation based on the shortest distance
        """
        # Set initial fun value and angle bounds
        ang_lo, ang_hi = 0, np.pi
        fun_lo = self.get_min_dist(ang_lo)
        fun_hi = self.get_min_dist(ang_hi)
        fun = fun_hi

        # Refine the orientation using a bisection method
        for _it in range(ori_attempts):

            # Return as soon as a good orientation is found
            if (fun > 0.8) and not self.short_dist():
                break

            # Compute the midpoint angle for bisection
            ang = (ang_lo + ang_hi) / 2
            fun = self.get_min_dist(ang)

            # Update based on the function value at the midpoint
            if fun_lo > fun_hi:
                ang_hi, fun_hi = ang, fun
            else:
                ang_lo, fun_lo = ang, fun

            if verbose: print("optimize_orientation_by_dist", _it, ang, fun)

        # Final adjustment
        if fun > fun_hi + 1e-4:
            self.orientation.change_orientation(ang_hi, update=True)
            fun = fun_hi
        elif fun > fun_lo + 1e-4:
            self.orientation.change_orientation(ang_lo, update=True)
            fun = fun_lo
        else:
            self.orientation.change_orientation(ang, update=True)

        return fun

    def get_energy(self, angle=None, wps=[], k1=1.0, r1=1.0, k2=2.0, k3=1.0, r2=6.0,
                   req=2.8, r_cut=2.0, verbose=False):
        """
        Compute the virtual energy based on Hbonding and repulsion

        Args:
            angle (float): rotation along a rotation axis
            wps (list): list of wps for consideration
        """
        # Set the orientation and compute the distance
        if angle is not None:
            matrix = self.orientation.change_orientation(angle, update=False)
        else:
            matrix = self.orientation.matrix

        coord_ref, _ = self._get_coords_and_species(first=True, unitcell=True, matrix=matrix)
        coord_ref = coord_ref @ self.lattice.matrix

        # Get total distances
        eng = 0
        d0, ids_0 = self.get_dists_auto(matrix=matrix, cutoff=8.0)
        d1, ids_1 = self.get_dists_WP(matrix=matrix, ignore=True, cutoff=8.0)
        d_total = np.append(d0, d1, axis=0)
        ids_total = np.append(ids_0, ids_1, axis=0)

        # initialize ids dicts
        acceptor_ids = {}; acceptor_ids[self.type] = self.molecule.active_sites[0]
        donor_ids = {}; donor_ids[self.type] = self.molecule.active_sites[1]
        H_ids = {}; H_ids[self.type] = self.molecule.active_sites[2]

        # Extract short distances ([x, y, z, d], [mol1, mol2, atom_id1, atom_id2])
        for wp2 in wps:
            d2, ids_2 = self.get_dists_WP2(wp2, matrix=matrix, ignore=True, cutoff=8.0)
            if len(d2) > 0:
                d_total = np.append(d_total, d2, axis=0)
                ids_total = np.append(ids_total, ids_2, axis=0)
                if wp2.type not in donor_ids.keys():
                    acceptor_ids[wp2.type] = wp2.molecule.active_sites[0]
                    donor_ids[wp2.type] = wp2.molecule.active_sites[1]
                    H_ids[wp2.type] = wp2.molecule.active_sites[2]

        # Count only 1 contribution per acceptor
        for id1 in self.molecule.active_sites[0]:
            rA = coord_ref[id1]
            dAD, angle, dAH = self.get_dist_angle_AD(d_total, ids_total, id1, donor_ids, H_ids, rA)
            eng += k2 * (dAD-req) ** 2
            eng += k3 * (angle-np.pi) ** 2
            #print(dAD, dAH, angle)
            if dAH < r_cut: eng -= k1 * dAH * np.exp(r_cut - dAH) - 1.0
            if verbose:
                print('Hbond AD', id1, dAD, dAH, np.degrees(angle), eng)

        # Count only 1 contribution per donor
        for i, id1 in enumerate(self.molecule.active_sites[1]):
            rD = coord_ref[id1]
            rH = coord_ref[self.molecule.active_sites[2][i]]

            dAD, angle, dAH = self.get_dist_angle_DA(d_total, ids_total, id1, acceptor_ids, rD, rH)
            eng += k2 * (dAD-req) ** 2
            eng += k3 * (angle-np.pi) ** 2
            if dAH < r_cut: eng -= k1 * dAH * np.exp(r_cut - dAH) - 1.0
            #print(dAD, dAH, angle)
            if verbose:
                print('Hbond DA', id1, dAD, dAH, np.degrees(angle), eng)

        # Count repulsion
        ds = d_total[:, 3]
        ds = ds[ds < r_cut]
        eng += k1 * (ds * np.exp(r_cut - ds) - 1.0).sum()
        #print(ds)
        if verbose:
            print('Repulsion', d_total[d_total < r_cut], eng)

        return eng

    def get_dist_angle_AD(self, d_total, ids_total, A_id, D_ids, H_ids, rA):
        """
        Args:
            d_total ():
            c_total ():
            A_id (int): index of acceptor atom
            D_ids ()
        """

        # Get satisfied
        rows0 = np.where((ids_total[:, 2]==A_id))[0]
        rows_D = []
        rows_H = []
        for row in rows0:
            m2 = ids_total[row, 1]
            if ids_total[row, 3] in D_ids[m2]:
                rows_D.append(row)
            elif ids_total[row, 3] in H_ids[m2]:
                rows_H.append(row)

        if len(rows_D) == 0:
            return 10.0, 0, 10.0
        else:
            myid = d_total[rows_D][:, 3].argmin()
            rD = d_total[rows_D][myid, :3]
            dAD = d_total[rows_D][myid, 3]

            # Find the shortest distance from H
            if len(rows_H) == 0:
                rH = rD + np.array([1, 0, 0])
            else:
                coords_H = d_total[rows_H]
                myid_H = coords_H[:, 3].argmin()
                rH = coords_H[myid_H][:3]

            rHD = rD - rH
            rHA = rA - rH
            dHD = np.linalg.norm(rHD)
            dHA = np.linalg.norm(rHA)
            if abs(dAD - np.linalg.norm(rD-rA)) > 1e-3:
                print("bug", dAD, rD, rA, np.linalg.norm(rD-rA))
                import sys; sys.exit()
            #print('debug dHA, dHD', dHA, dHD)
            cos = rHD @ rHA / (dHD * dHA)
            angle = np.arccos(np.clip(cos, -1.0, 1.0))

            return dAD, angle, dHA

    def get_dist_angle_DA(self, d_total, ids_total, D_id, A_ids, rD, rH):

        # Get satisfied
        rows0 = np.where((ids_total[:, 2]==D_id))[0]
        rows_A = []
        for row in rows0:
            m2 = ids_total[row, 1]
            if ids_total[row, 3] in A_ids[m2]:
                rows_A.append(row)

        if len(rows_A) == 0:
            return 10.0, 0, 10.0
        else:
            myid = d_total[rows_A][:, 3].argmin()
            rA = d_total[rows_A][myid, :3]
            dAD = d_total[rows_A][myid, 3]
            rHD = rD - rH
            rHA = rA - rH
            dHA = np.linalg.norm(rHA)
            dHD = np.linalg.norm(rHD)
            d_min = d_total[myid, 3]
            if abs(dAD - np.linalg.norm(rD-rA)) > 1e-3:
                print("bug", dAD, rD, rA, np.linalg.norm(rD-rA))
                import sys; sys.exit()
            #print('debug dHA, dHD', dHA, dHD)
            cos = rHD @ rHA / (dHD * dHA)
            angle = np.arccos(np.clip(cos, -1.0, 1.0))

            return dAD, angle, dHA

    def optimize_orientation_by_energy(self, wps=[], max_ax=20, max_ori=5, early_quit=3.0, verbose=False):
        """
        Iteratively optimize the orientation with the bisection method
        """
        for ax_trial in range(max_ax):

            # Select axis and compute the initial fun value and angle bounds
            # Select perpendicular????
            self.orientation.set_axis()
            ang_lo = 0 #self.orientation.angle
            ang_hi = np.pi #ang_lo + np.pi
            fun_lo = self.get_energy(wps=wps) #; print("call funlo", fun_lo)
            fun_hi = self.get_energy(ang_hi, wps=wps) #; print("call funhi", fun_hi)
            fun = fun_hi
            #if verbose: print("Init", ang_lo, fun_lo)

            # Refine the orientation using a bisection method
            for ori_trial in range(max_ori):

                # Compute the midpoint angle for bisection
                ang = (ang_lo + ang_hi) / 2
                fun = self.get_energy(ang, wps=wps)
                # Update based on the function value at the midpoint
                if fun_lo < fun_hi:
                    ang_hi, fun_hi = ang, fun
                else:
                    ang_lo, fun_lo = ang, fun
                #print("debug", ang_lo, ang_hi, fun_lo, fun_hi)
            # Finally pick the best one adjustment
            if fun > fun_hi + 1e-4:
                self.orientation.change_orientation(ang_hi, update=True)
                fun = fun_hi
            elif fun > fun_lo + 1e-4:
                self.orientation.change_orientation(ang_lo, update=True)
                fun = fun_lo
            else:
                self.orientation.change_orientation(ang, update=True)

            if verbose: print(f'Final {ax_trial:2d} {fun:.2f}')

            if fun <= early_quit:
                break

    def cut_lattice(self, ax, cut):
        """
        Cut lattice length on the given direction

        Args:
            ax (int): 0, 1, 2
            cut (float): the cut
        """
        paras = self.lattice.get_para()
        x0 = self.position[ax]
        x0 -= np.floor(x0)

        if x0 < 0.25:
            self.position[ax] = paras[ax] * x0 / (paras[ax]-cut)
        elif 0.25 <= x0 <= 0.75:
            self.position[ax] = (paras[ax] * x0 - 0.5 * cut) / (paras[ax]-cut)
        else:
            self.position[ax] = (paras[ax] * x0 - cut) / (paras[ax]-cut)
        #self.lattice.update_para(ax, -cut)

    def __str__(self):
        if not hasattr(self.wp, "site_symm"):
            self.wp.get_site_symmetry()

        self.angles = self.orientation.r.as_euler("zxy", degrees=True)
        formula = self.molecule.mol.formula.replace(" ", "")
        s = "{:12s} @ [{:7.4f} {:7.4f} {:7.4f}]  ".format(
            formula, *self.position)
        s += f"WP [{self.wp.get_label():s}] "
        s += "Site [{:}]".format(self.wp.site_symm.replace(" ", ""))
        if len(self.molecule.mol) > 1:
            s += " Euler [{:6.1f} {:6.1f} {:6.1f}]".format(*self.angles)

        return s

    def __repr__(self):
        return str(self)

    def _get_dof(self):
        """
        get the number of dof for the given wyckoff site:
        """
        dof = np.linalg.matrix_rank(self.wp.ops[0].rotation_matrix)
        self.dof = dof + 3
        if self.molecule.torsionlist is not None:
            self.dof += len(self.molecule.torsionlist)

    def get_bounds(self):
        """
        get the number of dof for the given structures:
        """
        self.bounds = []
        dof = np.linalg.matrix_rank(self.wp.ops[0].rotation_matrix)
        for i in range(dof):
            self.bounds.append([0.0, 1.0])
        self.bounds.append([-180.0, 180.0])
        self.bounds.append([-90.0, 90.0])
        self.bounds.append([-180.0, 180.0])
        if self.molecule.torsionlist is not None:
            for i in range(len(self.molecule.torsionlist)):
                self.bounds.append([0, 180.0])
        return self.bounds

    def save_dict(self):
        dict0 = {
            "position": self.position,
            "wp": self.wp.save_dict(),
            "molecule": self.molecule.save_str(),
            "orientation": self.orientation.save_dict(),
            "lattice": self.lattice.matrix,
            "lattice_type": self.lattice.ltype,
            "stype": self.type,
        }
        if self.molecule.torsionlist is not None:
            xyz, _ = self._get_coords_and_species(absolute=True, first=True)
            dict0["rotors"] = self.molecule.get_torsion_angles(xyz)

        return dict0

    @classmethod
    def load_dict(cls, dicts):
        """
        load the sites from a dictionary
        """
        from pyxtal.molecule import Orientation, pyxtal_molecule

        mol = pyxtal_molecule.load_str(dicts["molecule"])
        position = dicts["position"]
        orientation = Orientation.load_dict(dicts["orientation"])
        wp = Wyckoff_position.load_dict(dicts["wp"])
        lattice = Lattice.from_matrix(
            dicts["lattice"], ltype=dicts["lattice_type"])
        stype = dicts["stype"]
        return cls(mol, position, orientation, wp, lattice, stype)

    def encode(self):
        """
        transform dict to 1D vector
        [wp_id, x, y, z, or1, or2, or3, rotor1, rotor2, .etc]
        """
        if len(self.molecule.mol) > 1:
            xyz, _ = self._get_coords_and_species(absolute=True, first=True)
            # if len(xyz)==3: print("encode: \n", self.molecule.mol.cart_coords)
            rotor = self.molecule.get_torsion_angles(xyz)
            ori, _, reflect = self.molecule.get_orientation(xyz)
            return [self.wp.index, *list(self.position), *list(ori), *rotor, reflect]
        else:
            return [self.wp.index, *list(self.position), 0]

    def to_1D_dicts(self):
        """
        save the wp in 1D representation
        """
        xyz, _ = self._get_coords_and_species(absolute=True, first=True)
        dict0 = {"smile": self.molecule.smile}
        dict0["rotor"] = self.molecule.get_torsion_angles(xyz)
        dict0["orientation"], dict0["rmsd"], dict0["reflect"] = self.molecule.get_orientation(
            xyz)
        dict0["rotor"]
        # rdkit_mol = self.molecule.rdkit_mol(self.molecule.smile)
        # conf0 = rdkit_mol.GetConformer(0)
        # print(self.molecule.set_torsion_angles(conf0, angs))
        # print("save matrix"); print(self.orientation.r.as_matrix())
        # print("save angle"); print(self.orientation.r.as_euler('zxy', degrees=True))
        # print("angle"); print(dict0["orientation"])
        # self.molecule.get_center(xyz)
        dict0["center"] = self.position - np.floor(self.position)
        dict0["number"] = self.wp.number
        dict0["index"] = self.wp.index
        dict0["PBC"] = self.wp.PBC
        dict0["dim"] = self.wp.dim
        dict0["lattice"] = self.lattice.matrix
        dict0["lattice_type"] = self.lattice.ltype

        return dict0

    @classmethod
    def from_1D_dicts(cls, dicts):
        from pyxtal.molecule import Orientation, pyxtal_molecule

        mol = pyxtal_molecule(mol=dicts["smile"] + ".smi", fix=True)
        if len(mol.mol) > 1:
            if len(dicts["smile"]) > 1:
                conf = mol.rdkit_mol().GetConformer(0)
                if dicts["reflect"]:
                    mol.set_torsion_angles(conf, dicts["rotor"], False)
                xyz = mol.set_torsion_angles(
                    conf, dicts["rotor"], dicts["reflect"])
            else:
                # for H2O, use the standard one
                xyz = np.array(
                    [
                        [-0.00111384, 0.36313718, 0.0],
                        [-0.82498189, -0.18196256, 0.0],
                        [0.82609573, -0.18117463, 0.0],
                    ]
                )

            mol.reset_positions(xyz)
            matrix = R.from_euler(
                "zxy", dicts["orientation"], degrees=True).as_matrix()
            orientation = Orientation(matrix)
        else:
            orientation = Orientation(np.eye(3))

        g = dicts["hn"]
        index = int(dicts["index"])
        dim = dicts["dim"]
        wp = Wyckoff_position.from_group_and_index(g, index, dim, dicts["PBC"])
        lattice = Lattice.from_matrix(
            dicts["lattice"], ltype=dicts["lattice_type"])
        # np.dot(dicts["center"], lattice.inv_matrix)
        position = dicts["center"]
        position, wp, _ = wp.merge(position, lattice.matrix, 0.01)

        return cls(mol, position, orientation, wp, lattice)

    def show(self, id=None, **kwargs):
        """
        display WP on the notebook
        """
        from pyxtal.viz import display_molecular_site

        return display_molecular_site(self, id, **kwargs)

    def _get_coords_and_species(self, absolute=False, PBC=False,
                                first=False, unitcell=False,
                                matrix=None):
        """
        Used to generate coords and species for get_coords_and_species

        Args:
            absolute: return absolute or relative coordinates
            PBC: whether or not to add coordinates in neighboring unit cells,
            first: whether or not to extract the information from only the first site
            unitcell: whether or not to move the molecular center to the unit cell
            matrix: orientatin matrix

        Returns:
            atomic coords: a numpy array of atomic coordinates in the site
            species: a list of atomic species for the atomic coords
        """
        if matrix is None: matrix = self.orientation.matrix
        coord0 = self.molecule.mol.cart_coords.dot(matrix.T)
        wp_atomic_sites = []
        wp_atomic_coords = None

        for point_index, op2 in enumerate(self.wp.ops):
            # Obtain the center in absolute coords
            center_relative = op2.operate(self.position)
            if unitcell:
                center_relative -= np.floor(center_relative)
            center_absolute = np.dot(center_relative, self.lattice.matrix)

            # Rotate the molecule (Euclidean metric)
            # op2_m = self.wp.generators_m[point_index]
            op2_m = self.wp.get_euclidean_generator(
                self.lattice.matrix, point_index)
            rot = op2_m.affine_matrix[:3, :3].T
            # NOTE=====the euclidean_generator has wrong translation vectors,
            # but we don't care. This needs to be fixed later
            tmp = np.dot(coord0, rot)  # + tau

            # Add absolute center to molecule
            tmp += center_absolute
            tmp = tmp.dot(self.lattice.inv_matrix)
            wp_atomic_coords = tmp if wp_atomic_coords is None else np.append(
                wp_atomic_coords, tmp, axis=0)
            wp_atomic_sites.extend(self.symbols)
            if first:
                break

        if PBC:
            # Filter PBC of wp_atomic_coords
            wp_atomic_coords = filtered_coords(wp_atomic_coords, PBC=self.PBC)
            # Add PBC copies of coords
            m = create_matrix(PBC=self.PBC, omit=True)
            # Move [0,0,0] PBC vector to first position in array
            m2 = [[0, 0, 0], *list(m)]
            new_coords = np.vstack([wp_atomic_coords + v for v in m2])
            wp_atomic_coords = new_coords

        if absolute:
            wp_atomic_coords = wp_atomic_coords.dot(self.lattice.matrix)

        return wp_atomic_coords, wp_atomic_sites

    def get_coords_and_species(self, absolute=False, PBC=False, unitcell=False):
        """
        Lazily generates and returns the atomic coordinate and species for the
        Wyckoff position. Plugs the molecule into the provided orientation
        (with angle=0), and calculates the new positions.

        Args:
            absolute: return absolute or relative coordinates
            PBC: whether or not to add coordinates in neighboring unit cells
            unitcell: whether or not to move the molecule center to the unit cell

        Returns:
            coords: a np array of 3-vectors.
            species: a list of atomic symbols, e.g. ['H', 'H', 'O', 'H', 'H', 'O']
        """
        return self._get_coords_and_species(absolute, PBC, unitcell=unitcell)

    def perturbate(self, lattice, trans=0.1, rot=5):
        """
        Random perturbation of the molecular site

        Args:
            lattice: lattice vectors
            trans: magnitude of tranlation vectors (default: 0.1 A)
            rot: magnitude of rotation degree (default: 5.0)
        """
        dis = (np.random.sample(3) - 0.5).dot(lattice)
        dis /= np.linalg.norm(dis)
        dis *= trans
        self.translate(dis, True)
        if rot == "random":
            self.orientation.change_orientation()
        else:
            self.orientation.change_orientation(angle=rot / 180 * np.pi)

    def translate(self, disp=np.zeros(3), absolute=False):
        """
        To translate the molecule
        """
        disp = np.array(disp)
        if absolute:
            disp = disp.dot(self.lattice.inv_matrix)
        position = self.position + disp
        self.position = self.wp.project(position)

    def rotate(self, ax_id=0, ax_vector=None, angle=180):
        """
        To rotate the molecule
        Args:
            ax_id: the principle axis id
            ax_vector (float): 3-vector to define the axis
            angle (float): angle to rotate
        """
        p = self.orientation.r

        if ax_vector is not None:
            ax = ax_vector / np.linalg.norm(ax_vector)
        else:
            xyz = self.molecule.mol.cart_coords.dot(p.as_matrix().T)
            ax = self.molecule.get_principle_axes(xyz).T[ax_id]

        q = R.from_rotvec(ax * rad * angle)
        o = q * p
        self.orientation.r = o
        self.orientation.matrix = o.as_matrix()

    # def is_compatible_symmetry(self, tol=0.3):
    #    """
    #    Check if the molecular symmetry matches the site symmetry
    #    """
    #    mol = self.molecule.mol
    #    if len(mol) == 1 or self.wp.index==0:
    #        return True
    #    else:
    #        pga = PointGroupAnalyzer(mol, tol)
    #        for op in self.wp.symmetry_m[0]:
    #            if not pga.is_valid_op(op):
    #                return False
    #        return True

    def is_valid_orientation(self):
        pass

    def get_mol_object(self, id=0):
        """
        make the pymatgen molecule object

        Args:
            id: the index of molecules in the given site

        Returns:
            a molecule object
        """
        coord0 = self.molecule.mol.cart_coords
        coord0 = coord0 @ self.orientation.matrix.T

        # Obtain the center in absolute coords
        if not hasattr(self.wp, "generators"):
            self.wp.set_generators()

        if id <= len(self.wp.generators):
            # op = self.wp.generators[id]
            op = self.wp.get_euclidean_generator(self.lattice.matrix, id)
            center_relative = op.operate(self.position)
            center_relative -= np.floor(center_relative)
            center_absolute = np.dot(center_relative, self.lattice.matrix)

            # Rotate the molecule (Euclidean metric)
            # op_m = self.wp.generators_m[id]
            # rot = op_m.affine_matrix[0:3][:, 0:3].T
            # tau = op_m.affine_matrix[0:3][:, 3]
            op0 = self.wp.get_euclidean_generator(self.lattice.matrix, id)
            rot = op0.rotation_matrix.T
            tmp = np.dot(coord0, rot)
            # Add absolute center to molecule
            tmp += center_absolute
            return Molecule(self.symbols, tmp)
        else:
            raise ValueError("id is greater than the number of molecules")

    def show_molecule_in_box(self, id=0):
        """
        display molecule with box

        Args:
            id (int): molecular id
        """
        from pyxtal.viz import display_molecule

        mol = self.get_mol_object(id)
        cell, vertices, center = self.molecule.get_box_coordinates(
            mol.cart_coords)
        return display_molecule(mol, center, cell)

    def update(self, coords=None, lattice=None, absolute=False, update_mol=True):
        """
        After the geometry relaxation, the returned atomic coordinates
        maybe rescaled to [0, 1] bound. In this case, we need to refind
        the molecular coordinates according to the original neighbor list.
        If the list does not change, we return the new coordinates
        otherwise, terminate the calculation.

        Args:
            coords (float): input xyz coordinates for the given molecule
            absolute (bool): whether the input is absolute coordinates
            update_mol (bool): whether update the molecule (used for substitution)
        """
        from pyxtal.molecule import Orientation, compare_mol_connectivity

        try:
            from openbabel import openbabel, pybel
        except:
            import openbabel
            import pybel

        if lattice is not None:
            self.lattice = lattice

        if coords is None:
            coords, _ = self._get_coords_and_species(
                absolute=False, first=True, unitcell=True)
            absolute = False

        if not absolute:
            coords = coords.dot(self.lattice.matrix)

        # mol = Molecule(self.symbols, coords-np.mean(coords, axis=0))
        center = self.molecule.get_center(coords)
        mol = Molecule(self.symbols, coords - center)

        # match, _ = compare_mol_connectivity(mol, self.mol, True)
        # Update the orientation matrix
        match, _ = compare_mol_connectivity(mol, self.molecule.mol)
        if match:
            # position = np.mean(coords, axis=0).dot(self.lattice.inv_matrix)
            position = center.dot(self.lattice.inv_matrix)
            self.position = position - np.floor(position)
            if update_mol:
                self.orientation = Orientation(np.eye(3))
                self.molecule.mol = mol
            else:
                m1 = pybel.readstring("xyz", self.molecule.mol.to("xyz"))
                m2 = pybel.readstring("xyz", mol.to("xyz"))
                aligner = openbabel.OBAlign(True, False)
                aligner.SetRefMol(m1.OBMol)
                aligner.SetTargetMol(m2.OBMol)
                if aligner.Align():
                    print("RMSD: ", aligner.GetRMSD())
                    rot = np.zeros([3, 3])
                    for i in range(3):
                        for j in range(3):
                            rot[i, j] = aligner.GetRotMatrix().Get(i, j)
                    if abs(np.linalg.det(rot) - 1) < 1e-2:
                        self.orientation.matrix = rot
                        self.orientation.r = R.from_matrix(rot)
                    else:
                        raise ValueError("rotation matrix is wrong")
        else:
            import pickle

            with open("wrong.pkl", "wb") as f:
                pickle.dump([mol, self.molecule.mol], f)
                mol.to(filename="Wrong.xyz", fmt="xyz")
                self.molecule.mol.to(filename="Ref.xyz", fmt="xyz")
            raise ValueError("molecular connectivity changes! Exit")
        # todo check if connectivty changed

    def _create_matrix(self, center=False, ignore=False):
        """
        Used for calculating distances in lattices with periodic boundary
        conditions. When multiplied with a set of points, generates additional
        points in cells adjacent to and diagonal to the original cell

        Args:
            center:
            ignore:

        Returns:
            A numpy array of matrices which can be multiplied by a set of
            coordinates
        """
        abc = [self.lattice.a, self.lattice.b, self.lattice.c]

        # QZ: This should be based on the occupation of current molecule
        if hasattr(self, "ijk_lists"):
            ijk_lists = self.ijk_lists
        else:

            def get_ijk_range(pbc, abc_val, ignore, radius):
                if not pbc:
                    return [0]
                if not ignore and abc_val > 50.0 and radius < 10:
                    return [0]
                if abc_val < 8.5:
                    return list(range(-3, 4))
                return [-1, 0, 1]

            ijk_lists = [get_ijk_range(
                self.PBC[idx], abc[idx], ignore, self.radius) for idx in range(3)]

        matrix = [[0, 0, 0]] if center else []
        matrix += [
            [i, j, k] for i in ijk_lists[0] for j in ijk_lists[1] for k in ijk_lists[2] if [i, j, k] != [0, 0, 0]
        ]

        # In case a,b,c are all greater than 20
        if len(matrix) == 0:
            matrix = [[1, 0, 0]]
        return np.array(matrix, dtype=float)

    def get_distances(self, coord1, coord2, m2=None, center=True, ignore=False):
        """
        Compute the distance matrix between the central molecule (coord1) and
        neighboring molecules (coord2) under the periodic boundary condition.

        Args:
            coord1 (numpy array): Fractional coordinates of the central molecule.
                                Shape: (m1, 3), where m1 is the number of atoms
            coord2 (numpy array): Fractional coordinates of the neighboring molecules.
                                Shape: (N2*m2, 3), where N2 is the number of atoms
                                and m2 is the number of atoms in each neighboring molecule.
            m2 (int, optional): N_atoms in each neighboring molecule. If not provided,
                                it's assumed to be equal m1.
            center (bool, optional): If `True`, count self-image of the reference molecule
            ignore (bool, optional): If `True`, ignores some periodic boundary conditions.

        Returns:
            distance matrix: [m1*m2*pbc, m1, m2]
            coord2 under PBC: [pbc, m2, 3]
        """
        m1 = len(coord1)
        if m2 is None: m2 = m1
        N2 = int(len(coord2) / m2) # Number of molecule 2

        # peridoic images
        m = self._create_matrix(center, ignore)  # PBC matrix
        coord2 = np.vstack([coord2 + v for v in m])
        N = N2 * len(m) # Number of PBC images

        # absolute xyz
        coord1 = coord1 @ self.lattice.matrix
        coord2 = coord2 @ self.lattice.matrix

        d = cdist(coord1, coord2)
        d = d.reshape(m1, N, m2).transpose(1, 0, 2)
        coord2 = coord2.reshape([N, m2, 3])
        return d, coord2

    def get_dists_auto(self, matrix=None, ignore=False, cutoff=None):
        """
        Compute the distances between the periodic images
        N: number of atoms in the molecule
        M:

        Args:
            ignore (bool, optional): If `True`, ignores some periodic boundary conditions.
            cutoff (): if not None, reduce the distance matrix

        Returns:
            a distance matrix (M, N, N)
            list of molecular xyz (M, N, 3)
        """
        coord1, _ = self._get_coords_and_species(first=True, unitcell=True, matrix=matrix)
        ds, coords = self.get_distances(coord1, coord1, center=False, ignore=ignore)
        if cutoff is None:
            return ds, coords
        else:
            return self.extract_short_distances(ds, coords, cutoff, self.type, self.type)

    def extract_short_distances(self, ds, coords, cutoff, label1, label2):

        indices = np.where( ds < cutoff )
        d1s = np.zeros((len(indices[0]), 4))
        ids = np.zeros((len(indices[0]), 4), dtype=int)
        for i in range(len(indices[0])):
            m, n1, n2 = indices[0][i], indices[1][i], indices[2][i]
            d1s[i, :3] = coords[m, n2, :]
            d1s[i, 3] = ds[m, n1, n2]
            ids[i] = [label1, label2, n1, n2]
        return d1s, ids


    def get_dists_WP(self, matrix=None, ignore=False, idx=None, cutoff=None):
        """
        Compute the distances within the WP site

        Returns:
            a distance matrix (M, N, N)
            list of molecular xyz (M, N, 3)
        """
        m_length = len(self.symbols)
        coords, _ = self._get_coords_and_species(unitcell=True, matrix=matrix)
        coord1 = coords[:m_length]  # 1st molecular coords
        if idx is None:
            coord2 = coords[m_length:]
        else:
            coord2 = coords[m_length *(idx): m_length * (idx + 1)]
        ds, coords = self.get_distances(coord1, coord2, ignore=ignore)
        if cutoff is None:
            return ds, coords
        else:
            return self.extract_short_distances(ds, coords, cutoff, self.type, self.type)

    def get_dists_WP2(self, wp2, matrix=None, ignore=False, cutoff=None):
        """
        Compute the distances w.r.t to other WP site

        Returns:
            a distance matrix (M, N, N)
            list of molecular xyz (M, N, 3)
        """
        coord1, _ = self._get_coords_and_species(first=True, unitcell=True, matrix=matrix)
        coord2, _ = wp2._get_coords_and_species(unitcell=True)
        ds, coords = self.get_distances(coord1, coord2, len(wp2.numbers), ignore=ignore)
        if cutoff is None:
            return ds, coords
        else:
            return self.extract_short_distances(ds, coords, cutoff, self.type, wp2.type)

    def short_dist(self):
        """
        Check if the atoms are too close within the WP.

        Returns:
            True or False
        """
        tols_matrix = self.tols_matrix
        # Check periodic images
        d, _ = self.get_dists_auto()
        if np.min(d) < np.max(tols_matrix):
            tols = np.min(d, axis=0)
            if (tols < tols_matrix).any():
                return True

        if self.wp.multiplicity > 1:
            d, _ = self.get_dists_WP()
            if np.min(d) < np.max(tols_matrix):
                tols = np.min(d, axis=0)  # N*N matrix
                if (tols < tols_matrix).any():
                    return True

        return False

    def short_dist_with_wp2(self, wp2, tm=Tol_matrix(prototype="molecular")):
        """
        Check whether or not the molecules of two wp sites overlap. Uses
        ellipsoid overlapping approximation to check.

        Args:
            wp2: the 2nd wp sites
            tm: a Tol_matrix object (or prototype string) for distance checking
        Returns:
            True or False
        """

        # Get coordinates for both mol_sites
        c1, _ = self.get_coords_and_species()
        c2, _ = wp2.get_coords_and_species()
        m_length1 = len(self.numbers)
        m_length2 = len(wp2.numbers)

        # choose wp with bigger molecule as the center
        if len(c2) <= len(c1):
            coord1 = c1[:m_length1]
            coord2 = c2  # rest molecular coords
            tols_matrix = self.molecule.get_tols_matrix(wp2.molecule, tm)
            m2 = m_length2
        else:
            coord1 = c2[:m_length2]
            coord2 = c1
            tols_matrix = wp2.molecule.get_tols_matrix(self.molecule, tm)
            m2 = m_length1

        # compute the distance matrix
        d, _ = self.get_distances(
            coord1 - np.floor(coord1), coord2 - np.floor(coord2), m2)
        # print("short dist", len(c1), len(c2), d.min())

        if np.min(d) < np.max(tols_matrix):
            tols = np.min(d, axis=0)
            if (tols < tols_matrix).any():
                return False
        return True

    def get_neighbors_auto(self, factor=1.1, max_d=4.0, ignore_E=True, detail=False, etol=-5e-2):
        """
        Find neighboring molecules within a given distance threshold.

        The function identifies neighboring molecules within PBC and computes the shortest
        distances and (optionally) interaction energies. it returns detailed information
        about neighboring molecule pairs, distances, and energies.

        Args:
            factor (float, optional): Scaling factor for distance (default is 1.1).
            max_d (float, optional): Maximum distance for neighbors (default is 4.0 Å).
            ignore_E (bool, optional): Skips energy calculations (default is `True`).
            detail (bool, optional): Returns detailed eng, molecular pairs and
                                     distances instead of only shortest distances.
            etol (float, optional): Energy tolerance for filtering pairs in detailed
                                    mode (default is -5e-2).

        Returns:
            If `detail is True`:
                engs (list): List of interaction energies for valid molecular pairs.
                pairs (list): List of tuples with neighbor molecules and positions.
                dists (list): List of distances between neighboring molecular pairs.
            If `detail is False`:
                min_ds (list): List of shortest distances between central and neighbors.
                neighs (list): List of neighboring molecular coordinates with PBC.
                Ps (list): List of Wyckoff position multiplicities or translations.
                engs (list): List of energies, if energy calculation is not skipped.
        """

        # Compute the mol_center in Cartesian coordinate
        position = self.position - np.floor(self.position)
        mol_center = np.dot(position, self.lattice.matrix)

        # Atomic numbers for atoms in the central molecule
        numbers = self.molecule.mol.atomic_numbers

        # Get fractional coordinates for the central molecule
        coord1 = self._get_coords_and_species(first=True, unitcell=True)[0]

        # Initialize tolerance matrix for intermolecular distances
        tm = Tol_matrix(prototype="vdW", factor=factor)
        tols_matrix = self.molecule.get_tols_matrix(tm=tm)

        # Initialize coefficient matrix for energy calculations if needed
        coef_matrix = None
        if not ignore_E:
            coef_matrix = self.molecule.get_coefs_matrix()
            if coef_matrix is not None:
                A = coef_matrix[:, :, 0]
                B = coef_matrix[:, :, 1]
                C = coef_matrix[:, :, 2]

        # Initialize lists for results
        min_ds = []
        neighs = []
        Ps = []
        engs = []
        pairs = []
        dists = []

        # Find neighbors under PBC
        d, coord2 = self.get_dists_auto(ignore=True)

        # Loop through each neighboring molecule
        for i in range(d.shape[0]):
            if np.min(d[i]) < max_d and (d[i] < tols_matrix).any():
                if coef_matrix is not None:
                    if detail:
                        eng = A * np.exp(-B * d[i]) - C / (d[i] ** 6)
                        ids = np.where(eng < etol)
                        for id in range(len(ids[0])):
                            n1, n2 = numbers[ids[0][id]], numbers[ids[1][id]]
                            if 1 not in [n1, n2]:
                                pos = coord2[i][ids[1][id]] - mol_center
                                pairs.append((n2, pos))
                                engs.append(eng[ids[0][id], ids[1][id]])
                                dists.append(d[i][ids[0][id], ids[1][id]])
                    else:
                        eng0 = A * np.exp(-B * d[i]) - C / (d[i] ** 6)
                        engs.append(eng0.sum())
                else:
                    engs.append(None)
                    if detail:
                        ids = np.where(d[i] < max_d)
                        for id in range(len(ids[0])):
                            n1, n2 = numbers[ids[0][id]], numbers[ids[1][id]]
                            if 1 not in [n1, n2]:
                                pos = coord2[i][ids[1][id]] - mol_center
                                pairs.append((n2, pos))
                                dists.append(np.linalg.norm(pos))

                tmp = d[i] / tols_matrix
                _d = tmp[tmp < 1.0]
                id = np.argmin(tmp.flatten())
                d[i].flatten()[id]
                min_ds.append(min(_d) * factor)
                neighs.append(coord2[i])
                Ps.append(0)

        # Handle Wyckoff position multiplicities (if applicable)
        if self.wp.multiplicity > 1:
            for idx in range(1, self.wp.multiplicity):
                P = 0 if self.wp.is_pure_translation(idx) else 1
                d, coord2 = self.get_dists_WP(ignore=True, idx=idx)
                for i in range(d.shape[0]):
                    if np.min(d[i]) < max_d and (d[i] < tols_matrix).any():
                        if coef_matrix is not None:
                            if detail:
                                eng = A * np.exp(-B * d[i]) - C / (d[i] ** 6)
                                ids = np.where(eng < etol)
                                for id in range(len(ids[0])):
                                    n1, n2 = numbers[ids[0][id]
                                                     ], numbers[ids[1][id]]
                                    if 1 not in [n1, n2]:
                                        pos = coord2[i][ids[1]
                                                        [id]] - mol_center
                                        pairs.append((n2, pos))
                                        engs.append(
                                            eng[ids[0][id], ids[1][id]])
                                        dists.append(
                                            d[i][ids[0][id], ids[1][id]])

                            else:
                                eng0 = A * np.exp(-B * d[i]) - C / (d[i] ** 6)
                                engs.append(eng0.sum())
                        else:
                            engs.append(None)
                            if detail:
                                # print('OMMMM', d[i].min())
                                ids = np.where(d[i] < max_d)
                                # for id in zip(*ids):
                                for id in range(len(ids[0])):  # zip(*ids):
                                    n1, n2 = numbers[ids[0][id]
                                                     ], numbers[ids[1][id]]
                                    if 1 not in [n1, n2]:  # != [1, 1]:
                                        pos = coord2[i][ids[1]
                                                        [id]] - mol_center
                                        pairs.append((n2, pos))
                                        dists.append(np.linalg.norm(pos))

                        tmp = d[i] / tols_matrix
                        _d = tmp[tmp < 1]
                        id = np.argmin(tmp.flatten())
                        d[i].flatten()[id]
                        min_ds.append(min(_d) * factor)
                        neighs.append(coord2[i])
                        Ps.append(P)

        # Return results based on the detail flag
        if detail:
            return engs, pairs, dists
        else:
            return min_ds, neighs, Ps, engs

    def get_neighbors_wp2(self, wp2, factor=1.1, max_d=4.0, ignore_E=True, detail=False, etol=-5e-2):
        """
        Find the neigboring molecules from a 2nd wp site

        Returns
            min_ds: list of shortest distances
            neighs: list of neighboring molecular xyzs
        """

        tm = Tol_matrix(prototype="vdW", factor=factor)
        m_length1 = len(self.numbers)
        m_length2 = len(wp2.numbers)

        # Get coordinates for both mol_sites
        c1, _ = self.get_coords_and_species()
        c2, _ = wp2.get_coords_and_species()

        coord1 = c1[:m_length1]
        coord2 = c2  # rest molecular coords
        tols_matrix = self.molecule.get_tols_matrix(wp2.molecule, tm)
        coef_matrix = None
        if not ignore_E:
            coef_matrix = self.molecule.get_coefs_matrix(wp2.molecule)
            if coef_matrix is not None:
                A = coef_matrix[:, :, 0]
                B = coef_matrix[:, :, 1]
                C = coef_matrix[:, :, 2]

        # compute the distance matrix
        d, coord2 = self.get_distances(coord1, coord2, m_length2, ignore=True)
        min_ds = []
        neighs = []
        engs = []
        dists = []
        pairs = []

        for i in range(d.shape[0]):
            if np.min(d[i]) < max_d and (d[i] < tols_matrix).any():
                if coef_matrix is not None:
                    if detail:
                        eng = A * np.exp(-B * d[i]) - C / (d[i] ** 6)
                        ids = np.where(eng < etol)
                        for id in zip(*ids):
                            tmp1, tmp2 = coord1[id[0]], coord2[i][id[1]]
                            pairs.append((tmp1 + tmp2) / 2)
                            engs.append(eng[id])
                            dists.append(d[i][id])
                    else:
                        eng0 = A * np.exp(-B * d[i]) - C / (d[i] ** 6)
                        engs.append(eng0.sum())
                else:
                    engs.append(None)
                tmp = d[i] / tols_matrix
                _d = tmp[tmp < 1]
                id = np.argmin(tmp.flatten())
                d[i].flatten()[id]
                min_ds.append(min(_d) * factor)
                neighs.append(coord2[i])

        if detail:
            return engs, pairs, dists
        else:
            return min_ds, neighs, engs

    def get_ijk_lists(self, value=None):
        """
        Get the occupatation in the unit cell for the generating molecule
        This can be used to estimate the supercell size for finding neighbors

        Returns
            PBC
        """

        # Get coordinates for both mol_sites
        if value is None:
            ijk_lists = []
            c1, _ = self.get_coords_and_species(absolute=False)
            for id in range(3):
                if self.PBC[id]:
                    m1 = np.min(c1[:, id])
                    m2 = np.max(c1[:, id])
                    max_id = int(np.ceil(2 * m2 - m1))
                    min_id = int(np.floor(-m2))
                    ijk_lists.append(list(range(min_id - 1, max_id + 1)))
                else:
                    ijk_lists.append([0])
            self.ijk_lists = ijk_lists
        else:
            self.ijk_lists = value

    def to_atom_site(self, specie=1):
        """
        transform it to the mol_sites, i.e., to build a molecule on
        the current WP
        """
        dicts = {}
        dicts["specie"] = specie
        dicts["position"] = self.position
        dicts["hn"] = self.wp.hall_number
        dicts["index"] = self.wp.index
        return atom_site.load_dict(dicts)


if __name__ == "__main__":

    from pyxtal.symmetry import Wyckoff_position as WP

    wp = WP.from_group_and_letter(225, 'a')
    coordinate = [0.25, 0.25, 0.25]
    specie = 6
    atom_site_instance = atom_site(wp=wp,
                                   coordinate=coordinate,
                                   specie=specie)

    # Print the created instance
    print(atom_site_instance)
