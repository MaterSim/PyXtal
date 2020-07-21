"""
Module defining custom classes for symmetric structures. We store the Wyckoff position
information needed to re-generate a randomly generated crystal, or for optimizing a structure
without breaking symmetry.
"""
# Imports
# ------------------------------
# External libraries
from pymatgen.core.structure import Structure
import numpy as np

# Define functions
# ------------------------------
class Xstruct(Structure):
    """
    Class for storing an atomic crystal and its symmetry information.
    The Wyckoff_site's are stored in wyckoff_sites. This gives access to
    the generating coordinates for each Wyckoff position in the crystal.


    """

    def __init__(self, wyckoff_sites, lattice, group):
        self.wyckoff_sites = wyckoff_sites
        self.group = group
        # Extract coords and species
        species = []
        coords = []
        for ws in wyckoff_sites:
            species += [ws.specie] * ws.multiplicity
            for c in ws.coords:
                coords.append(c)
        coords = np.array(coords)
        Structure.__init__(self, lattice.matrix, species, coords)

    @classmethod
    def from_random_crystal(self, rc):
        """
        Initialize from a random_crystal object
        """
        if rc.valid:
            return Xstruct(rc.wyckoff_sites, rc.lattice, rc.group)
        else:
            return None

    def symmetrized_coords(self, coords):
        start_index = 0
        new_coords = []
        for ws in self.struc.wyckoff_sites:
            # Get coordinates associated with WP
            original_coords = coords[start_index : start_index + ws.multiplicity]
            # Re-generate the points from the first generating point
            gen_coord = coords[start_index]
            wp_coords0 = apply_ops(gen_coord, ws.wp.generators_m)
            # Calculate the PBC translation applied to each point
            translations = original_coords - wp_coords0
            frac_translations = np.dot(
                translations, np.linalg.inv(self.struc.lattice_matrix)
            )
            frac_translations = np.round(frac_translations)
            cart_translations = np.dot(frac_translations, self.struc.lattice_matrix)
            # Subtract translations from original coords
            translated_coords = original_coords - cart_translations
            # Apply inverse WP operations to get generating_points
            inverse_ops = ws.wp.inverse_generators_m
            generating_points = apply_ops_diagonal(translated_coords, inverse_ops)
            # Find the average of the generating points
            average_point = np.sum(generating_points, axis=0) / ws.multiplicity
            # Convert to fractional coordinates
            average_point = np.dot(
                average_point, np.linalg.inv(self.struc.lattice_matrix)
            )
            # Apply the WP operations to the averaged generating point
            wp_coords = apply_ops(average_point, ws.wp)
            # Convert back to Cartesian coordintes
            wp_coords = np.dot(wp_coords, self.struc.lattice_matrix)
            # Re-apply the cartesian translations
            wp_coords = wp_coords + cart_translations
            if len(new_coords) == 0:
                new_coords = wp_coords
            else:
                new_coords = np.vstack([new_coords, wp_coords])
            start_index += ws.multiplicity
        return new_coords

    def symmetrized_force(self, coords):
        start_index = 0
        new_coords = []
        for ws in self.struc.wyckoff_sites:
            # Get coordinates associated with WP
            original_coords = coords[start_index : start_index + ws.multiplicity]
            # Re-generate the points from the first generating point
            gen_coord = coords[start_index]
            wp_coords0 = apply_ops(gen_coord, ws.wp.generators_m)
            # Apply inverse WP operations to get generating_points
            inverse_ops = ws.wp.inverse_generators_m
            generating_points = apply_ops_diagonal(wp_coords0, inverse_ops)
            # Find the average of the generating points
            average_point = np.sum(generating_points, axis=0) / ws.multiplicity
            # Convert to fractional coordinates
            average_point = np.dot(
                average_point, np.linalg.inv(self.struc.lattice_matrix)
            )
            # For forces, we do not apply translational WP operations, so we truncate the ops
            matrices = np.array([op.affine_matrix[:3, :3] for op in ws.wp])
            # Apply the truncated WP operations to the averaged generating point
            wp_coords = np.dot(average_point, matrices)
            # Convert back to Cartesian coordintes
            wp_coords = np.dot(wp_coords, self.struc.lattice_matrix)
            if len(new_coords) == 0:
                new_coords = wp_coords
            else:
                new_coords = np.vstack([new_coords, wp_coords])
            start_index += ws.multiplicity
        return new_coords

    def symmetrized_stress(self, stress):
        # Make the matrix lower-diagonal
        for (i, j) in [(1, 0), (2, 1), (2, 2)]:
            stress[i][j] = stress[i][j] + stress[j][i]
            stress[j][i] = 0
        # Normalize the matrix
        snm = self.struc.lattice.stress_normalization_matrix
        m2 = np.multiply(stress, snm)
        # Normalize the on-diagonal elements
        indices = self.struc.lattice.stress_indices
        if len(indices) == 2:
            total = 0
            for index in indices:
                total += stress[index]
            value = total ** 0.5
            for index in inices:
                m2[index] = value
        elif len(indices) == 3:
            total = 0
            for index in indices:
                total += stress[index]
            value = np.cbrt(total)
            for index in inices:
                m2[index] = value
        return m2


class Mstruct(Structure):
    def __init__(self, mol_sites, lattice, group):
        self.mol_sites = mol_sites
        self.group = group
        # Extract coords and species
        species = []
        coords = []
        for ms in mol_sites:
            coords1, species1 = ms.get_coords_and_species()
            species += species1
            for c in coords1:
                coords.append(c)
        coords = np.array(coords)
        Structure.__init__(self, lattice.matrix, species, coords)

    @classmethod
    def from_molecular_crystal(self, mc):
        """
        Initialize from a molecular_crystal object
        """
        if mc.valid:
            return Mstruct(mc.mol_generators, mc.lattice, mc.group)
        else:
            return None


if __name__ == "__main__":
    from pyxtal.crystal import random_crystal

    c = random_crystal(225, ["C"], [4], 1)
    x = Xstruct.from_random_crystal(c)
    print(x.group)

    from pyxtal.molecular_crystal import molecular_crystal

    m = molecular_crystal(20, ["H2O"], [8], 1)
    x = Mstruct.from_molecular_crystal(m)
    print(x.group)
