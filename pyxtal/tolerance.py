import os

import numpy as np

from pyxtal.database.element import Element


class Tol_matrix:
    """
    Class for variable distance tolerance checking. Used within random_crystal
    and molecular_crystal to verify whether atoms are too close. Stores a matrix
    of atom-atom pair tolerances. Note that the matrix's indices correspond to
    atomic numbers, with the 0th entries being 0 (there is no atomic number 0).

    Args:
        prototype: a string representing the type of radii to use
            (`atomic`, `molecular`, `vdW` or `metallic`)
        factor: a float to scale the distances by.
        tuples: a list or tuple of tuples, which define custom tolerance values.
            Each tuple should be of the form (specie1, specie2, value), where
            value is the tolerance in Angstroms, and specie1 and specie2 can be
            strings, integers, Element objects, or pymatgen Specie objects.
            Custom values may also be set using set_tol
    """

    def __init__(self, *tuples, prototype="atomic", factor=1.0):
        f = factor
        self.prototype = prototype
        if prototype == "atomic":
            f *= 0.5
            attrindex = 5
            self.radius_type = "covalent"
        elif prototype == "molecular":
            attrindex = 5
            self.radius_type = "covalent"
            f *= 1.2
        elif prototype == "metallic":
            attrindex = 7
            self.radius_type = "metallic"
            f *= 0.5
        elif prototype == "vdW":
            attrindex = 6
            self.radius_type = "vdW"
        else:
            self.radius_type = "N/A"

        self.f = f
        H = Element("H")
        m = [[0.0] * (len(H.elements_list) + 1)]
        for tup1 in H.elements_list:
            m.append([0.0])
            for tup2 in H.elements_list:
                # Get the appropriate atomic radii
                val1 = tup1[5] if tup1[attrindex] is None else tup1[attrindex]
                # FIXME this is suspicious but matches exactly the original code
                val2 = (tup1[5] if tup2[5] is not None else None) if tup2[attrindex] is None else tup2[attrindex]

                if val1 is not None and val2 is not None:
                    m[-1].append(f * (val1 + val2))
                else:
                    # If no radius is found for either atom, set tolerance to None
                    m[-1].append(None)
        self.matrix = np.array(m)  # A symmetric matrix between specie pairs
        self.custom_values = []  # A list of tuples storing customized pair tolerances

        try:
            for tup in tuples:
                self.set_tol(*tup)
        except Exception as err:
            msg = "Error: Cannot not set custom tolerance value(s).\n"
            msg += "All entries should be entered using the following form:\n"
            msg += "(specie1, specie2, value), where the value is in Angstrom."
            raise RuntimeError(msg) from err

        self.radius_list = []
        for i in range(len(self.matrix)):
            if i == 0:
                continue
            x = self.get_tol(i, i)
            self.radius_list.append(x)

    def get_tol(self, specie1, specie2):
        """
        Returns the tolerance between two species.

        Args:
            specie1/2: atomic number (int or float), name (str), symbol (str),
                an Element object, or a pymatgen Specie object

        Returns:
            the tolerance between the provided pair of atomic species
        """
        if self.prototype == "single value":
            return self.matrix[0][0]
        index1 = Element.number_from_specie(specie1)
        index2 = Element.number_from_specie(specie2)
        if index1 is not None and index2 is not None:
            return self.matrix[index1][index2]
        else:
            return None

    def set_tol(self, specie1, specie2, value):
        """
        Sets the distance tolerance between two species.

        Args:
            specie1/2: atomic number (int or float), name (str), symbol (str),
                an Element object, or a pymatgen Specie object
            value:
                the tolerance (in Angstroms) to set to
        """
        index1 = Element.number_from_specie(specie1)
        index2 = Element.number_from_specie(specie2)
        if index1 is None or index2 is None:
            return
        self.matrix[index1][index2] = float(value)
        if index1 != index2:
            self.matrix[index2][index1] = float(value)
        if (index1, index2) not in self.custom_values and (
            index2,
            index1,
        ) not in self.custom_values:
            larger = max(index1, index2)
            smaller = min(index1, index2)
            self.custom_values.append((smaller, larger))

    @classmethod
    def from_matrix(self, matrix, prototype="atomic", factor=1.0, begin_with=0):
        """
        Given a tolerance matrix, returns a Tol_matrix object. Matrix indices
        correspond to the atomic number (with 0 pointing to Hydrogen by default).
        For atoms with atomic numbers not included in the matrix, the default
        value (specified by prototype) will be used, up to element 96. Note that
        if the matrix is asymmetric, only the value below the diagonal will be used.

        Args:
            matrix: a 2D matrix or list of tolerances between atomic species pairs. The
                indices correspond to atomic species (see begin_with variable description)
            prototype: a string representing the type of radii to use
                ("atomic", "molecular")
            factor: a float to scale the distances by. A smaller value means a smaller
                tolerance for distance checking
            begin_with: the index which points to Hydrogen within the matrix. Default 0

        Returns:
            a Tol_matrix object
        """
        np.array(matrix)
        tups = []
        for i, row in enumerate(matrix):
            for j, _value in enumerate(row):
                if j > i:
                    continue
                tups.append((i + 1 - begin_with, j + 1 - begin_with, matrix[i][j]))
        return Tol_matrix(*tups, prototype=prototype, factor=factor)

    @classmethod
    def from_radii(self, radius_list, prototype="atomic", factor=1.0, begin_with=0):
        """
        Given a list of atomic radii, returns a Tol_matrix object. For atom-atom pairs, uses
        the average radii of the two species as the tolerance value. For atoms with atomic
        numbers not in the radius list, the default value (specified by prototype) will be
        used, up to element 96.

        Args:
            radius_list: a list of atomic radii (in Angstroms), beginning with Hydrogen
            prototype: a string representing the type of radii to use
                ("atomic", "molecular")
            factor: a float to scale the distances by. A smaller value means a smaller
                tolerance for distance checking
            begin_with: the index which points to Hydrogen within the list. Default 0

        Returns:
            a Tol_matrix object
        """
        tups = []
        f = factor * 0.5
        for i, r1 in enumerate(radius_list):
            for j, r2 in enumerate(radius_list):
                if j > i:
                    continue
                tups.append((i + 1 - begin_with, j + 1 - begin_with, f * (r1 + r2)))
        return Tol_matrix(*tups, prototype=prototype, factor=factor)

    @classmethod
    def from_single_value(self, value):
        """
        Creates a Tol_matrix which only has a single tolerance value. Using
        `get_tol` will always return the same value.

        Args:
            value: the tolerance value to use

        Returns:
            a Tol_matrix object
        """
        tm = Tol_matrix()
        tm.prototype = "single value"
        tm.matrix = np.array([[value]])
        tm.custom_values = [(1, 1)]
        tm.radius_type = "N/A"
        return tm

    def __getitem__(self, index):
        Element.number_from_specie(index)
        return self.matrix[index]

    def __str__(self):
        s = "\n--Tol_matrix class object--"
        s += "\nPrototype: " + str(self.prototype)
        s += "\nAtomic radius type: " + str(self.radius_type)
        s += "\nRadius scaling factor: " + str(self.f)
        if self.prototype == "single value":
            s += "\nCustom tolerance value: " + str(self.matrix([0][0]))
        else:
            if self.custom_values == []:
                s += "\nCustom tolerance values: None"
            else:
                s += "\nCustom tolerance values:"
                for tup in self.custom_values:
                    name1 = str(Element(tup[0]).short_name)
                    name2 = str(Element(tup[1]).short_name)
                    s += f"\n{name1:s}-{name2:s}: {self.get_tol(tup[0], tup[1]):6.3f}"
        return s

    def to_file(self, filename=None):
        """
        Creates a file with the given filename.

        Args:
            filename: the file path

        Returns:
            Nothing. Creates a file at the specified path
        """
        if filename is None:
            filename = "custom_tol_matrix"
        # Check if filename already exists
        # If it does, add a new number to end of filename
        if os.path.exists(filename):
            i = 1
            while True:
                outdir = filename + "_" + str(i)
                if not os.path.exists(outdir):
                    break
                i += 1
                if i > 10000:
                    return "Cannot create file: too many files already created."
        else:
            outdir = filename
        try:
            np.save(filename, [self])
            return "Output file to " + outdir + ".npy"
        except:
            return "Error: Could not save Tol_matrix to file."

    @classmethod
    def from_file(self, filename):
        try:
            tm = np.load(filename)[0]
            if type(tm) == Tol_matrix:
                return tm
            else:
                raise RuntimeError("invalid file for Tol_matrix: ", filename)
        except Exception as err:
            raise RuntimeError("Could not load Tol_matrix from file: ", filename) from err


if __name__ == "__main__":
    from pyxtal.molecule import pyxtal_molecule

    for p in ["atomic", "molecular", "vdW", "metallic"]:
        tm = Tol_matrix(prototype=p)
        print(p, tm.get_tol("C", "H"))
        m = pyxtal_molecule("aspirin", tm=tm)
        print(m.tols_matrix)
