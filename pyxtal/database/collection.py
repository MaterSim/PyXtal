import json
import os.path as op

from pymatgen.core.structure import Molecule


class Collection:
    """Collection of molecular data.
    Used for obtaining pymatgen objects from a small database file.

    Example of use:

    >>> from pyxtal.database.collection import Collection
    >>> test=Collection('molecules')
    >>> test['H2O']
    Molecule Summary
    Site: O (0.0000, 0.0000, 0.0000)
    Site: H (0.2774, 0.8929, 0.2544)
    Site: H (0.6068, -0.2383, -0.7169)
    >>> list(test)
    ['C60', 'H2O', 'CH4', 'NH3', 'benzene', 'naphthalene', 'anthracene', 'tetracene', 'pentacene', 'coumarin', 'resorcinol', 'benzamide', 'aspirin', 'ddt', 'lindane', 'glycine', 'glucose', 'ROY']

    Args:
        name: the type of collection to get. Defaults to "molecules"
    """

    def __init__(self, name="molecules"):
        """Create a collection lazily.

        Will read data from json file when needed.

        A collection can be iterated over to get the Atoms objects and indexed
        with names to get individual members.

        Attributes:

        name: str
            Name of collection.
        data: object
            Pymetgen molecule object
        filename: str
            Location of json file.
        """

        self.name = name
        self._data = {}
        self.filename = op.join(op.dirname(__file__), name + ".json")
        with open(self.filename) as f:
            self.content = json.load(f)

    def __getitem__(self, name):
        self._read(name)
        if len(self._data) == 0:
            names = ""
            for dct in self.content:
                names += dct["name"] + ", "
            msg = name + " is not supported\n"
            msg += "Available molecules are:\n"
            msg += names
            raise NameError(msg)
        else:
            return self._data

    def __iter__(self):
        for dct in self.content:
            yield dct["name"]

    def _read(self, name):
        if self.name == "molecules":
            """
            read the data by name and convert it to pymatgen format
            """
            for dct in self.content:
                if dct["name"].lower() == name.lower():
                    pos = dct["xyz"]
                    symbols = dct["elements"]
                    self._data = Molecule(symbols, pos)
        elif self.name == "clusters":
            for dct in self.content:
                if dct["name"] == int(name):
                    self._data = dct

    def show_names(self):
        names = []
        for dct in self.content:
            names.append(dct["name"])
        print(names)
