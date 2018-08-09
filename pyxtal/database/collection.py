from pymatgen.core.structure import Molecule
import json
import os.path as op

class Collection:
    """Collection of molecular data.
    
    Example of use:

    >>> from pyxtal.database.collection import Collection
    >>> test=Collection('molecules')
    >>> test['H2O']
    Molecule Summary
    Site: O (0.0000, 0.0000, 0.0000)
    Site: H (0.2774, 0.8929, 0.2544)
    Site: H (0.6068, -0.2383, -0.7169)
    >>> list(test)
    ['H2O', 'CH4']
    """

    def __init__(self, name='molecules'):
        """Create a collection lazily.

        Will read data from json file when needed.

        A collection can be iterated over to get the Atoms objects and indexed
        with names to get individual members.

        Attributes:

        name: str
            Name of collection.
        data: dict
            Data dictionary.
        filename: str
            Location of json file.
        """

        self.name = name
        self._data = {}
        self.filename = op.join(op.dirname(__file__), name + '.json')
        with open(self.filename,"r") as f:
            self.content = json.load(f)

    def __getitem__(self, name):
        self._read(name)
        return self._data

    def __iter__(self):
        for dct in self.content:
            yield dct['name']

    def _read(self, name):
        """read the data by name and convert it to pymatgen format"""

        for dct in self.content:
            if dct['name'].lower() == name.lower():
                pos = dct['xyz']
                symbols = dct['elements']
                self._data = Molecule(symbols, pos)
                
