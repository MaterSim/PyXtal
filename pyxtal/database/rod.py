"""
Module for handling Rod groups for generation of 1d crystals. A Rod group is
the combination of a crystallographic point group with 1 direction of periodic
translational symmetry. There are 75 Rod groups, each of which corresponds to
a 3d space group, with a possible permutation of the axes.
"""

class rodgroup:
    """
    Class for storing Rod groups. Used for the generation of 1d crystals.
    Each layer group corresponds to the point group of a 3d space group. Thus, we
    can access the symmetry information using code for standard space groups. The
    primary class, rodgroup, allows access to the Rod group number, space group
    number, and space group symbol.

    Args:
        input_value: The Rod group number or symbol
    """
    def __init__(self, input_value):

        # list of layer group number, symbol, space group number and permutation
        #TODO: Add Rod group information
        self.group_list = [
            (1,     'P1',       1,      [1,2,3,1,2,3,3]),
            (2,     'P-1',      2,      [1,2,3,1,2,3,3])]

        self.input = str(input_value)
        self.error = False
        if self.input.isdigit():
            self.lg = int(self.input)
            """The Rod group number (between 1 and 75)"""
        else:
            for i, e1 in enumerate(self.group_list):
                if e1[1] == self.input:
                    self.lg = e1[0]
                    break
        if (self.lg is not None) and (0<self.lg<76):
            self.symbol=self.group_list[self.lg-1][1]
            """The Hermann-Mauguin symbol of the Rod group"""
            self.sgnumber = self.group_list[self.lg-1][2]
            """The international space group number (between 1 and 230)"""
            self.permutation=self.group_list[self.lg-1][3]
        else:
            self.error = True
            print('Error:   unable to find the Rod group, check your input: ', self.input)

    def print_all(self):
        if self.error is False:
            print('Rod group number: ', self.lg)
            print('Rod group symbol: ', self.symbol)
            print('Space group number: ', self.sgnumber)
            print('Permutation:        ', self.permutation)
