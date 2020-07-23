from pyxtal.constants import pyxtal_verbosity
from warnings import warn


def printx(text, priority=1):
    """
    Custom printing function based on verbosity.

    Args:
        text: string to be passed to print
        priority: the importance of printing the message
            0: Critical; must be printed
            1: Warning; unexpected error but program functioning
            2: Info; useful but not necessary print out
            3: Debug; detailed information for debugging

    Returns:
        Nothing
    """
    if priority <= 1:
        warn(text)
        return
    else:
        if priority <= pyxtal_verbosity:
            print(text)
