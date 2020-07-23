from structure import *


def so(string):
    return SymmOp.from_xyz_string(string)


def extras(sg):
    """
    Returns the additional symmetry sites ('+1/2,+1/2,+1/2', etc.) for a given spacegroup
    """
    symbol = sg_symbol_from_int_number(sg)
    letter = symbol[0]
    if letter == "P":
        return [so("x, y, z")]
    elif letter == "A":
        return [so("x, y, z"), so("x, y+1/2, z+1/2")]
    elif letter == "C":
        return [so("x, y, z"), so("x+1/2, y+1/2, z")]
    elif letter == "I":
        return [so("x, y, z"), so("x+1/2, y+1/2, z+1/2")]
    elif letter in ["R"]:
        return [so("x, y, z"), so("x+2/3, y+1/3, z+1/3"), so("x+1/3, y+2/3, z+2/3")]
    elif letter in ["F"]:
        return [
            so("x, y, z"),
            so("x, y+1/2, z+1/2"),
            so("x+1/2, y, z+1/2"),
            so("x+1/2, y+1/2, z"),
        ]
    else:
        return "Error: Could not determine lattice type"


new_wyckoffs_all = [None]
for sg in range(1, 231):
    new_wyckoffs_all.append([])
    ext = extras(sg)
    old_wyckoffs = get_wyckoffs(sg)
    for wp in old_wyckoffs:
        new_wyckoffs_all[-1].append([])
        short_length = int(len(wp) / len(ext))
        short_wp = wp[:short_length]
        for e in ext:
            for point in short_wp:
                new_wyckoffs_all[-1][-1].append((e * point).as_xyz_string())

