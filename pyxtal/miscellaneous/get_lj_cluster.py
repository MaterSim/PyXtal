import numpy as np
from scipy.spatial.distance import cdist
import json
from monty.serialization import MontyEncoder
import requests
from pymatgen.symmetry.analyzer import PointGroupAnalyzer
from pymatgen.core import Molecule

black_list = {
    "6": "Oh",  # takes long time
    "13": "Ih",  # takes long time
    "86": "C1",
    "92": "C3v",  # pymatgen returns Td
    "94": "C1",
    "95": "C1",
    "96": "C1",
    "97": "C1",
    "127": "C2v",  # pymatgen returns C124v
}


def parse_symmetry(pos):
    mol = Molecule(["C"] * len(pos), pos)
    pga = PointGroupAnalyzer(mol)
    return pga.sch_symbol


def get_lj_from_url(N, address="http://doye.chem.ox.ac.uk/jon/structures/LJ/points/"):
    url_address = address + str(N)
    data_str = requests.get(url_address).text
    pos = parse_url_text(data_str)
    energy = LJ(pos)
    pos = np.reshape(pos, (N, 3))
    if str(N) in black_list.keys():
        sym = black_list[str(N)]
    else:
        sym = parse_symmetry(pos)

    cluster = {
        "name": N,
        "position": pos,
        "energy": energy,
        "url": url_address,
        "pointgroup": sym,
    }
    return cluster


def parse_url_text(data_str):
    x_array = []
    text = data_str.split("\n")
    for line in text:
        [x_array.append(float(i)) for i in line.split()]
    return np.array(x_array)


def LJ(pos):
    """
    Calculate the total energy
    input:
    pos: N*3 array which represents the atomic positions
    output
    E: the total energy
    """
    N_atom = int(len(pos) / 3)
    pos = np.reshape(pos, (N_atom, 3))
    distance = cdist(pos, pos, "euclidean")
    ids = np.triu_indices(N_atom)
    distance = distance[ids]
    distance = distance[distance > 0]
    r6 = np.power(distance, 6)
    r12 = np.multiply(r6, r6)
    return np.sum(4 * (1 / r12 - 1 / r6))


"""
Read the LJ clusters from the cambridge database
"""
clusters = []
for i in range(3, 151):
    cluster = get_lj_from_url(i)
    clusters.append(cluster)
    print(i, cluster["energy"], cluster["pointgroup"])

dumped = json.dumps(clusters, cls=MontyEncoder, indent=2)
with open("clusters.json", "w") as f:
    f.write(dumped)
