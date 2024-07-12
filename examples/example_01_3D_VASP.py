"""
This script:
1. Generates random crystal structures
2. Performs multiple steps of optimization with ASE-VASP

Requirements:
- ASE must be installed
- VASP must be callable from ASE
"""

import warnings
from time import time
from typing import TYPE_CHECKING

import numpy as np
from ase.db import connect

from pyxtal import pyxtal
from pyxtal.interface.vasp import optimize

if TYPE_CHECKING:
    from ase.atoms import Atoms

warnings.filterwarnings("ignore")

N_STRUCTURES: int = 10
ELEMENTS: dict[str, list[int]] = {"C": [2, 4]}
OPTIMIZATION_LEVELS: list[int] = [0, 2]  # Add 3 if needed
CALCULATION_DIR: str = "Calc"
DB_FILENAME: str = "C-VASP.db"

rng = np.random.default_rng()


def generate_random_crystal() -> pyxtal:
    """Generate a random crystal structure."""
    while True:
        sg = rng.integers(2, 230)
        species = []
        num_ions = []
        for ele, count in ELEMENTS.items():
            species.append(ele)
            if len(count) == 2:
                num = rng.integers(count[0], count[1])
                num_ions.append(num)
            else:
                num_ions.append(count[0])

        crystal = pyxtal()
        try:
            crystal.from_random(3, sg, species, num_ions, force_pass=True)
        except pyxtal.msg.Comp_CompatibilityError:
            continue

        if crystal.valid:
            return crystal


def optimize_and_save(crystal: pyxtal, index: int) -> None:
    """Optimize the crystal structure and save results to database."""
    struc, energy, cputime, error = optimize(crystal, path=CALCULATION_DIR, levels=OPTIMIZATION_LEVELS)

    if not error:
        ase_atoms: Atoms = struc.to_ase()
        energy_per_atom = energy / len(ase_atoms)
        cputime_minutes = cputime / 60.0

        with connect(DB_FILENAME) as db:
            kvp = {
                "spg": struc.group.symbol,
                "dft_energy": energy_per_atom,
            }
            db.write(ase_atoms, key_value_pairs=kvp)

        print(
            f"{index:3d} {crystal.group.symbol:12s} -> {struc.group.symbol:12s}"
            f" {energy_per_atom:6.3f} eV/atom {cputime_minutes:6.2f} min"
        )


def main() -> None:
    """Main function to generate and optimize crystal structures."""
    for i in range(N_STRUCTURES):
        t0 = time()
        crystal = generate_random_crystal()
        optimize_and_save(crystal, i)
        print(f"Total time for structure {i}: {time() - t0:.2f} s")


if __name__ == "__main__":
    main()
