import torch

# --- Optional fix for PyTorch weights_only change (harmless if unused) ---
from torch.serialization import add_safe_globals
add_safe_globals([slice])
# ------------------------------------------------------------------------

from ase.io import read, write
from fairchem.core import pretrained_mlip, FAIRChemCalculator
import signal
from time import time
import numpy as np
from ase.constraints import FixSymmetry
from ase.filters import UnitCellFilter
from ase.optimize.fire import FIRE
import logging
import os
from ase.atoms import Atoms

# If you ever re-enable MACE/ANI you can restore imports here.
# For now, we only care about FAIRChem.
_cached_mace_mp = None  # keep this so get_calculator won't crash if 'MACE' is ever passed.
_cached_uma = None

def get_calculator(calculator):
    """
    Return an ASE calculator instance.

    Supported strings:
      - 'FAIRChem'
      - (optionally) 'ANI', 'MACE', 'MACEOFF' if you re-enable those blocks.
    Or you can pass an ASE calculator instance directly.
    """
    global _cached_mace_mp, _cached_uma

    if isinstance(calculator, str):
        if calculator == "UMA":
            if _cached_uma is None:
                predictor = pretrained_mlip.get_predict_unit("uma-s-1p1", device="cpu")
                _cached_uma = FAIRChemCalculator(predictor, task_name="omc")
            calc = _cached_uma

        elif calculator == "ANI":
            import torchani
            calc = torchani.models.ANI2x().ase()

        elif calculator == "MACE":
            from mace.calculators import mace_mp
            if _cached_mace_mp is None:
                _cached_mace_mp = mace_mp(
                    model="small",
                    dispersion=True,
                    device="cpu",
                )
            calc = _cached_mace_mp

        elif calculator == "MACEOFF":
            from mace.calculators import mace_off
            calc = mace_off(model="medium", device="cpu")

        else:
            raise ValueError(f"Unknown calculator: {calculator}")
    else:
        # already an ASE calculator instance
        calc = calculator

    return calc


def ASE_relax(
    struc,
    calculator="MACE",
    opt_lat=True,
    step=300,
    fmax=0.5,
    logfile=None,
    max_time=15.0,
    label="ase",
):
    """
    ASE optimizer used by pyxtal/DFS.

    Args:
        struc: ASE Atoms object or a pyxtal object (with .to_ase()).
        calculator: string ('FAIRChem', 'MACE', 'ANI', 'MACEOFF') or ASE calculator.
        opt_lat (bool): optimize lattice (cell) or not.
        step (int): maximum FIRE steps.
        fmax (float): force convergence criterion (eV/Å).
        logfile (str or None): FIRE log file.
        max_time (float): wall time limit in minutes.
        label (str): label for logging.

    Returns:
        ASE Atoms object if successful, otherwise None.
    """

    def handler(signum, frame):
        raise TimeoutError("Optimization timed out")

    logger = logging.getLogger()
    timeout = int(max_time * 60)  # seconds

    # Start wall-clock timeout
    signal.signal(signal.SIGALRM, handler)
    signal.alarm(timeout)

    # Convert to ASE Atoms if needed
    if isinstance(struc, Atoms):
        atoms = struc
    elif hasattr(struc, "to_ase"):
        atoms = struc.to_ase(resort=False)
    else:
        raise TypeError("ASE_relax expects an ASE Atoms or a pyxtal object with .to_ase().")

    step_init = min([30, int(step / 2)])
    _fmax = 1e5

    try:
        calc = get_calculator(calculator)#; print("The current Path:", os.getcwd())
        atoms.set_calculator(calc)
        atoms.set_constraint(FixSymmetry(atoms))

        if opt_lat:
            ecf = UnitCellFilter(atoms)
            dyn = FIRE(ecf, a=0.1, logfile=logfile) if logfile is not None else FIRE(ecf, a=0.1)
        else:
            dyn = FIRE(atoms, a=0.1, logfile=logfile) if logfile is not None else FIRE(atoms, a=0.1)

        # First stage
        dyn.run(fmax=fmax, steps=step_init)
        forces = atoms.get_forces()
        _fmax = np.sqrt((forces**2).sum(axis=1).max())

        # Second stage (only if not completely insane)
        if _fmax < 1e3 and step > step_init:
            dyn.run(fmax=fmax, steps=step - step_init)
            forces = atoms.get_forces()
            _fmax = np.sqrt((forces**2).sum(axis=1).max())
            if _fmax > 100:
                atoms = None
        else:
            atoms = None

    except TimeoutError:
        logger.warning(f"Warning {label} timed out after {timeout} seconds.")
        atoms = None

    except TypeError:
        logger.warning(f"Warning {label} spglib error in getting the lattice")
        atoms = None

    finally:
        # Always cancel alarm to avoid affecting later code
        signal.alarm(0)

    return atoms