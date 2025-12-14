# --- Optional fix for PyTorch weights_only change (harmless if unused) ---
from torch.serialization import add_safe_globals
add_safe_globals([slice])
# ------------------------------------------------------------------------

import signal
import numpy as np
from ase.constraints import FixSymmetry
from ase.filters import UnitCellFilter
from ase.optimize.fire import FIRE
import logging
from ase.atoms import Atoms

_cached_mace = None
_cached_uma = None
_cached_ani = None

def get_calculator(calculator):
    """
    Return an ASE calculator instance.

    Supported strings:
      - 'FAIRChem', 'ANI', 'MACE', 'MACEOFF' if you re-enable those blocks.
    Or you can pass an ASE calculator instance directly.
    """
    global _cached_mace, _cached_uma, _cached_ani

    if isinstance(calculator, str):
        if calculator == "UMA":
            if _cached_uma is None:
                from fairchem.core import pretrained_mlip, FAIRChemCalculator
                predictor = pretrained_mlip.get_predict_unit("uma-s-1p1")
                _cached_uma = FAIRChemCalculator(predictor,
                                                 task_name="omc")
            calc = _cached_uma

        elif calculator == "ANI":
            if _cached_ani is None:
                import torchani
                _cached_ani = torchani.models.ANI2x().ase()
            calc = _cached_ani

        elif calculator == "MACE":
            if _cached_mace is None:
                from mace.calculators import mace_mp
                _cached_mace = mace_mp(model="small", dispersion=True)
            calc = _cached_mace

        elif calculator == "MACEOFF":
            if _cached_mace is None:
                from mace.calculators import mace_off
                _cached_mace = mace_off(model="medium")#, device="cpu")
            calc = _cached_mace

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
        atoms.calc = get_calculator(calculator)#; print("The current Path:", os.getcwd())
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
