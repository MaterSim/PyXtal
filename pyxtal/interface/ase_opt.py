# --- Temporary fix for PyTorch 2.6 weights_only change ---
import torch as _torch
import functools as _functools
_original_torch_load = _torch.load
@_functools.wraps(_original_torch_load)
def _patched_torch_load(*args, **kwargs):
    kwargs.setdefault('weights_only', False)
    return _original_torch_load(*args, **kwargs)
_torch.load = _patched_torch_load
# ----------------------------------------------------------

import signal
import numpy as np
import os
import contextlib
from ase.constraints import FixSymmetry
from ase.filters import UnitCellFilter
from ase.optimize.fire import FIRE
import logging
from ase.atoms import Atoms

_calc_cache = {}

DEFAULT_MODELS = {
    "MACE": "small",
    "MACEOFF": "medium",
    "ORB": "conservative-inf-omat",
    "ORB-V3": "conservative-inf-omat",
    "ORBV3": "conservative-inf-omat",
}

QUICK_MODELS = {
    "MACE": "small",
    "MACEOFF": "small",
    "ORB": "direct-20-omat",
    "ORB-V3": "direct-20-omat",
    "ORBV3": "direct-20-omat",
}

@contextlib.contextmanager
def _suppress_stdout_stderr():
    """Redirect stdout and stderr to /dev/null."""
    with open(os.devnull, 'w') as devnull:
        with contextlib.redirect_stdout(devnull), contextlib.redirect_stderr(devnull):
            yield

def resolve_model(calculator, model=None, quick=False):
    """Return the model name to use for a calculator."""
    if quick:
        return QUICK_MODELS.get(calculator, DEFAULT_MODELS.get(calculator))
    if model is None:
        return DEFAULT_MODELS.get(calculator)
    return model


def get_calculator(calculator, model=None, quick=False):
    """
    Return an ASE calculator instance.

    Supported strings:
      - 'UMA', 'ANI', 'MACE', 'MACEOFF', 'ORB' / 'ORB-V3' / 'ORBV3'
    Or you can pass an ASE calculator instance directly.

    model: optional model variant (e.g. MACE 'small'/'medium'/'large',
           ORB 'conservative-inf-omat'/'direct-20-omat'/'direct-inf-omat').
    quick: use a cheaper/faster model preset for testing.
    """
    if not isinstance(calculator, str):
        return calculator

    model = resolve_model(calculator, model=model, quick=quick)
    cache_key = (calculator, model)
    if cache_key in _calc_cache:
        return _calc_cache[cache_key]

    with _suppress_stdout_stderr():
        if calculator == "UMA":
            from fairchem.core import pretrained_mlip, FAIRChemCalculator
            predictor = pretrained_mlip.get_predict_unit("uma-s-1p1")
            calc = FAIRChemCalculator(predictor, task_name="omc")

        elif calculator == "ANI":
            import torchani
            calc = torchani.models.ANI2x().ase()

        elif calculator == "MACE":
            from mace.calculators import mace_mp
            calc = mace_mp(model=model, dispersion=True)

        elif calculator == "MACEOFF":
            from mace.calculators import mace_off
            calc = mace_off(model=model)

        elif calculator in ("ORB-V3", "ORBV3", "ORB"):
            try:
                from orb_models.forcefield import pretrained
                from orb_models.forcefield.calculator import ORBCalculator
            except ImportError as e:
                raise ImportError(
                    "Please install ORB-V3 from https://github.com/orbital-materials/orb-models "
                    "or pip install orb-models to use the 'ORB-V3' calculator. Error: " + str(e)
                )
            orb_loaders = {
                "conservative-inf-omat": pretrained.orb_v3_conservative_inf_omat,
                "direct-20-omat": pretrained.orb_v3_direct_20_omat,
                "direct-inf-omat": pretrained.orb_v3_direct_inf_omat,
            }
            if model not in orb_loaders:
                raise ValueError(
                    f"Unknown ORB model '{model}'. Choose from: {', '.join(orb_loaders)}"
                )
            device = "cpu"
            # ORB defaults to torch.compile on CPU, which is slow and can fail on macOS.
            orbff = orb_loaders[model](
                device=device,
                precision="float32-high",
                compile=False,
            )
            orb_calc_kwargs = {"device": device}
            if model.startswith("direct"):
                orb_calc_kwargs["conservative"] = False
            if model == "direct-20-omat":
                orb_calc_kwargs["max_num_neighbors"] = 20
            calc = ORBCalculator(orbff, **orb_calc_kwargs)

        else:
            raise ValueError(f"Unknown calculator: {calculator}")

    _calc_cache[cache_key] = calc
    return calc


def ASE_relax(
    struc,
    calculator="MACE",
    model=None,
    quick=False,
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
        atoms.calc = get_calculator(calculator, model=model, quick=quick)
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
