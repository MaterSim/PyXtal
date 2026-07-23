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
try:
    import fcntl
except ImportError:
    fcntl = None
from ase.constraints import FixSymmetry
from ase.filters import UnitCellFilter
from ase.optimize.fire import FIRE
from ase.atoms import Atoms

_calc_cache = {}

DEFAULT_MODELS = {
    "MACE": "medium",
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


@contextlib.contextmanager
def _model_load_lock():
    """Serialize MACE model downloads across worker processes."""
    if fcntl is None:
        yield
        return
    lock_dir = os.path.expanduser("~/.cache/mace")
    os.makedirs(lock_dir, exist_ok=True)
    lock_path = os.path.join(lock_dir, ".model_load.lock")
    with open(lock_path, "a") as lock_file:
        fcntl.flock(lock_file.fileno(), fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(lock_file.fileno(), fcntl.LOCK_UN)


def _clear_mace_model_cache():
    cache_dir = os.path.expanduser("~/.cache/mace")
    if not os.path.isdir(cache_dir):
        return
    for name in os.listdir(cache_dir):
        if name.endswith(".model"):
            os.remove(os.path.join(cache_dir, name))


def _is_corrupt_model_error(exc):
    if isinstance(exc, EOFError):
        return True
    if isinstance(exc, RuntimeError):
        msg = str(exc)
        return "PytorchStreamReader" in msg or "failed finding central directory" in msg
    return False


_MACEOFF_CACHE = {
    "small": "MACE-OFF23_small.model",
    "medium": "MACE-OFF23_medium.model",
    "large": "MACE-OFF23_large.model",
}

_MACE_CACHE = {
    "small": "46jrkm3v",
    "medium": "5yyxdm76",
    "large": "5f5yavf3",
}


def _mace_cache_path(calculator, model):
    cache_dir = os.path.expanduser("~/.cache/mace")
    if calculator == "MACEOFF":
        fname = _MACEOFF_CACHE.get(model, _MACEOFF_CACHE["medium"])
    elif calculator == "MACE":
        fname = _MACE_CACHE.get(model, _MACE_CACHE["medium"])
    else:
        return None
    return os.path.join(cache_dir, fname)


def _mace_cache_ready(calculator, model, min_bytes=1_000_000):
    path = _mace_cache_path(calculator, model)
    return bool(path and os.path.isfile(path) and os.path.getsize(path) >= min_bytes)


def _download_mace_model(calculator, model):
    """Download MACE/MACEOFF weights to disk (caller must hold the download lock)."""
    with _suppress_stdout_stderr():
        if calculator == "MACEOFF":
            from mace.calculators import mace_off
            mace_off(model=model)
        elif calculator == "MACE":
            from mace.calculators import mace_mp
            mace_mp(model=model, dispersion=True)
        else:
            raise ValueError(f"Unknown calculator: {calculator}")


def ensure_mace_model_cached(calculator, model=None, quick=False, min_bytes=1_000_000):
    """Ensure MACE/MACEOFF weights exist on disk without keeping a calculator loaded."""
    if calculator not in ("MACE", "MACEOFF"):
        return None
    model = resolve_model(calculator, model=model, quick=quick)
    if _mace_cache_ready(calculator, model, min_bytes=min_bytes):
        return _mace_cache_path(calculator, model)
    with _model_load_lock():
        if _mace_cache_ready(calculator, model, min_bytes=min_bytes):
            return _mace_cache_path(calculator, model)
        _download_mace_model(calculator, model)
    return _mace_cache_path(calculator, model)

def resolve_model(calculator, model=None, quick=False):
    """Return the model name to use for a calculator."""
    if quick:
        return QUICK_MODELS.get(calculator, DEFAULT_MODELS.get(calculator))
    if model is None:
        return DEFAULT_MODELS.get(calculator)
    return model


def _import_orb_calculator():
    try:
        from orb_models.forcefield.inference.calculator import ORBCalculator
        return ORBCalculator, 3
    except ImportError:
        pass
    try:
        from orb_models.forcefield.calculator import ORBCalculator
        return ORBCalculator, 2
    except ImportError as e:
        raise ImportError(
            "Please install ORB from https://github.com/orbital-materials/orb-models "
            "or pip install orb-models to use the ORB calculator. Error: " + str(e)
        )


def _resolve_orb_loader(pretrained, model):
    aliases = {
        "conservative-inf-omat": "orb_v3_conservative_inf_omat",
        "direct-20-omat": "orb_v3_direct_20_omat",
        "direct-inf-omat": "orb_v3_direct_inf_omat",
    }
    if hasattr(pretrained, model):
        return getattr(pretrained, model)
    alias = aliases.get(model)
    if alias and hasattr(pretrained, alias):
        return getattr(pretrained, alias)
    orb_models = getattr(pretrained, "ORB_PRETRAINED_MODELS", {})
    key = model.replace("_", "-")
    if key in orb_models:
        return orb_models[key]
    return None


def _build_orb_calculator(model):
    """Build an ORB ASE calculator, supporting orb-models v2 and v3 APIs."""
    from orb_models.forcefield import pretrained

    ORBCalculator, _orb_api = _import_orb_calculator()
    device = "cpu"

    if hasattr(pretrained, "orb_v3_conservative_inf_omat"):
        loader = _resolve_orb_loader(pretrained, model)
        if loader is None:
            known = sorted(
                set(
                    [
                        "conservative-inf-omat",
                        "direct-20-omat",
                        "direct-inf-omat",
                        *getattr(pretrained, "ORB_PRETRAINED_MODELS", {}).keys(),
                    ]
                )
            )
            raise ValueError(
                f"Unknown ORB model '{model}'. Choose from: {', '.join(known)}"
            )
        loaded = loader(
            device=device,
            precision="float32-high",
            compile=False,
        )
        if isinstance(loaded, tuple):
            orbff, atoms_adapter = loaded
            calc_kwargs = {"device": device, "atoms_adapter": atoms_adapter}
            if model == "direct-20-omat":
                calc_kwargs["max_num_neighbors"] = 20
            return ORBCalculator(orbff, **calc_kwargs)

        orb_calc_kwargs = {"device": device}
        if model.startswith("direct"):
            orb_calc_kwargs["conservative"] = False
        if model == "direct-20-omat":
            orb_calc_kwargs["max_num_neighbors"] = 20
        return ORBCalculator(loaded, **orb_calc_kwargs)

    from orb_models.forcefield.atomic_system import SystemConfig

    v3_to_v2 = {
        "conservative-inf-omat": "orb-d3-v2",
        "direct-20-omat": "orb-d3-xs-v2",
        "direct-inf-omat": "orb-v2",
    }
    orb_models = getattr(pretrained, "ORB_PRETRAINED_MODELS", {})
    loader_key = v3_to_v2.get(model, model)
    if loader_key not in orb_models:
        raise ValueError(
            f"Unknown ORB model '{model}'. Choose from: {', '.join(sorted(set(v3_to_v2) | set(orb_models)))}"
        )
    orbff = orb_models[loader_key](device=device)
    max_num_neighbors = 20 if model == "direct-20-omat" else 20
    system_config = SystemConfig(radius=10.0, max_num_neighbors=max_num_neighbors)
    return ORBCalculator(
        orbff,
        device=device,
        system_config=system_config,
        brute_force_knn=False,
    )


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

    def _build_calculator():
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
                calc = _build_orb_calculator(model)

            else:
                raise ValueError(f"Unknown calculator: {calculator}")
        return calc

    use_model_lock = calculator in ("MACE", "MACEOFF")
    last_exc = None
    for attempt in range(2):
        try:
            if use_model_lock and not _mace_cache_ready(calculator, model):
                with _model_load_lock():
                    if not _mace_cache_ready(calculator, model):
                        _download_mace_model(calculator, model)
            calc = _build_calculator()
            break
        except Exception as exc:
            last_exc = exc
            if attempt == 0 and use_model_lock and _is_corrupt_model_error(exc):
                with _model_load_lock():
                    _clear_mace_model_cache()
                continue
            raise
    else:
        raise last_exc

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
        logfile (str or None): FIRE log file. ``None`` discards FIRE step logs
            (writes to ``os.devnull``).
        max_time (float): wall time limit in minutes.
        label (str): label for logging.

    Returns:
        ASE Atoms object if successful, otherwise None.
    """

    def handler(signum, frame):
        raise TimeoutError("Optimization timed out")

    def _warn(msg):
        # Avoid logging.getLogger(): in multiprocessing workers the root
        # logger may point at a closed stream, which produces
        # "--- Logging error ---" instead of a useful message.
        try:
            print(msg, flush=True)
        except Exception:
            pass

    timeout = int(max_time * 60)  # seconds
    # Discard FIRE per-step output unless the caller passed a real path.
    fire_logfile = logfile if logfile is not None else os.devnull

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
        try:
            atoms.set_constraint(FixSymmetry(atoms))
        except (TypeError, AttributeError, RuntimeError) as exc:
            # spglib can return None for some lattices; continue without symmetry constraint
            _warn(
                f"Warning {label} FixSymmetry failed ({exc}); "
                "relaxing without symmetry constraint"
            )

        if opt_lat:
            ecf = UnitCellFilter(atoms)
            dyn = FIRE(ecf, a=0.1, logfile=fire_logfile)
        else:
            dyn = FIRE(atoms, a=0.1, logfile=fire_logfile)

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
        # Cancel before any I/O so a second alarm cannot interrupt cleanup.
        signal.alarm(0)
        _warn(f"Warning {label} timed out after {timeout} seconds.")
        atoms = None

    finally:
        # Always cancel alarm to avoid affecting later code
        signal.alarm(0)

    return atoms
