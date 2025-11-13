import os
import subprocess
import re

import numpy as np
from ase import Atoms
from pyxtal import pyxtal
from pyxtal.lattice import Lattice


class LATTECalc:
    """
    Minimal LATTE calculator for PyXtal / ASE.

    Assumes the current working directory (workdir) contains:
      - LATTE executable (exe, e.g. './LATTE_DOUBLE')
      - TBparam/ directory with control.in etc.
      - MDcontroller
      - bl/ directory (will be created if missing)

    It does NOT change directories or create extra folders.
    """

    def __init__(
        self,
        struc,
        exe="./LATTE_DOUBLE",
        workdir=".",
        log_name="latte.log",
        relax=None,   # currently only used as a flag, not editing control.in
        timeout=3600,
        control_relax=None,    # int 0/1, if None => don’t touch RELAX
        control_maxiter=None,  # int, if None => don’t touch MAXITER
        control_rlxftol=None,  # float, if None => don’t touch RLXFTOL
        control_restart=None,  # int 0/1, if None => don’t touch RESTART
        ):
        # keep paths absolute to be safe
        self.workdir = os.path.abspath(workdir)
        self.exe = os.path.abspath(os.path.join(self.workdir, exe))
        self.log_name = log_name
        self.log_path = os.path.join(self.workdir, self.log_name)
        self.relax = relax if relax is not None else self._read_relax_from_control()
        self.timeout = timeout
        self.control_relax   = control_relax
        self.control_maxiter = control_maxiter
        self.control_rlxftol = control_rlxftol
        self.control_restart = control_restart

        # --- handle structure (pyxtal or ASE Atoms) ---
        if isinstance(struc, pyxtal):
            self.pyxtal = struc
            self.species = struc.species
            struc = struc.to_ase(resort=False)
        else:
            self.pyxtal = None

        if isinstance(struc, Atoms):
            self.lattice = Lattice.from_matrix(struc.cell)
            self.cart_coords = struc.get_positions()  # Å
            self.sites = struc.get_chemical_symbols()
        else:
            raise NotImplementedError("LATTECalc only supports pyxtal or ASE Atoms objects")

        self.structure = struc

        # results
        self.energy = None
        self.energy_per_atom = None
        self.forces = None
        self.stress = None
        self.optimized = False
        self.cputime = 0.0
        self.error = False

    def run(self):
        """
        Run LATTE once:
          1) write bl/inputblock.dat
          2) run LATTE executable
          3) parse log for total energy
        """
        self._write_inputblock()
        self._apply_control_overrides()
        self._execute()
        self._read_output()
        # optional later: self._read_restart_geometry()

    def to_ase(self):
        #Return ASE Atoms for the (possibly) relaxed structure.
        #Right now we just keep the original positions; deal later

        cell = self.lattice.matrix
        positions = self.cart_coords
        return Atoms(self.sites, positions=positions, cell=cell, pbc=True)

    def to_pyxtal(self):
        """
        Convert final structure back to pyxtal.
        """
        ase_atoms = self.to_ase()
        for tol in [1e-2, 1e-3, 1e-4, 1e-5]:
            try:
                s = pyxtal()
                s.from_seed(ase_atoms, tol=tol)
                return s
            except Exception:
                pass
        return None

    def _write_inputblock(self):
        """
        Write bl/inputblock.dat in the format LATTE expects:
          Natoms
          3 lines: lattice vectors
          Natoms lines: element x y z (Å)
        """
        bl_dir = os.path.join(self.workdir, "bl")
        if not os.path.exists(bl_dir):
            os.makedirs(bl_dir)

        natoms = len(self.sites)
        cell = self.lattice.matrix

        inputblock_path = os.path.join(bl_dir, "inputblock.dat")
        with open(inputblock_path, "w") as f:
            f.write(f"{natoms:d}\n")
            for i in range(3):
                f.write(f"{cell[i,0]: .10f} {cell[i,1]: .10f} {cell[i,2]: .10f}\n")
            for elem, pos in zip(self.sites, self.cart_coords):
                f.write(f"{elem:s} {pos[0]: .10f} {pos[1]: .10f} {pos[2]: .10f}\n")
    def _apply_control_overrides(self):
        """
        Update the RELAX line in TBparam/control.in with:
        RELAX= <int> RELAXTYPE= SD MAXITER= <int> RLXFTOL= <float>
    
        If any of these parameters are missing (None), keep the old values.
        """
        control_path = os.path.join(self.workdir, "TBparam", "control.in")
        if not os.path.exists(control_path):
            print(f"Warning: {control_path} not found; cannot apply overrides.")
            return
    
        # Read existing file
        with open(control_path, "r") as f:
            lines = f.readlines()
    
        new_lines = []
        relax_found = False
        restart_found = False
        # Defaults (will be replaced if line already has them)
        relax_val   = self.control_relax
        relax_type  = "SD"
        maxiter_val = self.control_maxiter
        rlxftol_val = self.control_rlxftol
        restart_val = self.control_restart
    
        for line in lines:
            if re.match(r"^\s*RELAX\s*=", line):
                relax_found = True
    
                # Extract old values if not provided
                old_relax = re.search(r"RELAX\s*=\s*([0-9]+)", line)
                old_maxit = re.search(r"MAXITER\s*=\s*([0-9]+)", line)
                old_ftol  = re.search(r"RLXFTOL\s*=\s*([-0-9.Ee]+)", line)
                old_type  = re.search(r"RELAXTYPE\s*=\s*(\w+)", line)
    
                if relax_val   is None and old_relax: relax_val   = int(old_relax.group(1))
                if maxiter_val is None and old_maxit: maxiter_val = int(old_maxit.group(1))
                if rlxftol_val is None and old_ftol:  rlxftol_val = float(old_ftol.group(1))
                if old_type: relax_type = old_type.group(1)
    
                # Construct new formatted line
                relax_val   = 1 if relax_val is None else int(relax_val)
                maxiter_val = 300 if maxiter_val is None else int(maxiter_val)
                rlxftol_val = 1e-3 if rlxftol_val is None else float(rlxftol_val)
    
                new_line = f"RELAX= {relax_val:d} RELAXTYPE= {relax_type} MAXITER= {maxiter_val:d} RLXFTOL= {rlxftol_val:.3f}\n"
                new_lines.append(new_line)
            # ---- RESTART block ----
            elif re.match(r"^\s*RESTART\s*=", line):
                restart_found = True
                old_restart = re.search(r"RESTART\s*=\s*([0-9]+)", line)
                if restart_val is None and old_restart:
                    restart_val = int(old_restart.group(1))
                restart_val = 1 if restart_val is None else int(restart_val)
                new_lines.append(f"RESTART= {restart_val:d}\n")
            else:
                new_lines.append(line)
    
        if not relax_found:
            # Append it at the end if missing
            relax_val   = 1 if relax_val is None else int(relax_val)
            maxiter_val = 300 if maxiter_val is None else int(maxiter_val)
            rlxftol_val = 1e-3 if rlxftol_val is None else float(rlxftol_val)
            new_lines.append(f"\nRELAX= {relax_val:d} RELAXTYPE= {relax_type} MAXITER= {maxiter_val:d} RLXFTOL= {rlxftol_val:.3f}\n")

        # Append RESTART if missing
        if not restart_found:
            restart_val = 1 if restart_val is None else int(restart_val)
            new_lines.append(f"RESTART= {restart_val:d}\n")

        with open(control_path, "w") as f:
            f.writelines(new_lines)

    def _read_relax_from_control(self):
        """
        Read RELAX flag from TBparam/control.in and return True/False.
        Default to True if file not found or unreadable.
        """
        control_path = os.path.join(self.workdir, "TBparam", "control.in")
        if not os.path.exists(control_path):
            return True
        with open(control_path, "r") as f:
            for line in f:
                if re.match(r"^\s*RELAX\s*=", line):
                    try:
                        val = int(line.split("=")[1].strip().split()[0])
                        return bool(val)
                    except Exception:
                        return True
        return True


    def _execute(self):
        """
        Run LATTE in self.workdir, logging stdout+stderr to self.log_path.
        """
        if not os.path.exists(self.exe):
            print(f"LATTE executable not found at {self.exe}")
            self.error = True
            return

        try:
            with open(self.log_path, "w") as fout:
                subprocess.run(
                    [self.exe],
                    stdout=fout,
                    stderr=subprocess.STDOUT,
                    timeout=self.timeout,
                    check=True,
                    cwd=self.workdir,
                )
        except subprocess.TimeoutExpired:
            print(f"LATTE timed out after {self.timeout} s")
            self.error = True
        except subprocess.CalledProcessError as e:
            print(f"LATTE failed with return code {e.returncode}")
            self.error = True

    def _read_output(self):
        """
        Parse total energy from the log file.
        """
        if self.error:
            return

        if not os.path.exists(self.log_path):
            print("LATTE log file not found")
            self.error = True
            return

        Tot_energy = re.compile(r"FREE ENERGY\s*=\s*([-0-9.Ee]+)")

        last_e = None
        with open(self.log_path, "r") as f:
            for line in f:
                m = Tot_energy.search(line)
                if m:
                    last_e = float(m.group(1))

        if last_e is None:
            print("Could not find total energy in LATTE log")
            self.error = True
            return

        self.energy = last_e
        self.energy_per_atom = self.energy / len(self.sites)
        self.optimized = self.relax

    # optional stub for later
    def _read_restart_geometry(self):
        """
        Stub: implement later if you want to read relaxed coords from bl/restart.dat.
        """
        restart_path = os.path.join(self.workdir, "bl", "restart.dat")
        if not os.path.exists(restart_path):
            return
        # parse restart.dat here and update self.lattice & self.cart_coords
        # once you know the exact format
        pass


def latte_single_optimize(
    struc,
    exe="./LATTE_DOUBLE",
    workdir=".",
    log_name="latte.log",
    relax=None,                 # if None: auto-read RELAX from control.in
    timeout=3600,
    control_relax=None,         # int 0/1 -> sets RELAX=
    control_maxiter=None,       # int     -> sets MAXITER=
    control_rlxftol=None,       # float   -> sets RLXFTOL=
    control_restart=None,       # int 0/1 -> sets RESTART=
):
    """
    Returns: (pyxtal_struct, energy_per_atom, cputime, error_flag)
    Currently cputime is always 0.0 (LATTE doesn't give us an easy CPU time).
    """
    calc = LATTECalc(
        struc,
        exe=exe,
        workdir=workdir,
        log_name=log_name,
        relax=relax,
        timeout=timeout,
        control_relax=control_relax,
        control_maxiter=control_maxiter,
        control_rlxftol=control_rlxftol,
        control_restart=control_restart,
    )
    calc.run()

    if calc.error:
        print("LATTE error in latte_single_optimize")
        return None, None, 0.0, True
    else:
        final_pyxtal = calc.to_pyxtal() if calc.pyxtal is None else calc.to_pyxtal()
        return final_pyxtal, calc.energy_per_atom, calc.cputime, False

