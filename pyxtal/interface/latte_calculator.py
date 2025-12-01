from ase.io import read, write
from ase.optimize import FIRE
from ase.constraints import UnitCellFilter
from ase.calculators.calculator import Calculator, all_changes

from collections import Counter
from pathlib import Path
import numpy as np
import os
import subprocess
import re

# ===================== 1. LATTE Calculator =====================

class LATTECalculator(Calculator):
    """
    Minimal ASE calculator that:
       1) writes bl/inputblock.dat
       2) runs LATTE
       3) reads energy, stress, forces from
              latte.log
              lastsystem.cfg  (final frame only)

    Requirements in workdir:
       ./LATTE_DOUBLE
       TBparam/control.in
       MDcontroller
       bl/ (auto-created)
    """

    implemented_properties = ["energy", "forces", "stress"]
    
    def __init__(self, exe="./LATTE_DOUBLE", workdir=".", timeout=3600, **kwargs):
        super().__init__(**kwargs)
        self.workdir = os.path.abspath(workdir)
        self.exe = os.path.abspath(os.path.join(self.workdir, exe))
        self.timeout = timeout
        self.logfile = os.path.join(self.workdir, "latte.log")
        self.cfgfile = os.path.join(self.workdir, "lastsystem.cfg")
        
    def calculation_required(self, atoms, properties):
        return True
        
    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        if atoms is None:
            atoms = self.atoms
        else:
            self.atoms = atoms.copy()

        # clear previous results
        self.results = {}

        # write input for *current* geometry
        self._write_inputblock(self.atoms)
        self._set_single_point_mode()   # RELAX=0 in control.in

        # run LATTE and parse outputs
        self._run_latte()

        E = self._read_energy_from_log()
        F = self._read_forces_from_cfg()
        S = self._read_stress_from_log()

        self.results["energy"] = E
        self.results["forces"] = F
        if S is not None:
            self.results["stress"] = S

    # Write inputblock.dat
    def _write_inputblock(self, atoms):
        bl = os.path.join(self.workdir, "bl")
        os.makedirs(bl, exist_ok=True)
        fpath = os.path.join(bl, "inputblock.dat")

        cell = atoms.cell.array
        pos = atoms.get_positions()
        syms = atoms.get_chemical_symbols()

        with open(fpath, "w") as f:
            f.write(f"{len(atoms)}\n")
            for v in cell:
                f.write(f"{v[0]:.10f} {v[1]:.10f} {v[2]:.10f}\n")
            for s, r in zip(syms, pos):
                f.write(f"{s} {r[0]:.10f} {r[1]:.10f} {r[2]:.10f}\n")

    # Force single-point mode in control.in
    def _set_single_point_mode(self):
        cpath = os.path.join(self.workdir, "TBparam", "control.in")
        if not os.path.exists(cpath):
            print("WARNING: No control.in found.")
            return

        lines = open(cpath).read().splitlines()
        new = []
        updated_relax = False

        for line in lines:
            if line.strip().startswith("RELAX"):
                updated_relax = True
                new.append("RELAX= 0 RELAXTYPE= SD MAXITER= 10 RLXFTOL= 0.1\n")
            else:
                new.append(line + "\n")

        if not updated_relax:
            new.append("\nRELAX= 0 RELAXTYPE= SD MAXITER= 10 RLXFTOL= 0.1\n")

        open(cpath, "w").write("".join(new))

    # Run LATTE
    def _run_latte(self):
        with open(self.logfile, "w") as f:
            subprocess.run(
                [self.exe],
                cwd=self.workdir,
                stdout=f, stderr=subprocess.STDOUT,
                timeout=self.timeout, check=True
            )

    # Read total energy from latte.log
    def _read_energy_from_log(self):
        lines = open(self.logfile).read().splitlines()
        E = None
        pattern = re.compile(r"FREE ENERGY\s*=\s*([-0-9.Ee]+)")
        for l in lines:
            m = pattern.search(l)
            if m:
                E = float(m.group(1))
        if E is None:
            raise RuntimeError("Energy not found in log!")
        return E

    # Read stress tensor from latte.log (GPa → eV/A³)
    def _read_stress_from_log(self):
        lines = open(self.logfile).read().splitlines()
        S = None
        for i, l in enumerate(lines):
            if "Stress tensor (GPa)" in l:
                m = []
                for j in range(1, 4):
                    parts = lines[i+j].replace("#","").split()
                    m.append([float(x) for x in parts[:3]])
                S = np.array(m)
                break

        if S is None:
            return None

        # convert GPa → eV/Å³ and flip sign (LATTE → ASE convention)
        conv = 0.00625
        return -(S*conv)

    # Read the LAST FRAME from lastsystem.cfg
    def _read_forces_from_cfg(self):
        if not os.path.exists(self.cfgfile):
            raise RuntimeError("lastsystem.cfg not found")

        lines = open(self.cfgfile).read().splitlines()

        # number of atoms
        natoms = None
        for l in lines:
            if "Number of particles" in l:
                natoms = int(l.split("=")[1])
                break

        if natoms is None:
            raise RuntimeError("Could not find 'Number of particles' in lastsystem.cfg")

        # H0 matrix
        H = np.zeros((3, 3))
        for l in lines:
            if "H0(" in l:
                tmp = l.replace("=", " ").replace("A", "").split()
                ij = tmp[0][3:-1]
                i, j = ij.split(",")
                H[int(i) - 1, int(j) - 1] = float(tmp[1])

        # collect all numeric atom-lines (≥7 floats)
        all_atoms = []
        for l in lines:
            p = l.split()
            if len(p) < 7:
                continue
            try:
                nums = [float(x) for x in p[:7]]
                all_atoms.append(nums)
            except ValueError:
                continue

        if len(all_atoms) < natoms:
            raise RuntimeError("Not enough atom lines found in lastsystem.cfg")

        # last natoms rows = final forces
        block = all_atoms[-natoms:]

        frac = np.array([b[:3] for b in block])
        forces = np.array([b[3:6] for b in block])

        coords_cart = frac @ H  # not used but kept for completeness

        return forces

# ===================== 2. Reorder CIF & fix formula =====================

def reorder_cif_and_fix_formula(infile, outfile, desired_order=("O", "H", "C")):
    """
    Read CIF, reorder atoms by desired_order, write new CIF,
    and fix _chemical_formula_sum line. Returns reordered Atoms.
    """
    atoms = read(infile)

    # Reorder atoms
    symbols = np.array(atoms.get_chemical_symbols())
    indices = []
    for el in desired_order:
        idx = np.where(symbols == el)[0]
        indices.extend(idx.tolist())
    atoms_reordered = atoms[indices]

    # Write reordered CIF
    write(outfile, atoms_reordered, format="cif")

    # Count elements
    counts = Counter(atoms_reordered.get_chemical_symbols())
    parts = []
    for el in desired_order:
        if el in counts:
            n = counts[el]
            parts.append(f"{el}{n}" if n != 1 else el)
    custom_formula = " ".join(parts)

    # Patch _chemical_formula_sum
    text = Path(outfile).read_text()
    lines = text.splitlines()

    new_lines = []
    formula_found = False
    for line in lines:
        if line.lower().startswith("_chemical_formula_sum"):
            new_lines.append(f"_chemical_formula_sum    '{custom_formula}'")
            formula_found = True
        else:
            new_lines.append(line)

    if not formula_found:
        inserted = []
        inserted_flag = False
        for line in new_lines:
            inserted.append(line)
            if (not inserted_flag) and line.lower().startswith("data_"):
                inserted.append(f"_chemical_formula_sum    '{custom_formula}'")
                inserted_flag = True
        new_lines = inserted

    Path(outfile).write_text("\n".join(new_lines))

    print("Wrote reordered CIF to:", outfile)
    print("Formula set to:", custom_formula)
    return atoms_reordered

# ===================== 3. Stress logger =====================

def log_stress(atoms, filename="stress_log.txt"):
    """
    Append latest stress tensor from atoms.calc.results["stress"]
    to filename without forcing a new LATTE evaluation.
    """
    calc = atoms.calc
    S = calc.results.get("stress", None)
    if S is None:
        # fall back to ASE call if needed (this will trigger a calc)
        S = atoms.get_stress(voigt=False)

    S = np.array(S)

    first_time = not os.path.exists(filename)
    with open(filename, "a") as f:
        if first_time:
            f.write("# xx yy zz yz xz xy\n")
        if S.shape == (3, 3):
            xx, yy, zz = S[0, 0], S[1, 1], S[2, 2]
            yz, xz, xy = S[1, 2], S[0, 2], S[0, 1]
            f.write(f"{xx: .6f} {yy: .6f} {zz: .6f} {yz: .6f} {xz: .6f} {xy: .6f}\n")
        else:
            # assume already 6-component
            f.write(" ".join(f"{x: .6f}" for x in S) + "\n")

# ===================== 4. Main workflow =====================

if __name__ == "__main__":
    # --- user settings ---
    infile = "bddc_exp.cif"               # Original CIF
    rearranged_cif = "rearranged_bddc.cif"  # Reordered output CIF
    latte_workdir = "."               # Folder with LATTE_DOUBLE, TBparam, MDcontroller

    # 1) Reorder CIF + fix formula, get reordered Atoms
    desired_order = ("O", "H", "C")   #(you can see the order through running with original one for 1 step and observe from "lastsystem.cfg")
    atoms = reorder_cif_and_fix_formula(infile, rearranged_cif, desired_order)

    # 2) Attach LATTE calculator
    calc = LATTECalculator(exe="./LATTE_DOUBLE", workdir=latte_workdir)
    atoms.set_calculator(calc)

    # 3) Cell + atomic relaxation with FIRE
    ucf = UnitCellFilter(atoms)
    opt = FIRE(ucf, dt=0.2, maxmove=0.2, logfile="latte_opt.log")

    # Attach stress logger (uses closure over 'atoms')
    opt.attach(lambda: log_stress(atoms, filename="stress_log.txt"), interval=10)

    opt.run(fmax=0.1)

    # 4) Save relaxed structure
    write("Relaxed_structure.cif", atoms)
    print("Relaxed structure written to Relaxed_structure.cif")

