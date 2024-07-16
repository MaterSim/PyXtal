import os
from time import time

import pymatgen.analysis.structure_matcher as sm

from pyxtal import pyxtal
from pyxtal.interface.charmm import CHARMM
from pyxtal.interface.dftb import DFTB, DFTB_relax
from pyxtal.interface.gulp import GULP_OC as GULP_relax
from pyxtal.msg import ReadSeedError
from pyxtal.representation import representation
from pyxtal.util import ase2pymatgen


class benchmark:
    """
    Automated benchmark for experimental crystals with
        - VASP
        - DFTB
        - ANI
        - GULP
        - CHARMM
    """

    def __init__(self, struc, smiles, clean=True, work_dir="tmp", **kwargs):
        """
        Args:
            struc: pymatgen structure or filename
            smiles: list of smile codes
            skf_dir (optional): DFTB parameter directory
            charmm_info (optional): list of atom labels
            charmm_prm (optional): charmm prm file
            charmm_rtf (optional): charmm rtf file
            gulp_info (optional): gulp charge/label info
            clean: whether or not clean up tmp files
        """
        self.clean = clean
        self.valid = True
        self.work_dir = work_dir
        if not os.path.exists(work_dir):
            os.makedirs(work_dir)

        for i in range(len(smiles)):
            if not smiles[i].endswith(".smi"):
                smiles[i] = smiles[i] + ".smi"
        self.smiles = smiles

        xtal = pyxtal(molecular=True)
        try:
            try:
                xtal.from_seed(struc, molecules=smiles)
            except ReadSeedError:
                xtal.from_seed(struc, molecules=smiles, add_H=True)

            # Check if it needs to transform to subgroup representation
            if xtal.has_special_site():
                xtal = xtal.to_subgroup()
            # initialize dicts
            self.xtal = {}
            self.rep = {}
            self.time = {}
            self.energy = {}
            self.pmg = {}
            self.ase = {}
            self.diff = {}
            self.Z = sum(xtal.numMols)
            self.xtal["reference"] = xtal
            self.rep["reference"] = representation.from_pyxtal(xtal).x
            self.ase["reference"] = xtal.to_ase(resort=True)
            self.pmg["reference"] = xtal.to_pymatgen()
            self.pmg["reference"].remove_species("H")  # remove hydrogen
        except ReadSeedError:
            print("Fail to read crystal")
            self.valid = False

        # set up the optional parameters
        if "skf_dir" in kwargs:
            self.skf_dir = kwargs.pop("skf_dir")

        if "charmm_info" in kwargs:
            self.charmm_info = kwargs.pop("charmm_info")

        if "charmm_prm" in kwargs:
            self.charmm_prm = kwargs.pop("charmm_prm")

        if "charmm_rtf" in kwargs:
            self.charmm_rtf = kwargs.pop("charmm_rtf")

        if "gulp_info" in kwargs:
            self.gulp_info = kwargs.pop("gulp_info")

    def calc(self, calculator, show=True):
        """
        execute calculation
        """

        cwd = os.getcwd()
        os.chdir(self.work_dir)
        if calculator.find("dftb") > -1:
            tmp = calculator.split("_")
            disp = tmp[-1] if len(tmp) > 1 else "D3"
            self.dftb(disp, show)
        elif calculator == "vasp":
            self.vasp(show)
        elif calculator == "ani":
            self.ani(show)
        elif calculator == "charmm":
            self.charmm(show)
        elif calculator == "gulp":
            self.gulp(show)
        else:
            raise KeyError("unknow calculator", calculator)

        os.chdir(cwd)

    def dftb(self, disp, show=True, logfile="ase.log"):
        """
        DFTB calculation
        """
        if not hasattr(self, "skf_dir"):
            raise KeyError("skf_dir is not defined for DFTB calculation")
        t0 = time()
        ase = self.ase["reference"].copy()
        # ase = dftb_relax(ase, self.skf_dir, kresol=0.08, logfile=logfile)
        ase, _ = DFTB(ase, self.skf_dir, mode="relax", kresol=0.08, disp=disp)
        ase, _ = DFTB(ase, self.skf_dir, mode="vc-relax", step=300, kresol=0.06, disp=disp)
        ase = DFTB_relax(ase, self.skf_dir, opt_cell=True, kresol=0.06, disp=disp, logfile=logfile)
        xtal = pyxtal(molecular=True)
        pmg = ase2pymatgen(ase)
        xtal.from_seed(pmg, molecules=self.smiles)

        pmg.remove_species("H")
        calc = "dftb_" + disp

        self.xtal[calc] = xtal
        self.pmg[calc] = pmg
        self.energy[calc] = ase.get_potential_energy()
        self.rep[calc] = representation.from_pyxtal(xtal).x
        self.time[calc] = time() - t0
        if show:
            self.summary(calc)

    def vasp(self, show=True):
        """
        VASP calculation
        """

        t0 = time()
        ase, energy = vasp_relax(self.ase["reference"].copy())
        ase, energy = vasp_relax(ase, opt_cell=True)
        xtal = pyxtal(molecular=True)
        pmg = ase2pymatgen(ase)
        xtal.from_seed(pmg, molecules=self.smiles)

        pmg.remove_species("H")
        self.xtal["vasp"] = xtal
        self.pmg["vasp"] = pmg
        self.rep["vasp"] = representation.from_pyxtal(xtal).x
        self.energy["vasp"] = energy
        self.time["vasp"] = time() - t0
        if show:
            self.summary("vasp")

    def ani(self, show=True, logfile="ase-ani.log"):
        """
        Torch ANI calculation
        """

        t0 = time()
        # self.ase.write('ani.cif', format='cif')
        ase = self.ase["reference"].copy()
        ase = ani_relax(ase, logfile=logfile)
        ase = ani_relax(ase, opt_cell=True, logfile=logfile, max_time=10.0)
        ase = ani_relax(ase, opt_cell=True, logfile=logfile, max_time=10.0)

        # ase.write('ani_final.cif', format='cif'); import sys; sys.exit()
        xtal = pyxtal(molecular=True)
        self.ase["ani"] = ase
        pmg = ase2pymatgen(ase)
        try:
            xtal.from_seed(pmg, molecules=self.smiles)
            pmg.remove_species("H")
            self.xtal["ani"] = xtal
            self.pmg["ani"] = pmg
            self.rep["ani"] = representation.from_pyxtal(xtal).x
            self.energy["ani"] = ase.get_potential_energy()
            self.time["ani"] = time() - t0
            if show:
                self.summary("ani")
            # print(xtal); xtal.to_file("test.cif"); import sys; sys.exit()
        except ReadSeedError:
            print("Molecular form is broken after relaxation")

    def charmm(self, show=True, steps=None):
        """
        CHARMM-GAFF
        """
        if steps is None:
            steps = [2000, 3000]
        struc = self.xtal["reference"].copy()

        t0 = time()
        calc = CHARMM(struc, "ben", steps=steps, atom_info=self.charmm_info, debug=True)
        calc.run(clean=self.clean)  # clean=False); import sys; sys.exit()
        struc = calc.structure

        pmg = struc.to_pymatgen()
        pmg.remove_species("H")
        self.xtal["charmm"] = struc
        self.pmg["charmm"] = pmg
        self.rep["charmm"] = representation.from_pyxtal(struc).x
        self.energy["charmm"] = calc.structure.energy
        self.time["charmm"] = time() - t0
        if show:
            self.summary("charmm")

    def gulp(self, show=True, step=None, stepmx=None):
        """
        GULP-GAFF
        """
        if stepmx is None:
            stepmx = [0.001, 0.005, 0.02]
        if step is None:
            step = [400, 400, 1000]
        struc = self.xtal["reference"].copy()
        g_info = self.gulp_info
        t0 = time()
        calc = GULP_relax(struc, "ben", opt="conv", steps=step[0], stepmx=stepmx[0], atom_info=g_info)
        calc.run(clean=self.clean)  # ; print(os.getcwd()); import sys; sys.exit()
        if not calc.optimized:
            raise RuntimeError("GULP calculation is wrong")
        struc = calc.structure
        calc = GULP_relax(
            struc,
            "ben",
            opt="conp",
            steps=step[1],
            stepmx=stepmx[1],
            atom_info=g_info,
            dump="1.cif",
        )
        calc.run(clean=self.clean)
        if not calc.optimized:
            raise RuntimeError("GULP calculation is wrong")
        struc = calc.structure
        calc = GULP_relax(struc, "ben", opt="conp", steps=step[2], stepmx=stepmx[2], atom_info=g_info)
        # print(struc)
        calc.run(clean=self.clean)  # , pause=True); import sys; sys.exit()
        struc = calc.structure
        if not calc.optimized:
            raise RuntimeError("GULP calculation is wrong")
        pmg = struc.to_pymatgen()
        pmg.remove_species("H")
        self.xtal["gulp"] = struc
        self.pmg["gulp"] = pmg
        self.rep["gulp"] = representation.from_pyxtal(struc).x
        self.energy["gulp"] = calc.structure.energy
        self.time["gulp"] = time() - t0
        if show:
            self.summary("gulp")

    def charmm_md(self):
        raise NotImplementedError

    def summary(self, calc=None):
        calcs = self.time.keys() if calc is None else [calc]

        pmg0 = self.pmg["reference"]
        v0 = pmg0.volume
        matcher = sm.StructureMatcher(ltol=0.3, stol=0.3, angle_tol=10)
        for calc in calcs:
            time = self.time[calc] / 60
            rep = representation(self.rep[calc], self.smiles)
            dv = (self.pmg[calc].volume - v0) / v0
            strs = f"{rep.to_string(eng=self.energy[calc] / self.Z):s} "
            strs += f"{calc:8s} {time:6.2f} {dv:6.3f}"
            rmsd = matcher.get_rms_dist(pmg0, self.pmg[calc])
            if rmsd is not None:
                strs += f"{rmsd[0]:6.3f}{rmsd[1]:6.3f}"
                self.diff[calc] = [dv, rmsd[0], rmsd[1]]
            else:
                self.diff[calc] = [dv, None, None]
            print(strs)


if __name__ == "__main__":
    import warnings

    from pyxtal.db import database

    warnings.filterwarnings("ignore")

    work_dir = "tmp"
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)

    # if 'DFTB_PREFIX' in os.environ.keys():
    #    skf_dir = os.environ['DFTB_PREFIX'] + '3ob-3-1/'
    # else:
    #    raise RuntimeError("Cannot find DFTB_PREFIX in the environment")
    skf_dir = None

    db = database("pyxtal/database/test.db")
    code = "ACSALA"  #'BENZEN'
    row = db.get_row(code)
    xtal = db.get_pyxtal(code)

    c_info = row.data["charmm_info"]
    with open(work_dir + "/pyxtal.prm", "w") as prm:
        prm.write(c_info["prm"])
    with open(work_dir + "/pyxtal.rtf", "w") as rtf:
        rtf.write(c_info["rtf"])
    g_info = row.data["gulp_info"]

    pmg = xtal.to_pymatgen()
    smi = row.mol_smi.split(".")
    ben = benchmark(
        pmg,
        smi,
        charmm_info=c_info,
        gulp_info=g_info,
        work_dir=work_dir,
        skf_dir=skf_dir,
        clean=False,
    )

    rep = representation(ben.rep["reference"], smi)
    print(rep.to_string() + " reference")

    # for calc in ['charmm', 'ani', 'gulp']: #, 'dftb_TS']:
    # for calc in ['dftb_TS']:
    # for calc in ['ani', 'ani', 'ani', 'ani']:
    for calc in ["charmm", "gulp"]:  # , 'dftb_TS']:
        ben.calc(calc, show=True)
    print("=========================================")
"""
qzhu@cms-2:~/GitHub/HT-OCSP$ python htocsp/benchmark.py
 29 0  9.20  7.29  6.69 1 0.00 0.75 0.00  111.8  -43.7   19.4 1  reference
 29 0  9.20  7.44  6.68 1 0.00 0.75 1.00  -72.0  -42.7   20.6 0      -1.444 charmm     0.15  0.019 0.030 0.032
 61 0  6.63  7.42  9.36 1 1.00 0.00 1.00   -7.5   42.5 -108.3 0   -6318.357 ani        0.24  0.027 0.008 0.010
 29 0  9.24  7.35  6.64 1 0.00 0.75 1.00  -70.6  -43.2   19.4 0      -0.305 gulp       0.60  0.006 0.024 0.028
 61 0  6.56  7.12  9.11 1 0.00 1.00 0.00   -7.8   43.2 -109.0 0    -341.479 dftb_D3    2.31 -0.051 0.017 0.022
 61 0  6.70  7.24  9.20 1 0.00 0.00 0.00 -142.9   43.3  -70.8 0    -341.147 dftb_TS    0.82 -0.005 0.011 0.013
"""
