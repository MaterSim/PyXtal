"""
Common utlities for global search
"""

import os
import warnings
from time import time

import numpy as np
from random import choice
from ase import units

from pyxtal import pyxtal
from pyxtal.interface.ase_opt import ASE_relax
from pyxtal.interface.charmm import CHARMM
from pyxtal.optimize.benchmark import benchmark
from pyxtal.representation import representation
from pyxtal.symmetry import Group, Hall
from pyxtal.util import ase2pymatgen
from pyxtal.XRD import Similarity, pxrd_refine

warnings.filterwarnings("ignore")

def check_stable_structure(xtal, c_info, w_dir, job_tag, skip_ani, optimizer, disps=[0.5, 5.0], random=False, verbose=False):
    """
    Check the stability of input xtal based on lattice mutation
    """
    comp = xtal.get_1D_comp()
    disp_cell, disp_ang = disps[0], disps[1]
    res = optimizer(xtal, c_info, w_dir, job_tag, skip_ani=skip_ani)
    smiles = [m.smile for m in xtal.molecules]#; print(smiles)
    if res is not None:
        xtal0, eng0 = res["xtal"], res["energy"] #/sum(xtal.numMols)
        rep0 = xtal0.get_1D_representation()
        if verbose: print("Optim ", rep0.to_string(eng0/sum(xtal.numMols)))
        cell0 = np.array(xtal.lattice.encode())

        update = False
        for i, c in enumerate(cell0):
            if i <= 2:
                disps = [-disp_cell, disp_cell]
            else:
                disps = [-disp_ang, disp_ang]
            for disp in disps:
                cell = cell0.copy()
                cell[i] += disp
                x = [[xtal0.group.hall_number] + cell.tolist()]
                wps = [site.encode() for site in xtal0.mol_sites]
                # Add random perturbation
                if random:
                    for wp in wps:
                        for id in range(4, len(wp)-1):
                            wp[id] += 10.0*(np.random.random()-0.5)
                x.extend(wps)#; print(x)
                rep1 = representation(x, smiles); print("Init  ", rep1.to_string())
                xtal1 = rep1.to_pyxtal(composition=comp)#; print(xtal1.lattice)
                res = optimizer(xtal1, c_info, w_dir, job_tag, skip_ani=skip_ani)
                if res is not None:
                    xtal2, eng = res["xtal"], res["energy"] #/sum(xtal1.numMols)
                    if eng < eng0 + 1e-4:
                        xtal0, eng0 = xtal2, eng
                        if verbose:
                            rep2 = xtal0.get_1D_representation()
                            print("Update", rep2.to_string(eng0/sum(xtal0.numMols)))
                        update = True
                        #return xtal0, eng0, True
        return xtal0, eng0, update
    else:
        raise RuntimeError("Error in optimization")


def mutator(xtal, smiles, opt_lat, ref_pxrd=None, dr=0.125, random_state=None):
    """A random mutation."""
    rng = np.random.default_rng(random_state)
    # perturb cell
    comp = xtal.get_1D_comp()
    x = xtal.get_1D_representation().x
    if opt_lat:
        if ref_pxrd is None:
            sg = x[0][0]
            disp_cell = rng.uniform(-1.0, 1.0, len(x[0]) - 1)
            if sg <= 2:  # no change in angles to prevent some bad angles
                disp_cell[3:] = 0
            x[0][1:] *= 1 + dr * disp_cell

            # flip the inclination angle
            if 3 <= sg <= 15 and rng.random() > 0.7 and abs(90 - x[0][-1]) < 15:
                x[0][-1] = 180 - x[0][-1]
        else:
            thetas = [ref_pxrd[0][0], min([35, ref_pxrd[0][-1]])]
            _, x[0][1:], _ = pxrd_refine(xtal, ref_pxrd, thetas)

    # perturb molecules
    for i in range(1, len(x)):
        disp_mol = rng.uniform(-1.0, 1.0, len(x[i]) - 1)
        x[i][:-1] *= 1 + dr * disp_mol
        # change the orientation and torsions
        for j in range(3, len(x[i]) - 1):
            rad_num = rng.random()
            if rad_num < 0.25:
                x[i][j] += rng.choice([45.0, 90.0])
            elif rad_num < 0.5:
                x[i][j] *= -1
    try:
        struc = representation(x, smiles).to_pyxtal(composition=comp)
        for i, molecule in enumerate(xtal.molecules):
            if hasattr(molecule, 'active_sites'):
                struc.molecules[i].active_sites = molecule.active_sites
        return struc
    except:
        print(xtal)
        print("x", x)
        print("smiles", smiles)
        print("comp", comp)
        xtal.to_file("bug.cif")
        print("is_valid_matrix\n", xtal.lattice.get_matrix())
        print("cell_para", xtal.lattice.get_para(degree=True))
        print(x[0])
        #raise RuntimeError("Problem occurs in mutation_lattice")
        return None

def randomizer(
    smiles,
    sgs,
    comp,
    lattice=None,
    block=None,
    num_block=None,
    torsions=None,
    molecules=None,
    sites=None,
    use_hall=False,
    factor=1.1,
    random_state=None,
):
    """
    A random structure generation engine

    Args:
        smiles: e.g. `['CCCCC', 'CC']`
        sgs: e.g. `[2, 4, 14]`
        comp: e.g. `[0.5, 0.5]`
        lattice: pyxtal.Lattice object
        block:
        num_block:
        torsions:
        molecules: pre-specified pyxtal_molecule object

    Returns:
        PyXtal object
    """
    rng = np.random.default_rng(random_state)
    mols = [smi + ".smi" for smi in smiles] if molecules is None else [choice(m) for m in molecules]
    sg = rng.choice(sgs)
    wp = Group(sg, use_hall=True)[0] if use_hall else Group(sg)[0]
    mult = len(wp)
    numIons = [int(c * mult) for c in comp]

    # speed up generation for general positions
    if sites is None and comp[0] >= 1:
        letter = wp.letter
        sites = []
        for c in comp:
            sites.append([str(mult) + letter] * int(c))

    while True:
        xtal = pyxtal(molecular=True)
        if use_hall:
            hn = sg
        else:
            perm = sg > 15
            # For specical site or high symmetry, we only do standard_setting
            hn = Hall(sg).hall_default if min(comp) < 1 or sg > 142 else rng.choice(Hall(sg, permutation=perm).hall_numbers)
        xtal.from_random(
            3,
            hn,
            mols,
            numIons,
            factor=factor,
            block=block,
            num_block=num_block,
            lattice=lattice,
            force_pass=True,
            torsions=torsions,
            sites=sites,
            use_hall=True,
            #random_state=random_state,
        )
        if xtal.valid:
            break

    if xtal.has_special_site():
        try:
            xtal = xtal.to_subgroup()
        except:
            print(xtal)
            print(xtal.to_file())
            raise ValueError("Error in making subgroup")
    return xtal


def optimizer(
    struc,
    atom_info,
    workdir,
    tag = "job_0",
    opt_lat = True,
    calculators = None,
    max_time = 180,
    skip_ani = False,
    pre_opt = False,
):
    """
    Structural relaxation for each individual pyxtal structure.

    Args:
        struc: pyxtal
        workdir: working directory
        calculators: e.g., `['CHARMM', 'GULP']`

    Returns:
        a dictionary with xtal, energy and time
    """
    # Perform pre-relaxation
    if pre_opt:
        struc.optimize_lattice_and_rotation()

    if calculators is None:
        calculators = ["CHARMM"]

    cwd = os.getcwd()
    t0 = time()
    os.chdir(workdir)

    results = None
    for i, calculator in enumerate(calculators):
        if calculator == "CHARMM":
            if i == 0:
                calc = CHARMM(struc, tag, steps=[1000], atom_info=atom_info, debug=True)
                calc.run() #clean=False); import sys; sys.exit()
                # print("CCCCCC", calc.optimized); import sys; sys.exit()
                #print(calc.error)
                #print(calc.struc.lattice, calc.error, calc.structure.energy)
                if calc.error:
                    os.chdir(cwd)
                    return None
                else:
                    if not calc.optimized:  #
                        calc = CHARMM(struc, tag, steps=[500], atom_info=atom_info)
                        calc.run()
                        if calc.error:
                            os.chdir(cwd)
                            return None

            # print(calc.structure)
            steps = [1000, 1000] if opt_lat else [2000]
            calc = CHARMM(calc.structure, tag, steps=steps, atom_info=atom_info)
            calc.run()  # print("Debug", calc.optimized); import sys; sys.exit()

            # only count good struc
            if calc.error:
                os.chdir(cwd)
                return None
            else:
                calc = CHARMM(calc.structure, tag, steps=steps, atom_info=atom_info)
                calc.run()#clean=False)

                if calc.error:
                    os.chdir(cwd)
                    return None
                else:
                    if calc.optlat:  # with bad inclination angles
                        calc = CHARMM(calc.structure, tag, steps=steps, atom_info=atom_info)
                        calc.run()#clean=False)
                        if calc.error:
                            os.chdir(cwd)
                            return None

                # Check if there exists a 2nd FF model for better energy ranking
                if os.path.exists("pyxtal1.prm"):
                    calc = CHARMM(
                        calc.structure,
                        tag,
                        prefix="pyxtal1",
                        steps=[2000],
                        atom_info=atom_info,
                    )
                    calc.run()
                    if calc.error:
                        os.chdir(cwd)
                        return None

        if calc.error:
            os.chdir(cwd)
            return None
        else:
            struc = calc.structure
            struc.resort()
    #print(struc.lattice, calc.structure.energy)
    os.chdir(cwd)

    # density should not be too small
    if not skip_ani:
        stress_tol = 10.0 if len(struc.mol_sites[0].molecule.mol) < 10 else 5.0
        if (
            struc.energy < 9999
            and struc.lattice.is_valid_matrix()
            and struc.check_distance()
            and 0.25 < struc.get_density() < 3.0
        ):
            s = struc.to_ase()
            s = ASE_relax(s, 'ANI', step=50, fmax=0.1, logfile="ase.log")
            eng = s.get_potential_energy()
            stress = max(abs(s.get_stress())) / units.GPa  # print(id, eng, stress)

            t = time() - t0
            if t > max_time:  # struc.to_pymatgen().density < 0.9:
                try:
                    print("!!!!! Long time in ani calculation", t)
                    print(struc.get_1D_representation().to_string())
                    struc.optimize_lattice()
                except:
                    print("Trouble in optLat")
                    return None
            elif stress < stress_tol:
                results = {}
                results["xtal"] = struc
                if eng is not None:
                    results["energy"] = eng  # /sum(struc.numMols)
                else:
                    results["energy"] = 10000
                    return None
                results["time"] = time() - t0
            else:
                print(f"stress is wrong {stress:6.2f}")
                return None
                # print(struc)
                # print(struc.to_file())
                # import sys; sys.exit()
    else:
        results = {}
        results["xtal"] = struc
        results["energy"] = calc.structure.energy
        results["time"] = time() - t0

    return results


def optimizer_par(
    xtals,
    ids,
    mutates,
    job_tags,
    randomizer,
    optimizer,
    smiles,
    block,
    num_block,
    atom_info,
    workdir,
    sgs,
    comp,
    lattice,
    torsions,
    molecules,
    sites,
    ref_pmg,
    matcher,
    ref_pxrd,
    use_hall,
    skip_ani,
    check_stable,
    pre_opt,
):
    """
    A routine used for parallel structure optimization

    Args:
        xtals: list of xtals
        ids: list of structure ids
    """
    results = []
    for i in range(len(ids)):
        xtal = xtals[i]
        id = ids[i]
        mutate = mutates[i]
        job_tag = job_tags[i]
        xtal, match = optimizer_single(
            xtal,
            id,
            mutate,
            job_tag,
            randomizer,
            optimizer,
            smiles,
            block,
            num_block,
            atom_info,
            workdir,
            sgs,
            comp,
            lattice,
            torsions,
            molecules,
            sites,
            ref_pmg,
            matcher,
            ref_pxrd,
            use_hall,
            skip_ani,
            check_stable,
            pre_opt,
        )
        results.append((id, xtal, match))
    return results


def optimizer_single(
    xtal,
    id,
    mutate,
    job_tag,
    randomizer,
    optimizer,
    smiles,
    block,
    num_block,
    atom_info,
    workdir,
    sgs,
    comp,
    lattice,
    torsions,
    molecules,
    sites,
    ref_pmg,
    matcher,
    ref_pxrd,
    use_hall,
    skip_ani,
    check_stable,
    pre_opt,
):
    """
    A routine used for individual structure optimization

    Args:
        xtal:
        id: structure id
        randomizer:
        optimizer:
    """

    # 1. Obtain the structure model
    opt_lat = lattice is None
    if xtal is None:
        xtal = randomizer(
            smiles,
            sgs,
            comp,
            lattice,
            block,
            num_block,
            torsions,
            molecules,
            sites,
            use_hall,
        )
        tag = "Random  "
    else:
        if mutate:
            tag = "Mutation"
            xtal = mutator(xtal, smiles, opt_lat, None)#ref_pxrd)
        else:
            tag = "QRandom "

    # 2. Optimization
    if xtal is None:
        res = None
    else:
        res = optimizer(xtal, atom_info, workdir, job_tag, opt_lat, skip_ani=skip_ani, pre_opt=pre_opt)

    # 3. Check match w.r.t the reference
    match = False
    if res is not None:
        xtal, eng = res["xtal"], res["energy"]
        if check_stable and eng < 9999.:
            xtal, eng, status = check_stable_structure(xtal, atom_info, workdir, job_tag, skip_ani, optimizer)
            if status:
                xtal, eng, status = check_stable_structure(xtal, atom_info, workdir, job_tag, skip_ani, optimizer)

        rep = xtal.get_1D_representation()
        N = sum(xtal.numMols)
        strs = rep.to_string(None, eng / N, tag)  # print(strs)

        if ref_pmg is not None:
            pmg_s1 = xtal.to_pymatgen()
            pmg_s1.remove_species("H")
            # Prevent the false call of matcher call
            try:
                rmsd = matcher.get_rms_dist(ref_pmg, pmg_s1)
            except:
                rmsd = None
            if rmsd is not None:
                # Further refine the structure
                match = True
                str1 = f"Match {rmsd[0]:6.2f} {rmsd[1]:6.2f} {eng / N:12.3f} "
                if not skip_ani:
                    xtal, eng1 = refine_struc(xtal, smiles, ASE_relax)
                    str1 += f"Full Relax -> {eng1 / N:12.3f}"
                    eng = eng1
                print(str1)

        elif ref_pxrd is not None:
            thetas = [ref_pxrd[0][0], ref_pxrd[0][-1]]
            xrd = xtal.get_XRD(thetas=thetas)
            p1 = xrd.get_profile(res=0.15, user_kwargs={"FWHM": 0.25})
            match = Similarity(p1, ref_pxrd, x_range=thetas).value
            strs += f" {match:.3f}"

        xtal.energy = eng
        print(f"{id:3d} " + strs)#; print(xtal, eng); import sys; sys.exit()
        return xtal, match
    else:
        return None, match


def refine_struc(xtal, smiles, calculator):
    """
    refine the structure with the ML calculator

    Args:
        - xtal: pyxtal structure
        - calculator: ANI_relax or MACE_relax
    """
    s = xtal.to_ase()
    s = calculator(s, 'ANI', step=50, fmax=0.1, logfile="ase.log")
    s = calculator(s, 'ANI', step=250, opt_cell=True, logfile="ase.log")
    s = calculator(s, 'ANI', step=50, fmax=0.1, logfile="ase.log")
    eng1 = s.get_potential_energy()  # /sum(xtal.numMols)

    xtal = pyxtal(molecular=True)
    pmg = ase2pymatgen(s)
    mols = [smi + ".smi" for smi in smiles]
    xtal.from_seed(pmg, mols)
    return xtal, eng1


def compute_par(row, pmg, work_dir, skf_dir, queue, compute):
    """
    Args:
        xtal:
        queue
    """
    data = compute(row, pmg, work_dir, skf_dir)
    queue.put((row.id, row.csd_code, data))


def compute(row, pmg, work_dir, skf_dir, info=None):
    """
    perform the benchmark for a ase db row object

    Arg:
        row: ase db
    """
    data = {
        "charmm_info": {},
        "gulp_info": {},
        "ani_info": {},
        "dftb_D3_info": {},
        "dftb_TS_info": {},
    }

    if info is None:
        c_info = row.data["charmm_info"]
        g_info = row.data["gulp_info"]
    else:
        c_info = info["c_info"]
        g_info = info["g_info"]
    data["charmm_info"] = c_info
    data["gulp_info"] = g_info

    w_dir = work_dir + "/" + row.csd_code
    smiles = row.mol_smi.split(".")

    # Invoke benchmark
    ben = benchmark(
        pmg,
        smiles,
        work_dir=w_dir,
        skf_dir=skf_dir,
        charmm_info=c_info,
        gulp_info=g_info,
    )
    file1, file2 = "/pyxtal.prm", "/pyxtal.rtf"
    with open(w_dir + file1, "w") as prm:
        prm.write(c_info["prm"])
    with open(w_dir + file2, "w") as rtf:
        rtf.write(c_info["rtf"])

    # reference
    print("\n", row.csd_code, row.mol_smi, pmg.volume)
    rep = representation(ben.rep["reference"], smiles)
    print(rep.to_string() + " reference")

    for calc in ["charmm", "gulp", "ani", "dftb_D3", "dftb_TS"]:
        try:
            ben.calc(calc, show=True)
            key = calc + "_info"
            data[key]["rep"] = ben.rep[calc]
            data[key]["time"] = ben.time[calc]
            data[key]["diff"] = ben.diff[calc]
            data[key]["cif"] = ben.xtal[calc].to_file()
            data[key]["energy"] = ben.energy[calc] / sum(ben.xtal[calc].numMols)
        except:  # noqa: PERF203
            print("=====Calculation is wrong ", calc, row.csd_code)
            data[key]["energy"] = 100000

    return data

def load_reference_from_db(db_name, code=None):
    """
    Load the reference data from the db file

    Args:
        db_name (str): database path
        code (str): code name

    Returns:
        a list of args
    """
    from pyxtal.db import database
    db = database(db_name)
    if code is None:
        codes = db.codes
    else:
        codes = [code]

    args = []
    for code in codes:
        wdir = code
        row = db.get_row(code)
        xtal = db.get_pyxtal(code)
        smile, wt, spg = row.mol_smi, row.mol_weight, row.space_group.replace(" ", "")

        chm_info = None
        if 'charmm_info' in row.data.keys():
            print("prepare charmm input", wdir)
            os.makedirs(wdir, exist_ok=True)
            os.makedirs(wdir+'/calc', exist_ok=True)

            chm_info = row.data['charmm_info']
            prm = open(wdir+'/calc/pyxtal.prm', 'w'); prm.write(chm_info['prm']); prm.close()
            rtf = open(wdir+'/calc/pyxtal.rtf', 'w'); rtf.write(chm_info['rtf']); rtf.close()

        if xtal.has_special_site(): xtal = xtal.to_subgroup()

        pmg0 = xtal.to_pymatgen()
        sg = xtal.group.number
        comp = xtal.get_zprime(integer=True)
        N_torsion = xtal.get_num_torsions()
        lat = xtal.lattice
        tag = code.lower()
        args.append((smile, wdir, sg, tag, chm_info, comp, lat, pmg0, wt, spg, N_torsion))
    return args

if __name__ == "__main__":

    import pymatgen.analysis.structure_matcher as sm
    from pyxtal.db import database

    w_dir = "tmp"
    if not os.path.exists(w_dir):
        os.makedirs(w_dir)

    db = database("pyxtal/database/test.db")
    row = db.get_row("ACSALA")
    xtal1 = db.get_pyxtal("ACSALA")
    smile = row.mol_smi

    # Prepare prm
    c_info = row.data["charmm_info"]
    with open(w_dir + "/pyxtal.prm", "w") as prm:
        prm.write(c_info["prm"])
    with open(w_dir + "/pyxtal.rtf", "w") as rtf:
        rtf.write(c_info["rtf"])

    # Relax expt. xtal
    res = optimizer(xtal1, c_info, w_dir)
    xtal1, eng = res["xtal"], res["energy"]
    rep0 = xtal1.get_1D_representation()
    print(rep0.to_string(eng))

    # Redo it from the 1D. Rep.
    rep2 = representation.from_string(rep0.to_string(), [smile])
    xtal2 = rep2.to_pyxtal()
    res = optimizer(xtal2, c_info, w_dir)
    xtal2, eng = res["xtal"], res["energy"]
    rep0 = xtal2.get_1D_representation()
    print(rep0.to_string(eng))

    pmg1 = xtal1.to_pymatgen()
    pmg2 = xtal2.to_pymatgen()
    pmg1.remove_species("H")
    pmg2.remove_species("H")
    print(sm.StructureMatcher().fit(pmg1, pmg2))

    reps = [
        "81  9.71  6.19 14.25  84.7 1 0 0.21 0.44 0.12  169.6  -16.9  176.2   77.6    9.6   24.9 0",
        "81  8.38 10.06 11.10 107.8 1 0 0.26 0.42 0.31  118.8  -22.6 -111.9 -117.3    0.4   11.9 0",
        "82  9.37  7.92 12.13 111.0 1 0 0.29 0.34 0.10  155.5  -27.6 -161.1   74.6   10.7 -149.8 0",
    ]
    for rep in reps:
        rep = representation.from_string(rep, [smile])
        xtal1 = rep.to_pyxtal()
        check_stable_structure(xtal1, c_info, w_dir, skip_ani=True, optimizer=optimizer)
"""
 81 11.38  6.48 11.24  96.9 1 0 0.23 0.43 0.03  -44.6   25.0   34.4  -76.6   -5.2  171.5 0 -70594.48
 81 11.38  6.48 11.24  96.9 1 0 0.23 0.43 0.03  -44.6   25.0   34.4  -76.6   -5.2  171.5 0 -70594.48
True
"""
