"""
Global Optimizer
"""

from __future__ import annotations

from time import time
from typing import TYPE_CHECKING

import numpy as np
from scipy.stats import qmc

# import threading
import psutil
from numpy.random import Generator
from pymatgen.analysis.structure_matcher import StructureMatcher

from pyxtal.optimize.base import GlobalOptimize
from pyxtal.representation import representation
from pyxtal.lattice import Lattice

if TYPE_CHECKING:
    from pyxtal.molecule import pyxtal_molecule

def generate_qrs_cell(sampler, cell_bounds, ref_volume, ltype):
    """
    A routine to generate quasi random samples for lattice and wp
    """
    # Sample cell parameters
    min_vol, max_vol = 0.75*ref_volume, 2.5*ref_volume
    lb = [b[0] for b in cell_bounds]
    ub = [b[1] for b in cell_bounds]
    count = 0
    while True:
        count += 1
        sample = qmc.scale(sampler.random(), lb, ub)[0].tolist()
        lat = Lattice.from_1d_representation(sample, ltype)
        if min_vol < lat.volume < max_vol:
            #print(sample, ltype)
            return sample
        if count == 1000:
            raise ValueError("Cannot generate valid cell with 1000 attempts")

def generate_qrs_xtals(cell, wp_bounds, N_pop, smiles):
    """
    Get the qrs xtal samples
    """
    #cell = [81, 11.38,  6.48, 11.24,  96.9]
    m = max([int(np.log2(N_pop))+3, 9])
    xtals = []
    lb = [b[0] for b in wp_bounds]
    ub = [b[1] for b in wp_bounds]
    sampler_wp = qmc.Sobol(d=len(wp_bounds), scramble=False)
    sample_wp = sampler_wp.random_base2(m=m)
    sample_wp = qmc.scale(sample_wp, lb, ub)
    for wp0 in sample_wp:
        x = [cell, [0] + wp0.tolist() + [False]] #print(x)
        rep = representation(x, smiles)
        xtal = rep.to_pyxtal()
        if not xtal.has_special_site() and len(xtal.check_short_distances(r=0.8)) == 0:
            #print(rep, len(xtal.check_short_distances(r=0.8)))
            xtals.append(xtal)
            if len(xtals) == N_pop:
                return xtals
        #else:
        #    print(rep, len(xtal.check_short_distances(r=0.6)))
    print("DDDDDD", len(xtals), m)
    return xtals


class QRS(GlobalOptimize):
    """
    Quasi Monte Carlo

    Args:
        smiles (str): smiles string
        workdir (str): path of working directory
        sg (int or list): space group number or list of spg numbers
        tag (string): job prefix
        ff_opt (bool): activate on the fly FF mode
        ff_style (str): automated force style (`gaff` or `openff`)
        ff_parameters (str or list): ff parameter xml file or list
        reference_file (str): path of reference xml data for FF training
        N_gen (int): number of generation (default: `10`)
        N_pop (int): number of populations (default: `10`)
        fracs (list): fractions for each variation (default: `[0.5, 0.5, 0.0]`)
        N_cpu (int): number of cpus for parallel calculation (default: `1`)
        cif (str): cif file name to store all structure information
        block: block mode
        num_block: list of blocks
        compositions: list of composition, (default is [1]*Num_mol)
        lattice (bool): whether or not supply the lattice
        torsions: list of torsion angle
        molecules (list): list of pyxtal_molecule objects
        sites (list): list of wp sites, e.g., [['4a']]
        use_hall (bool): whether or not use hall number (default: False)
        skip_ani (bool): whether or not use ani or not (default: True)
        eng_cutoff (float): the cutoff energy for FF training
        E_max (float): maximum energy defined as an invalid structure
        verbose (bool): show more details
    """

    def __init__(
        self,
        smiles: str,
        workdir: str,
        sg: int | list,
        tag: str,
        info: dict[any, any] | None = None,
        ff_opt: bool = False,
        ff_style: str = "openff",
        ff_parameters: str = "parameters.xml",
        reference_file: str = "references.xml",
        ref_criteria: dict[any, any] | None = None,
        N_gen: int = 10,
        N_pop: int = 10,
        N_cpu: int = 1,
        cif: str | None = None,
        block: list[any] | None = None,
        num_block: list[any] | None = None,
        composition: list[any] | None = None,
        lattice: Lattice | None = None,
        torsions: list[any] | None = None,
        molecules: list[pyxtal_molecule] | None = None,
        sites: list[any] | None = None,
        use_hall: bool = False,
        skip_ani: bool = True,
        factor: float = 1.1,
        eng_cutoff: float = 5.0,
        E_max: float = 1e10,
        N_survival: int = 20,
        verbose: bool = False,
        random_state: int | None = None,
        max_time: float | None = None,
        matcher: StructureMatcher | None = None,
        early_quit: bool = True,
    ):

        self.N_gen = N_gen # Number of lattice points
        self.N_pop = N_pop # Number of wp varieties
        self.verbose = verbose
        # initialize other base parameters
        GlobalOptimize.__init__(
            self,
            smiles,
            workdir,
            sg,
            tag,
            info,
            ff_opt,
            ff_style,
            ff_parameters,
            reference_file,
            ref_criteria,
            N_cpu,
            cif,
            block,
            num_block,
            composition,
            lattice,
            torsions,
            molecules,
            sites,
            use_hall,
            skip_ani,
            factor,
            eng_cutoff,
            E_max,
            None, #random_state,
            max_time,
            matcher,
            early_quit,
        )

        strs = self.full_str()
        self.logging.info(strs)
        print(strs)

    def full_str(self):
        s = str(self)
        s += "\nMethod    : Deterministic Quasi-Random Sampling"
        s += f"\nGeneration: {self.N_gen:4d}"
        s += f"\nPopulation: {self.N_pop:4d}"
        # The rest base information from now on
        return s

    def run(self, ref_pmg=None, ref_eng=None, ref_pxrd=None):
        """
        The main code to run GA prediction

        Args:
            ref_pmg: reference pmg structure
            ref_eng: reference energy
            ref_pxrd: reference pxrd profile in 2D array

        Returns:
            (generation, np.min(engs), None, None, None, 0, len(engs))
        """
        self.ref_volumes = []
        if ref_pmg is not None:
            ref_pmg.remove_species("H")

        # Related to the FF optimization
        N_added = 0

        # To save for comparison
        current_survivals = [0] * self.N_pop # track the survivals
        hist_best_xtals = [None] * self.N_pop
        hist_best_engs = [self.E_max] * self.N_pop

        # lists for structure information
        current_reps = [None] * self.N_pop
        current_matches = [False] * self.N_pop if ref_pxrd is None else [0.0] * self.N_pop
        current_engs = [self.E_max] * self.N_pop
        current_tags = ["Random"] * self.N_pop

        for gen in range(self.N_gen):
            print(f"\nGeneration {gen:d} starts")
            self.logging.info(f"Generation {gen:d} starts")
            self.generation = gen + 1
            t0 = time()

            if gen > 0:
                cell = generate_qrs_cell(self.sampler,
                                         self.cell_bounds,
                                         self.ref_volumes[-1],
                                         self.ltype)
                cell = [self.hall_number] + cell
                current_xtals = generate_qrs_xtals(cell, self.wp_bounds, self.N_pop, self.smiles)
                strs = f"Cell parameters in Gen{gen:d}: "
                print(strs, cell, len(current_xtals))
            else:
                # 1st generation from random
                current_xtals = [None] * self.N_pop

            # Local optimization
            gen_results = self.local_optimization(gen, current_xtals, ref_pmg, ref_pxrd)

            # Summary and Ranking
            for id, res in enumerate(gen_results):
                (xtal, match) = res

                if xtal is not None:
                    current_xtals[id] = xtal
                    rep = xtal.get_1D_representation()
                    current_reps[id] = rep.x
                    current_engs[id] = xtal.energy / sum(xtal.numMols)
                    current_matches[id] = match
                    if current_engs[id] < hist_best_engs[id]:
                        hist_best_engs[id] = current_engs[id]
                        hist_best_xtals[id] = xtal

                    # Don't write bad structure
                    if self.cif is not None and xtal.energy < self.E_max:
                        if self.verbose:
                            print("Add qualified structure", id, xtal.energy)
                        with open(self.workdir + "/" + self.cif, "a+") as f:
                            label = self.tag + "-g" + str(gen) + "-p" + str(id)
                            f.writelines(xtal.to_file(header=label))
                        # else:
                        #    print("Neglect bad structure", id, xtal.energy)
                    self.engs.append(xtal.energy / sum(xtal.numMols))
                    # print(output)

            strs = f"Generation {gen:d} finishes"  # ; import sys; sys.exit()
            print(strs)
            self.logging.info(strs)
            t1 = time()

            # Store the best structures
            engs = current_engs
            count = 0
            xtals = []
            vols = []
            ids = np.argsort(engs)
            for id in ids:
                xtal = current_xtals[id]
                rep = current_reps[id]
                eng = current_engs[id]
                tag = current_tags[id]
                if self.new_struc(xtal, xtals):
                    xtals.append(xtal)
                    self.best_reps.append(rep)
                    d_rep = representation(rep, self.smiles)
                    strs = d_rep.to_string(None, eng, tag)
                    out = f"{gen:3d} {strs:s} Top"
                    if ref_pxrd is not None:
                        out += f" {current_matches[id]:6.3f}"
                    print(out)
                    count += 1
                    vols.append(xtal.lattice.volume)
                if count == 3:
                    break

            # update best volume
            self.ref_volumes.append(np.array(vols).mean())
            if gen == 0:
                best_xtal = xtals[0]
                self.cell_bounds = best_xtal.lattice.get_bounds(2.5, 25)
                self.ltype = best_xtal.lattice.ltype
                self.wp_bounds = best_xtal.mol_sites[0].get_bounds()
                self.hall_number = best_xtal.group.hall_number
                self.sampler = qmc.Sobol(d=len(self.cell_bounds), scramble=False)

            t2 = time()
            gen_out = f"Gen{gen:3d} time usage: "
            gen_out += f"{t1 - t0:5.1f}[Calc] {t2 - t1:5.1f}[Proc]"
            print(gen_out)

            # Save the reps for next move
            prev_xtals = current_xtals  # ; print(self.engs)
            self.min_energy = np.min(np.array(self.engs))
            self.N_struc = len(self.engs)

            # Update the FF parameters if necessary
            if self.ff_opt:
                N_max = min([int(self.N_pop * 0.6), 50])
                ids = np.argsort(engs)
                xtals = self.select_xtals(current_xtals, ids, N_max)
                print("Select Good structures for FF optimization", len(xtals))
                N_added = self.ff_optimization(xtals, N_added)

            else:
                success_rate = self.success_count(gen, current_xtals, current_matches, current_tags, engs, ref_pmg)
                gen_out = f"Success rate at Gen {gen:3d}: "
                gen_out += f"{success_rate:7.4f}%"
                self.logging.info(gen_out)
                print(gen_out)

                if self.early_termination(success_rate):
                    return success_rate

        return None

if __name__ == "__main__":
    import argparse
    import os

    from pyxtal.db import database

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-g",
        "--gen",
        dest="gen",
        type=int,
        default=10,
        help="Number of generation, optional",
    )
    parser.add_argument(
        "-p",
        "--pop",
        dest="pop",
        type=int,
        default=10,
        help="Population size, optional",
    )
    parser.add_argument("-n", "--ncpu", dest="ncpu", type=int, default=1, help="cpu number, optional")
    parser.add_argument("--ffopt", action="store_true", help="enable ff optimization")

    options = parser.parse_args()
    gen = options.gen
    pop = options.pop
    ncpu = options.ncpu
    ffopt = options.ffopt
    db_name, name = "pyxtal/database/test.db", "ACSALA"
    wdir = name
    os.makedirs(wdir, exist_ok=True)
    os.makedirs(wdir + "/calc", exist_ok=True)

    db = database(db_name)
    row = db.get_row(name)
    xtal = db.get_pyxtal(name)
    smile, wt, spg = row.mol_smi, row.mol_weight, row.space_group.replace(" ", "")
    chm_info = None
    if not ffopt:
        if "charmm_info" in row.data:
            # prepare charmm input
            chm_info = row.data["charmm_info"]
            with open(wdir + "/calc/pyxtal.prm", "w") as prm:
                prm.write(chm_info["prm"])
            with open(wdir + "/calc/pyxtal.rtf", "w") as rtf:
                rtf.write(chm_info["rtf"])
        else:
            # Make sure we generate the initial guess from ambertools
            if os.path.exists("parameters.xml"):
                os.remove("parameters.xml")

    # load reference xtal
    pmg0 = xtal.to_pymatgen()
    if xtal.has_special_site():
        xtal = xtal.to_subgroup()

    # GO run
    t0 = time()
    go = QRS(
        smile,
        wdir,
        xtal.group.number,
        name.lower(),
        info=chm_info,
        ff_style="openff",  #'gaff',
        ff_opt=ffopt,
        N_gen=gen,
        N_pop=pop,
        N_cpu=ncpu,
        cif="pyxtal.cif",
    )

    suc_rate = go.run(pmg0)
    print(strs + " in Gen {:d}\n".format(go.generation))

    if len(go.matches) > 0:
        best_rank = go.print_matches()
        mytag = f"True {best_rank:d}/{go.N_struc:d} Succ_rate: {suc_rate:7.4f}%"
    else:
        mytag = f"False 0/{go.N_struc:d}"

    eng = go.min_energy
    t1 = int((time() - t0)/60)
    strs = "Final {:8s} [{:2d}]{:10s} ".format(code, sum(xtal.numMols), spg)
    strs += "{:3d}m {:2d} {:6.1f}".format(t1, N_torsion, wt)
    strs += "{:12.3f} {:20s} {:s}".format(eng, mytag, smile)
    print(strs)
