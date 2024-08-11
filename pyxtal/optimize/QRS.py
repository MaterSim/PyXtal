"""
Global Optimizer
"""

from __future__ import annotations

from time import time
from typing import TYPE_CHECKING

import numpy as np
from scipy.stats import qmc

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

def generate_qrs_xtals(cell, wp_bounds, N_pop, smiles, comp, sampler_wp=None, d_tol=0.85):
    """
    Get the qrs xtal samples

    Args:
        cell (list): [hall, a, b, c]
        wp_bounds (list): [[wp0], [wp1], ...]
        N_pop (int): number of valid candidates
        smiles (list): []
        comp (list): [1]
        sampler_wp: sampler
        d_tol (float): short distance tolerance value
    """
    #cell = [81, 11.38,  6.48, 11.24,  96.9]
    xtals = []
    lb, ub = [], []
    seqs = []
    for wp_bound in wp_bounds:
        lb += [b[0] for b in wp_bound]
        ub += [b[1] for b in wp_bound]
        seqs.append(len(wp_bound))

    if sampler_wp is None:
        sampler_wp = qmc.Sobol(d=len(lb), scramble=False)

    m = max([int(np.log2(N_pop))+3, 9])

    for i in range(2**m):
        sample_wp = sampler_wp.random()#; print(sample_wp)
        sample_wp = qmc.scale(sample_wp, lb, ub)[0].tolist()
        x = [cell]
        prev = 0
        for seq in seqs:
            wp = [0] + sample_wp[prev:prev+seq] + [0]#; print('DDDD', prev, prev+seq, sample_wp[prev:prev+seq], wp)
            x.append(wp) #, [0] + wp0.tolist() + [False]] #print(x)
            prev = seq
        rep = representation(x, smiles)
        xtal = rep.to_pyxtal(composition=comp)
        if not xtal.has_special_site() and len(xtal.check_short_distances(r=d_tol)) == 0:
            #print("debug", rep)
            xtals.append(xtal)
            if len(xtals) == N_pop:
                return xtals
        #else:
        #    print(rep, len(xtal.check_short_distances(r=0.6)))
    return xtals


class QRS(GlobalOptimize):
    """
    Quasi-Random Sampling

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
        tag: str = 'test',
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
        early_quit: bool = False,
        check_stable: bool = False,
    ):

        # POPULATION parameters:
        self.N_gen = N_gen # Number of lattice points
        self.N_pop = N_pop # Number of wp varieties
        self.verbose = verbose
        self.name = 'QRS'

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
            check_stable,
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
                if self.lattice is not None:
                    cell = [self.hall_number] + self.lattice.encode()
                    current_xtals = generate_qrs_xtals(cell,
                                                       self.wp_bounds,
                                                       self.N_pop,
                                                       self.smiles,
                                                       self.composition,
                                                       self.sampler)
                else:
                    cell = generate_qrs_cell(self.sampler,
                                             self.cell_bounds,
                                             self.ref_volumes[-1],
                                             self.ltype)
                    cell = [self.hall_number] + cell
                    current_xtals = generate_qrs_xtals(cell,
                                                       self.wp_bounds,
                                                       self.N_pop,
                                                       self.smiles,
                                                       self.composition)
                strs = f"Cell parameters in Gen-{gen:d}: "
                print(strs, cell, self.ref_volumes[-1], len(current_xtals))
            else:
                 # 1st generation from random
                 current_xtals = [None] * self.N_pop

            # Local optimization
            gen_results = self.local_optimization(gen, current_xtals, ref_pmg, ref_pxrd, True)

            # Summary and Ranking
            current_xtals, current_matches, engs = self.gen_summary(gen,
                    t0, gen_results, current_xtals, current_tags, ref_pxrd)

            # update hist_best
            vols = []
            for id, xtal in enumerate(current_xtals):
                if xtal is not None:
                    vols.append(xtal.lattice.volume)

            # update best volume
            self.ref_volumes.append(np.array(vols).mean())
            if gen == 0:
                best_xtal = xtals[0]
                self.cell_bounds = best_xtal.lattice.get_bounds(2.5, 25)
                self.ltype = best_xtal.lattice.ltype
                self.wp_bounds = [site.get_bounds() for site in best_xtal.mol_sites]
                self.hall_number = best_xtal.group.hall_number
                if self.lattice is not None:
                    len_reps = sum(len(bound) for bound in self.wp_bounds)
                else:
                    len_reps = len(self.cell_bounds)
                self.sampler = qmc.Sobol(d=len_reps, scramble=False)

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
                if ref_pmg is not None:
                    success_rate = self.success_count(gen,
                                                      current_xtals,
                                                      current_matches,
                                                      current_tags,
                                                      ref_pmg)
                    gen_out = f"Success rate at Gen {gen:3d}: {success_rate:7.4f}%"
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
    if xtal.has_special_site(): xtal = xtal.to_subgroup()
    N_torsion = xtal.get_num_torsions()

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
        check_stable=True,
    )

    suc_rate = go.run(pmg0)
    print(f"CSD {name:s} in Gen {go.generation:d}")

    if len(go.matches) > 0:
        best_rank = go.print_matches()
        mytag = f"True {best_rank:d}/{go.N_struc:d} Succ_rate: {suc_rate:7.4f}%"
    else:
        mytag = f"False 0/{go.N_struc:d}"

    eng = go.min_energy
    t1 = int((time() - t0)/60)
    strs = "Final {:8s} [{:2d}]{:10s} ".format(name, sum(xtal.numMols), spg)
    strs += "{:3d}m {:2d} {:6.1f}".format(t1, N_torsion, wt)
    strs += "{:12.3f} {:20s} {:s}".format(eng, mytag, smile)
    print(strs)
