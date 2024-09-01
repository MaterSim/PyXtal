"""
WFS sampler
"""

from __future__ import annotations

from time import time
from typing import TYPE_CHECKING

import numpy as np
from numpy.random import Generator
from pymatgen.analysis.structure_matcher import StructureMatcher

from pyxtal.optimize.base import GlobalOptimize

if TYPE_CHECKING:
    from pyxtal.lattice import Lattice
    from pyxtal.molecule import pyxtal_molecule


class WFS(GlobalOptimize):
    """
    Standard Population algorithm

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
        matcher : structurematcher from pymatgen
        early_quit (bool): if quit the program early when the target is found
        use_mpi (bool): if use mpi
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
        fracs: list | None = None,
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
        verbose: bool = False,
        random_state: int | None = None,
        max_time: float | None = None,
        matcher: StructureMatcher | None = None,
        early_quit: bool = False,
        check_stable: bool = False,
        use_mpi: bool = False,
    ):
        if isinstance(random_state, Generator):
            self.random_state = random_state.spawn(1)[0]
        else:
            self.random_state = np.random.default_rng(random_state)

        # POPULATION parameters:
        if fracs is None:
            fracs = [0.6, 0.4]
        self.N_gen = N_gen
        self.N_pop = N_pop
        self.fracs = np.array(fracs)
        self.verbose = verbose
        self.name = 'WFS'

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
            random_state,
            max_time,
            matcher,
            early_quit,
            check_stable,
            use_mpi,
        )

        # Setup the stats [N_gen, Npop, (E, matches)]
        self.stats = np.zeros([self.N_gen, self.N_pop, 2])
        self.stats[:, :, 0] = self.E_max

        if self.rank == 0:
            strs = self.full_str()
            self.logging.info(strs)
            print(strs)

    def full_str(self):
        s = str(self)
        s += "\nMethod    : Stochastic Width First Sampling"
        s += f"\nGeneration: {self.N_gen:4d}"
        s += f"\nPopulation: {self.N_pop:4d}"
        s += "\nFraction  : {:4.2f} {:4.2f}".format(*self.fracs)
        # The rest base information from now on
        return s

    def run_mpi(self):
        """
        The mpi code to run WFS prediction

        Returns:
            success_rate or None
        """
        # Related to the FF optimization
        N_added = 0
        success_rate = 0

        for gen in range(self.N_gen):
            print(f"Rank {self.rank} entering generation {gen}")

            current_xtals = None

            if self.rank == 0:
                print(f"\nGeneration {gen:d} starts")
                self.logging.info(f"Generation {gen:d} starts")
                self.generation = gen + 1
                t0 = time()

                # Initialize
                current_xtals = [(None, "Random")] * self.N_pop

                # WFS update
                if gen > 0:
                    N_pops = [int(self.N_pop * i) for i in self.fracs]
                    count = N_pops[0]
                    for _sub_pop in range(N_pops[1]):
                        id = self._selTournament(engs)
                        current_xtals[count] = (prev_xtals[id][0], "Mutation")
                        count += 1

            # broadcast
            current_xtals = self.comm.bcast(current_xtals, root=0)
            #print(f"Rank {self.rank} after broadcast: current_xtals = {current_xtals}")

            # Local optimization
            gen_results = self.local_optimization(gen, current_xtals)

            prev_xtals = None
            if self.rank == 0:
                # pass results, summary_and_ranking
                current_xtals, matches, engs = self.gen_summary(gen,
                                    t0, gen_results, current_xtals)

                # Save the reps for next move
                prev_xtals = current_xtals  # ; print(self.engs)

            # broadcast
            prev_xtals = self.comm.bcast(prev_xtals, root=0)

            # Update the FF parameters if necessary
            if self.ff_opt:
                N_added = self.update_ff_paramters(
                    current_xtals, engs, N_added)
            else:
                if self.rank == 0:
                    if self.ref_pmg is not None:
                        success_rate = self.success_count(gen,
                                                          current_xtals,
                                                          matches)
                        if self.early_termination(success_rate):
                            return success_rate

                    elif self.ref_pxrd is not None:
                        self.count_pxrd_match(gen,
                                              current_xtals,
                                              matches)

        return success_rate

    def run_serial(self):
        """
        The serial/multiprocess code to run WFS prediction

        Returns:
            success_rate or None
        """
        # Related to the FF optimization
        N_added = 0
        success_rate = 0

        for gen in range(self.N_gen):

            print(f"\nGeneration {gen:d} starts")
            self.logging.info(f"Generation {gen:d} starts")
            self.generation = gen + 1
            t0 = time()

            # Initialize
            current_xtals = [(None, "Random")] * self.N_pop

            # WFS update
            if gen > 0:
                N_pops = [int(self.N_pop * i) for i in self.fracs]
                count = N_pops[0]
                for _sub_pop in range(N_pops[1]):
                    id = self._selTournament(engs)
                    current_xtals[count] = (prev_xtals[id][0], "Mutation")
                    count += 1

            # Local optimization
            gen_results = self.local_optimization(gen, current_xtals)

            # pass results, summary_and_ranking
            current_xtals, matches, engs = self.gen_summary(gen,
                                t0, gen_results, current_xtals)

            # Save the reps for next move
            prev_xtals = current_xtals  # ; print(self.engs)

            # Update the FF parameters if necessary
            if self.ff_opt:
                N_added = self.update_ff_paramters(
                    current_xtals, engs, N_added)
            else:
                if self.ref_pmg is not None:
                    success_rate = self.success_count(gen,
                                                      current_xtals,
                                                      matches)

                    if self.early_termination(success_rate):
                        return success_rate

                elif self.ref_pxrd is not None:
                    self.count_pxrd_match(gen,
                                          current_xtals,
                                          matches)

        return success_rate

    def _selTournament(self, fitness, factor=0.35):
        """
        Select the best individual among *tournsize* randomly chosen
        individuals, *k* times. The list returned contains
        references to the input *individuals*.
        """
        IDs = self.random_state.choice(len(fitness), size=int(
            len(fitness) * factor), replace=False)
        min_fit = np.argmin(fitness[IDs])
        return IDs[min_fit]


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
    parser.add_argument("-n", "--ncpu", dest="ncpu", type=int,
                        default=1, help="cpu number, optional")
    parser.add_argument("--ffopt", action="store_true",
                        help="enable ff optimization")

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
    smile, wt, spg = row.mol_smi, row.mol_weight, row.space_group.replace(
        " ", "")
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
    N_torsion = xtal.get_num_torsions()

    # GO run
    t0 = time()
    go = WFS(
        smile,
        wdir,
        xtal.group.number,
        name.lower(),
        info=chm_info,
        ff_style="openff",  # 'gaff',
        ff_opt=ffopt,
        N_gen=gen,
        N_pop=pop,
        N_cpu=ncpu,
        cif="pyxtal.cif",
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
