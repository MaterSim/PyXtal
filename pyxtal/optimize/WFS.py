"""
Global optimization using the Stochastic Width First Sampling (WFS) algorithm
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
        composition: list of composition, (default is [1]*Num_mol)
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
        pre_opt: bool = False,
        check: bool = True,
    ):
        if isinstance(random_state, Generator):
            self.random_state = random_state.spawn(1)[0]
        else:
            self.random_state = np.random.default_rng(random_state)

        # POPULATION parameters:

        self.check = check
        self.N_gen = N_gen
        self.N_pop = N_pop
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
            pre_opt,
        )

        if fracs is None:
            if len(self.sg) > 4:
                fracs = [0.8, 0.2]
            else:
                fracs = [0.6, 0.4]
        self.fracs = np.array(fracs)

        if self.rank == 0:
            strs = self.full_str()
            self.logging.info(strs)
            print(strs)
        print(f"Rank {self.rank} finish initialization {self.tag}")

    def full_str(self):
        s = str(self)
        s += "\nMethod    : Stochastic Width First Sampling"
        s += f"\nGeneration: {self.N_gen:4d}"
        s += f"\nPopulation: {self.N_pop:4d}"
        s += "\nFraction  : {:4.2f} {:4.2f}".format(*self.fracs)
        return s

    def _run(self, pool=None):
        """
        The main code to run WFS prediction

        Returns:
            success_rate or None
        """

        # Related to the FF optimization
        success_rate = 0
        print(f"Rank {self.rank} starts WFS in {self.tag}")

        for gen in range(self.N_gen):
            self.generation = gen
            cur_xtals = None

            if self.rank == 0:
                print(f"\nGeneration {gen:d} starts")
                self.logging.info(f"Generation {gen} starts")
                t0 = time()

                # Initialize
                cur_xtals = [(None, "Random")] * self.N_pop

                # WFS update
                if gen > 0:
                    for i in range(self.N_pop):
                        if self.check_stable and not self.stats[gen-1][i, -1]:
                            # If the previously kept structure has no improvement,
                            # reset it to Mutation
                            if gen >= 2 and not self.stats[gen-2][i, -1] and \
                            self.stats[gen-1][i, 0] + 1e-3 > self.stats[gen-2][i, 0]:
                                cur_xtals[i] = (prev_xtals[i][0], "Mutation")
                            else:
                                cur_xtals[i] = (prev_xtals[i][0], "Kept")
                        else:
                            if self.random_state.random() > self.fracs[0]:
                                # Mutation
                                id = self._selTournament(engs)
                                cur_xtals[i] = (prev_xtals[id][0], "Mutation")

                        # If the space group is 1, it is a random structure
                        if len(self.sg) > 0 and cur_xtals[i][0] is not None and cur_xtals[i][1] == "Mutation":
                            if cur_xtals[i][0].group.number == 1:
                                cur_xtals[i] = (None, "Random")

            # broadcast
            if self.use_mpi: cur_xtals = self.comm.bcast(cur_xtals, root=0)

            # Local optimization
            gen_results = self.local_optimization(cur_xtals, pool=pool)

            prev_xtals = None
            if self.rank == 0:
                # pass results, summary_and_ranking
                cur_xtals, matches, engs = self.gen_summary(t0,
                                        gen_results, cur_xtals)

                # Save the reps for next move
                prev_xtals = cur_xtals

            # Broadcast
            if self.use_mpi:
                prev_xtals = self.comm.bcast(prev_xtals, root=0)
                self.logging.info(f"Gen {gen} bcast in Rank {self.rank}")

            # Update the FF parameters if necessary
            if self.ff_opt:
                self.export_references(cur_xtals, engs)
            else:
                quit = False
                if self.rank == 0:
                    if self.ref_pmg is not None:
                        success_rate = self.success_count(cur_xtals, matches)
                        if self.early_termination(success_rate):
                            quit = True

                    elif self.ref_pxrd is not None:
                        self.count_pxrd_match(cur_xtals, matches)

                # quit the loop
                if self.use_mpi:
                    quit = self.comm.bcast(quit, root=0)
                    self.comm.Barrier()
                    self.logging.info(f"Gen {gen} Finish in Rank {self.rank}")
                # Ensure that all ranks exit
                if quit:
                    self.logging.info(f"Early Termination in Rank {self.rank}")
                    return success_rate

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

    @classmethod
    def load(cls, filename):
        """
        Load the status of the WFS object
        """
        from pyxtal.optimize.base import load_xml

        # Define the parameter names in the same order as load_xml returns them
        param_names = [
        'smiles', 'workdir', 'sg', 'tag', 'info', 'ff_opt', 'ff_style',
        'ff_parameters', 'reference_file', 'ref_criteria', 'N_gen',
        'N_pop', 'N_cpu', 'fracs', 'cif', 'block', 'num_block',
        'composition', 'lattice', 'torsions', 'molecules', 'sites',
        'use_hall', 'skip_ani', 'factor', 'eng_cutoff', 'E_max',
        'verbose', 'random_state', 'max_time', 'matcher', 'early_quit',
        'check_stable', 'use_mpi', 'pre_opt']

        # Convert tuple to dictionary
        args = dict(zip(param_names, load_xml(filename)))
        return cls(**args)


if __name__ == "__main__":
    import argparse
    import os
    from pyxtal.db import database

    parser = argparse.ArgumentParser()
    parser.add_argument("--gen", dest="gen", type=int,
                        default=10, help="Generation size")
    parser.add_argument("--pop", dest="pop", type=int,
                        default=10, help="Population size")
    parser.add_argument("--ncpu", dest="ncpu", type=int,
                        default=1, help="cpu number")
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
    smile, wt = row.mol_smi, row.mol_weight
    spg = row.space_group.replace(" ", "")
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
    print(f"CSD {name} in Gen {go.generation}")

    if len(go.matches) > 0:
        best_rank = go.print_matches()
        mytag = f"True {best_rank}/{go.N_struc} Succ_rate: {suc_rate:.4f}%"
    else:
        mytag = f"False 0/{go.N_struc}"

    eng = go.min_energy
    t1 = int((time() - t0)/60)
    strs = "Final {:8s} [{:2d}]{:10s} ".format(name, sum(xtal.numMols), spg)
    strs += "{:3d}m {:2d} {:6.1f}".format(t1, N_torsion, wt)
    strs += "{:12.3f} {:20s} {:s}".format(eng, mytag, smile)
    print(strs)
