"""
Global Optimizer
"""

from typing import List, Dict, Optional, Union

import numpy as np
from time import time
from random import sample
from concurrent.futures import ProcessPoolExecutor

from pyxtal.optimize.base import GlobalOptimize
from pyxtal.representation import representation
from pyxtal.molecule import pyxtal_molecule
from pyxtal.lattice import Lattice
from pyxtal.optimize.common import optimizer_single, optimizer_par


class GA(GlobalOptimize):
    """
    Standard Genetic algorithm

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
        sg: Union[int, list],
        tag: str,
        info: Optional[Dict[any, any]] = None,
        ff_opt: bool = False,
        ff_style: str = "openff",
        ff_parameters: str = "parameters.xml",
        reference_file: str = "references.xml",
        ref_criteria: Optional[Dict[any, any]] = None,
        N_gen: int = 10,
        N_pop: int = 10,
        fracs: list = [0.6, 0.4, 0.0],
        N_cpu: int = 1,
        cif: Optional[str] = None,
        block: Optional[List[any]] = None,
        num_block: Optional[List[any]] = None,
        composition: Optional[List[any]] = None,
        lattice: Optional[Lattice] = None,
        torsions: Optional[List[any]] = None,
        molecules: Optional[List[pyxtal_molecule]] = None,
        sites: Optional[List[any]] = None,
        use_hall: bool = False,
        skip_ani: bool = True,
        factor: float = 1.1,
        eng_cutoff: float = 5.0,
        E_max: float = 1e10,
        verbose: bool = False,
    ):
        # GA parameters:
        self.N_gen = N_gen
        self.N_pop = N_pop
        self.fracs = np.array(fracs)
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
        )

        print(self.full_str())

    def full_str(self):
        s = str(self)
        s += "\nMethod    : Genetic Algorithm"
        s += "\nGeneration: {:4d}".format(self.N_gen)
        s += "\nPopulation: {:4d}".format(self.N_pop)
        s += "\nFraction  : {:4.2f} {:4.2f} {:4.2f}".format(*self.fracs)
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
        if ref_pmg is not None:
            ref_pmg.remove_species("H")

        self.best_reps = []
        self.reps = []
        self.engs = []

        # Related to the FF optimization
        N_added = 0

        for gen in range(self.N_gen):
            print("\nGeneration {:d} starts".format(gen))
            self.generation = gen + 1
            t0 = time()

            # lists for structure information
            current_reps = [None] * self.N_pop
            if ref_pxrd is None:
                current_matches = [False] * self.N_pop
            else:
                current_matches = [0.0] * self.N_pop
            current_engs = [self.E_max] * self.N_pop
            current_xtals = [None] * self.N_pop
            current_tags = ["Random"] * self.N_pop

            # Set up the origin for new structures
            if gen > 0:
                N_pops = [int(self.N_pop * i) for i in self.fracs]
                count = N_pops[0]

                for sub_pop in range(N_pops[1]):
                    id = self._selTournament(engs)
                    xtal = prev_xtals[id]
                    current_tags[count] = "Mutation"
                    if xtal is None:
                        print(id)
                        print(len(engs))
                        print(engs)
                        print(len(prev_xtals))
                        print(prev_xtals)
                        raise ValueError("Problem in selection")
                    current_xtals[count] = xtal
                    count += 1

                # Not Working
                for sub_pop in range(N_pops[2]):
                    xtal1 = prev_xtals[self.selTournament(engs)]
                    xtal2 = prev_xtals[self.selTournament(engs)]
                    current_tags[count] = "Crossover"
                    current_xtals[count] = self._crossover(xtal1, xtal2)
                    count += 1

            # Local optimization
            args = [
                self.randomizer,
                self.optimizer,
                self.smiles,
                self.block,
                self.num_block,
                self.atom_info,
                self.workdir + "/" + "calc",
                self.sg,
                self.composition,
                self.lattice,
                self.torsions,
                self.molecules,
                self.sites,
                ref_pmg,
                self.matcher,
                ref_pxrd,
                self.use_hall,
                self.skip_ani,
            ]

            gen_results = [None] * len(current_xtals)
            if self.ncpu == 1:
                for pop in range(len(current_xtals)):
                    xtal = current_xtals[pop]
                    job_tag = self.tag + "-g" + str(gen) + "-p" + str(pop)
                    mutated = False if xtal is None else True
                    my_args = [xtal, pop, mutated, job_tag] + args
                    gen_results[pop] = optimizer_single(*tuple(my_args))

            else:
                # parallel process
                N_cycle = int(np.ceil(self.N_pop / self.ncpu))
                args_lists = []
                for i in range(self.ncpu):
                    id1 = i * N_cycle
                    id2 = min([id1 + N_cycle, len(current_xtals)])
                    # os.makedirs(folder, exist_ok=True)
                    ids = range(id1, id2)
                    job_tags = [
                        self.tag + "-g" + str(gen) + "-p" + str(id) for id in ids
                    ]
                    xtals = current_xtals[id1:id2]
                    mutates = [False if xtal is None else True for xtal in xtals]
                    my_args = [xtals, ids, mutates, job_tags] + args
                    args_lists.append(tuple(my_args))

                with ProcessPoolExecutor(max_workers=self.ncpu) as executor:
                    results = [executor.submit(optimizer_par, *p) for p in args_lists]
                    # loop each cpu
                    for result in results:
                        # loop each pop
                        for res in result.result():
                            (id, xtal, match) = res
                            gen_results[id] = (xtal, match)

            # Summary and Ranking
            for id, res in enumerate(gen_results):
                (xtal, match) = res

                if xtal is not None:
                    current_xtals[id] = xtal
                    rep = xtal.get_1D_representation()
                    current_reps[id] = rep.x
                    current_engs[id] = xtal.energy / sum(xtal.numMols)
                    current_matches[id] = match

                    # Don't write bad structure
                    if self.cif is not None:
                        if xtal.energy < 9999:
                            if self.verbose:
                                print("Add qualified structure", id, xtal.energy)
                            with open(self.workdir + "/" + self.cif, "a+") as f:
                                label = self.tag + "-g" + str(gen) + "-p" + str(id)
                                f.writelines(xtal.to_file(header=label))
                        # else:
                        #    print("Neglect bad structure", id, xtal.energy)
                    self.engs.append(xtal.energy / sum(xtal.numMols))
                    # print(output)

            print("Generation {:d} finishes".format(gen))  # ; import sys; sys.exit()
            t1 = time()

            # Apply Gaussian
            if ref_pxrd is None:
                engs = self._apply_gaussian(current_reps, current_engs)
            else:
                engs = self._apply_gaussian(
                    current_reps, -1 * np.array(current_matches)
                )

            # Store the best structures
            count = 0
            xtals = []
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
                    out = "{:3d} {:s} Top".format(gen, strs)
                    if ref_pxrd is not None:
                        out += " {:6.3f}".format(current_matches[id])
                    print(out)
                    count += 1
                if count == 3:
                    break

            t2 = time()
            gen_out = "Gen{:3d} time usage: ".format(gen)
            gen_out += "{:5.1f}[Calc] {:5.1f}[Proc]".format(t1 - t0, t2 - t1)
            print(gen_out)

            # Save the reps for next move
            prev_xtals = current_xtals  # ; print(self.engs)

            self.min_energy = np.min(np.array(self.engs))
            self.N_struc = len(self.engs)

            # Update the FF parameters if necessary
            #import sys; sys.exit()
            if self.ff_opt:
                N_max = min([int(self.N_pop * 0.6), 50])
                ids = np.argsort(engs)
                xtals = self.select_xtals(current_xtals, ids, N_max)
                print("Select Good structures for FF optimization", len(xtals))
                N_added = self.ff_optimization(xtals, N_added)


            elif match is not None:
                match = self.early_termination(current_xtals,
                                               current_matches,
                                               current_engs,
                                               current_tags,
                                               ref_pmg,
                                               ref_eng)
                print("Early termination")
                return

        return

    def _selTournament(self, fitness, factor=0.35):
        """
        Select the best individual among *tournsize* randomly chosen
        individuals, *k* times. The list returned contains
        references to the input *individuals*.
        """
        IDs = sample(set(range(len(fitness))), int(len(fitness) * factor))
        min_fit = np.argmin(fitness[IDs])
        return IDs[min_fit]

    def _crossover(self, x1, x2):
        """
        How to design this?
        """
        raise NotImplementedError


if __name__ == "__main__":
    from pyxtal.db import database
    import os
    import argparse

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
    parser.add_argument(
        "-n", "--ncpu", dest="ncpu", type=int, default=1, help="cpu number, optional"
    )
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
        if "charmm_info" in row.data.keys():
            # prepare charmm input
            chm_info = row.data["charmm_info"]
            prm = open(wdir + "/calc/pyxtal.prm", "w")
            prm.write(chm_info["prm"])
            prm.close()
            rtf = open(wdir + "/calc/pyxtal.rtf", "w")
            rtf.write(chm_info["rtf"])
            rtf.close()
        else:
            # Make sure we generate the initial guess from ambertools
            if os.path.exists("parameters.xml"):
                os.remove("parameters.xml")
    # load reference xtal
    pmg0 = xtal.to_pymatgen()
    if xtal.has_special_site():
        xtal = xtal.to_subgroup()

    # GA run
    t0 = time()
    ga = GA(
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

    match = ga.run(pmg0)
    if match is not None:
        eng = match["energy"]
        tmp = "{:}/{:}".format(match["rank"], ga.N_struc)
        mytag = "True[{:}] {:10s}".format(match["tag"], tmp)
        mytag += "{:5.2f}{:5.2f}".format(match["l_rms"], match["a_rms"])
    else:
        eng = ga.min_energy
        tmp = "0/{:}".format(ga.N_struc)
        mytag = "False   {:10s}".format(tmp)

    t1 = int((time() - t0) / 60)
    strs = "Final {:8s} {:8s} {:3d}m".format(name, spg, t1)
    strs += " {:3d}[{:2d}] {:6.1f}".format(ga.generation, ga.N_torsion, wt)
    strs += "{:12.3f} {:30s} {:s}".format(eng, mytag, smile)
    print(strs)
