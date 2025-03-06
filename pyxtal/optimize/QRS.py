"""
Global Optimizer base on Quasi-Random Sampling
"""
from __future__ import annotations
from time import time
from typing import TYPE_CHECKING

import numpy as np
from scipy.stats import qmc
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
            xtals.append((xtal, "QRandom"))
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
        use_mpi: bool = False,
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
        use_mpi: bool = False,
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
        s += "\nMethod    : Deterministic Quasi-Random Sampling"
        s += f"\nGeneration: {self.N_gen:4d}"
        s += f"\nPopulation: {self.N_pop:4d}"
        # The rest base information from now on
        return s

    def _run(self, pool=None):
        """
        The main code to run QRS prediction

        Returns:
            success_rate or None
        """
        self.ref_volumes = []
        success_rate = 0
        print(f"Rank {self.rank} starts QRS in {self.tag}")

        for gen in range(self.N_gen):
            self.generation = gen
            cur_xtals = None
            self.logging.info(f"Gen {gen} starts in Rank {self.rank}")

            if self.rank == 0:
                print(f"\nGeneration {gen:d} starts")
                self.logging.info(f"Generation {gen:d} starts")
                t0 = time()

                # Initialize
                cur_xtals = [(None, "Random")] * self.N_pop

                # QRS update
                if gen > 0:
                    if self.lattice is not None:
                        cell = [self.hall_number] + self.lattice.encode()
                        sampler = self.sampler
                    else:
                        cell = generate_qrs_cell(self.sampler,
                                                 self.cell_bounds,
                                                 self.ref_volumes[-1],
                                                 self.ltype)
                        cell = [self.hall_number] + cell
                        sampler = None

                    cur_xtals = generate_qrs_xtals(cell,
                                                   self.wp_bounds,
                                                   self.N_pop,
                                                   self.smiles,
                                                   self.composition,
                                                   sampler)
                    strs = f"Cell parameters in Gen-{gen:d}: "
                    print(strs, cell, self.ref_volumes[-1], len(cur_xtals))


            # Broadcast
            if self.use_mpi:
                cur_xtals = self.comm.bcast(cur_xtals, root=0)

            # Local optimization
            gen_results = self.local_optimization(cur_xtals, qrs=True, pool=pool)
            self.logging.info(f"Rank {self.rank} finishes local_opt")

            # Summary and Ranking
            quit = False
            if self.rank == 0:
                cur_xtals, matches, engs = self.gen_summary(t0,
                                        gen_results, cur_xtals)

                # update hist_best
                vols = []
                for id, (xtal, _) in enumerate(cur_xtals):
                    if xtal is not None:
                        vols.append(xtal.lattice.volume)

                # update best volume
                self.ref_volumes.append(np.array(vols).mean())
                if gen == 0:
                    best_xtal = cur_xtals[0][0]
                    self.cell_bounds = best_xtal.lattice.get_bounds(2.5, 25)
                    self.ltype = best_xtal.lattice.ltype
                    self.wp_bounds = [site.get_bounds() for site in best_xtal.mol_sites]
                    self.hall_number = best_xtal.group.hall_number
                    if self.lattice is not None:
                        len_reps = sum(len(bound) for bound in self.wp_bounds)
                    else:
                        len_reps = len(self.cell_bounds)
                    self.sampler = qmc.Sobol(d=len_reps, scramble=False)

                if self.ref_pmg is not None:
                    success_rate = self.success_count(cur_xtals, matches)

                    if self.early_termination(success_rate):
                        quit = True

                elif self.ref_pxrd is not None:
                    self.count_pxrd_match(cur_xtals, matches)

            if self.use_mpi:
                quit = self.comm.bcast(quit, root=0)
                self.comm.Barrier()
                self.logging.info(f"Gen {gen} Finish in Rank {self.rank}")

            if quit:
                self.logging.info(f"Early Termination in Rank {self.rank}")
                return success_rate

        return success_rate

if __name__ == "__main__":
    print("test")