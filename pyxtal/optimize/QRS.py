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
            x.append(wp) # print(x)
            prev = seq
        rep = representation(x, smiles)
        xtal = rep.to_pyxtal(composition=comp)
        if not xtal.has_special_site() and len(xtal.check_short_distances(r=d_tol)) == 0:
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
        mlp (str): MACE | UMA | ANI
        skip_mlp (bool): whether or not use ani or not (default: True)
        output_mlp (bool): whether or not output the ANI relaxed structure (default: True)
        eng_cutoff (float): the cutoff energy for FF training
        E_max (float): maximum energy defined as an invalid structure
        verbose (bool): show more details
        use_mpi (bool): whether or not use mpi for parallel calculation
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
        mlp: str = "MACE",
        skip_mlp: bool = True,
        output_mlp: bool = True,
        factor: float = 1.1,
        eng_cutoff: float = 5.0,
        E_max: float = 1e10,
        verbose: bool = False,
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
            mlp,
            skip_mlp,
            output_mlp,
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

        if self.rank == 0:
            strs = self.full_str()
            self.logging.info(strs)
            print(strs)

    def full_str(self):
        s = str(self)
        s += "\nMethod    : Deterministic Quasi-Random Sampling"
        s += f"\nGeneration: {self.N_gen:4d}"
        s += f"\nPopulation: {self.N_pop:4d}"
        return s

    def _init_qrs_params(self):
        """
        Bootstrap hall_number, wp_bounds, ltype, and the Sobol sampler from a
        single randomly generated crystal using the fixed lattice.  Called once
        before the main generation loop when self.lattice is not None, so that
        generate_qrs_xtals can be used from gen 0 onwards.
        """
        from pyxtal import pyxtal
        from pyxtal.symmetry import Group
        tmp = pyxtal(molecular=True)
        #print(self.smiles, self.sg, self.lattice); import sys; sys.exit()
        smiles = [smi + ".smi" for smi in self.smiles]
        sg = self.sg[0]
        wp = Group(sg, use_hall=True)[0] if self.use_hall else Group(sg)[0]
        mult = len(wp)
        numIons = [int(c * mult) for c in self.composition]
        for _ in range(20):
            tmp.from_random(3, sg, smiles, numIons,
                            lattice=self.lattice, use_hall=self.use_hall,
                            sites=self.sites,   
                            force_pass=True)  # skip structures with bad distances
            if tmp.valid:
                break
        self.hall_number = sg
        self.ltype = self.lattice.ltype
        self.wp_bounds = [site.get_bounds() for site in tmp.mol_sites]
        len_reps = sum(len(b) for b in self.wp_bounds)
        self.sampler = qmc.Sobol(d=len_reps, scramble=False)

    def _run(self, pool=None):
        """
        The main code to run QRS prediction

        Returns:
            success_rate or None
        """
        self.ref_volumes = []
        success_rate = 0
        print(f"Rank {self.rank} starts QRS in {self.tag}")

        # When the lattice is fixed, bootstrap QRS parameters before gen 0 so
        # that the Sobol grid can be used for every generation (including gen 0).
        if self.rank == 0 and self.lattice is not None:
            self._init_qrs_params()

        for gen in range(self.N_gen):
            self.generation = gen
            cur_xtals = None

            if self.rank == 0:
                print(f"\nGeneration {gen:d} starts")
                self.logging.info(f"Generation {gen:d} starts")
                t0 = time()

                if self.lattice is not None:
                    # Fixed-lattice mode: use QRS sampling for every generation.
                    # self.sampler is advanced progressively so each generation
                    # draws the next block of the Sobol sequence without repeats.
                    cell = [self.hall_number] + self.lattice.encode()
                    cur_xtals = generate_qrs_xtals(cell,
                                                   self.wp_bounds,
                                                   self.N_pop,
                                                   self.smiles,
                                                   self.composition,
                                                   self.sampler)
                    strs = f"Cell parameters in Gen-{gen:d}: "
                    print(strs, cell, len(cur_xtals))
                else:
                    # Variable-lattice mode: gen 0 uses random crystals to
                    # calibrate the cell bounds; gen > 0 uses QRS.
                    cur_xtals = [(None, "Random")] * self.N_pop
                    if gen > 0:
                        cell = generate_qrs_cell(self.sampler,
                                                 self.cell_bounds,
                                                 self.ref_volumes[-1],
                                                 self.ltype)
                        cell = [self.hall_number] + cell
                        cur_xtals = generate_qrs_xtals(cell,
                                                       self.wp_bounds,
                                                       self.N_pop,
                                                       self.smiles,
                                                       self.composition)
                        strs = f"Cell parameters in Gen-{gen:d}: "
                        print(strs, cell, self.ref_volumes[-1], len(cur_xtals))


            # Broadcast
            if self.use_mpi: cur_xtals = self.comm.bcast(cur_xtals, root=0)

            # Local optimization
            gen_results = self.local_optimization(cur_xtals, qrs=True, pool=pool)

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
                if gen == 0 and self.lattice is None:
                    # Variable-lattice mode only: initialise cell and WP bounds
                    # from the best structure found in the random gen-0 batch.
                    best_xtal = cur_xtals[0][0]
                    self.cell_bounds = best_xtal.lattice.get_bounds(2.5, 25)
                    self.ltype = best_xtal.lattice.ltype
                    self.wp_bounds = [site.get_bounds() for site in best_xtal.mol_sites]
                    self.hall_number = best_xtal.group.hall_number
                    self.sampler = qmc.Sobol(d=len(self.cell_bounds), scramble=False)

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
    from pyxtal.representation import representation
    smiles = ["CC(=O)OC1=CC=CC=C1C(=O)O"]
    string = "82 11.43  6.49 11.19 83.31 1  0 0.77  0.57  0.53 48.55 24.31 145.9 -77.85 -4.40 170.9 0"
    rep3 = representation.from_string(string, smiles)
    xtal_ref = rep3.to_pyxtal()
    print(xtal_ref)

    # ── Reproducibility check when the lattice is fixed ───────────────────────
    # When a lattice is provided to QRS, the cell is frozen and only Wyckoff
    # positions are sampled via a Sobol sequence (scramble=False).
    # Because Sobol(scramble=False) is fully deterministic, two independent
    # invocations with the same cell / wp_bounds / N_pop must produce
    # exactly the same list of structures in the same order.
    print("\n--- Reproducibility check (fixed lattice) ---")
    cell = [xtal_ref.group.hall_number] + xtal_ref.lattice.encode()
    wp_bounds = [site.get_bounds() for site in xtal_ref.mol_sites]
    N_pop = 8
    comp = [1]

    xtals_a = generate_qrs_xtals(cell, wp_bounds, N_pop, smiles, comp,
                                  sampler_wp=qmc.Sobol(d=sum(len(b) for b in wp_bounds), scramble=False))
    xtals_b = generate_qrs_xtals(cell, wp_bounds, N_pop, smiles, comp,
                                  sampler_wp=qmc.Sobol(d=sum(len(b) for b in wp_bounds), scramble=False))

    assert len(xtals_a) == len(xtals_b), \
        f"Different number of valid structures: {len(xtals_a)} vs {len(xtals_b)}"

    all_match = True
    for i, ((xa, _), (xb, _)) in enumerate(zip(xtals_a, xtals_b)):
        ra = xa.get_1D_representation().to_string()
        rb = xb.get_1D_representation().to_string()
        if ra != rb:
            print(f"  Structure {i}: MISMATCH\n  run-1: {ra}\n  run-2: {rb}")
            all_match = False

    if all_match:
        print(f"  PASSED — both runs produced {len(xtals_a)} identical structures.")