"""
A base class for global optimization including
- GA
- PSO
- BasinHopping
- QRS
"""

import os
from random import sample
from time import time
from typing import Optional, Union

import numpy as np
import pymatgen.analysis.structure_matcher as sm
from ost.parameters import ForceFieldParameters, compute_r2, get_lmp_efs

from pyxtal.lattice import Lattice
from pyxtal.molecule import find_rotor_from_smile, pyxtal_molecule
from pyxtal.optimize.common import optimizer, randomizer
from pyxtal.representation import representation
from pyxtal.util import new_struc


class GlobalOptimize:
    """
    Base-class for all global optimization methods

    Args:
        smiles (str): smiles string
        workdir (str): path of working directory
        sg (int or list): space group number or list of spg numbers
        tag (string): job prefix
        ff_opt (bool): activate on the fly FF mode
        ff_style (str): automated force style (`gaff` or `openff`)
        ff_parameters (str or list): ff parameter xml file or list
        reference_file (str): path of reference xml data for FF training
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
    """

    def __init__(
        self,
        smiles: str,
        workdir: str,
        sg: Union[int, list[int]],
        tag: str,
        info: Optional[dict[any, any]] = None,
        ff_opt: bool = False,
        ff_style: str = "openff",
        ff_parameters: str = "parameters.xml",
        reference_file: str = "references.xml",
        ref_criteria: Optional[dict[any, any]] = None,
        N_cpu: int = 1,
        cif: Optional[str] = None,
        block: Optional[list[any]] = None,
        num_block: Optional[list[any]] = None,
        composition: Optional[list[any]] = None,
        lattice: Optional[Lattice] = None,
        torsions: Optional[list[any]] = None,
        molecules: Optional[list[pyxtal_molecule]] = None,
        sites: Optional[list[any]] = None,
        use_hall: bool = False,
        skip_ani: bool = True,
        factor: float = 1.1,
        eng_cutoff: float = 5.0,
        E_max: float = 1e10,
    ):
        # Molecular information
        self.smile = smiles
        self.smiles = self.smile.split(".")  # list
        self.torsions = torsions
        self.molecules = molecules
        self.block = block
        self.num_block = num_block
        self.composition = [1] * len(self.smiles) if composition is None else composition
        self.N_torsion = 0
        for smi, comp in zip(self.smiles, self.composition):
            self.N_torsion += len(find_rotor_from_smile(smi)) * int(max([comp, 1]))

        # Crystal information
        self.sg = [sg] if type(sg) == int or type(sg) == np.int64 else sg
        self.use_hall = use_hall
        self.factor = factor
        self.sites = sites
        self.lattice = lattice
        self.opt_lat = lattice is None
        self.ref_criteria = ref_criteria
        self.eng_cutoff = eng_cutoff

        # Generation and Optimization
        self.workdir = workdir
        self.ncpu = N_cpu
        self.skip_ani = skip_ani
        self.randomizer = randomizer
        self.optimizer = optimizer

        self.ff_opt = ff_opt

        if info is not None:
            self.atom_info = info
            self.parameters = None
            self.ff_opt = False
        else:
            self.ff_style = ff_style
            self.ff_parameters = ff_parameters
            self.reference_file = reference_file
            self.parameters = ForceFieldParameters(self.smiles, style=ff_style, f_coef=1.0, s_coef=1.0, ncpu=self.ncpu)

            # Preload two set for FF parameters 1 for opt and 2 for refinement
            if type(self.ff_parameters) == list:
                assert len(self.ff_parameters) == 2
                for para_file in self.ff_parameters:
                    if not os.path.exists(para_file):
                        raise RuntimeError("File not found", para_file)
                params0, dic = self.parameters.load_parameters(self.ff_parameters[0])
                if "ff_style" in dic:
                    assert dic["ff_style"] == self.ff_style
                # print(params0)
                params1, dic = self.parameters.load_parameters(self.ff_parameters[1])
                if "ff_style" in dic:
                    assert dic["ff_style"] == self.ff_style
                # print(params1)
                self.prepare_chm_info(params0, params1)
            else:
                if os.path.exists(self.ff_parameters):
                    print("Preload the existing FF parameters from", self.ff_parameters)
                    params0, _ = self.parameters.load_parameters(self.ff_parameters)
                else:
                    print(
                        "No FF parameter file exists, using the default setting",
                        ff_style,
                    )
                    params0 = self.parameters.params_init.copy()
                    self.parameters.export_parameters(self.workdir + "/" + self.ff_parameters, params0)

                self.prepare_chm_info(params0)

        # Structure matcher
        self.matcher = sm.StructureMatcher(ltol=0.3, stol=0.3, angle_tol=10)

        # I/O stuff
        self.E_max = E_max
        self.tag = tag
        self.cif = cif
        if cif is not None:
            with open(self.workdir + "/" + cif, "w") as f:
                f.writelines(str(self))
        # print(self)

    def __str__(self):
        s = "\n-------Global Crystal Structure Prediction------"
        s += f"\nsmile     : {self.smile:s}"
        s += f"\nZprime    : {self.composition!s:s}"
        s += f"\nN_torsion : {self.N_torsion:d}"
        s += f"\nsg        : {self.sg!s:s}"
        s += f"\nncpu      : {self.ncpu:d}"
        s += f"\ndiretory  : {self.workdir:s}"
        s += f"\nopt_lat   : {self.opt_lat!s:s}\n"
        if self.cif is not None:
            s += f"cif       : {self.cif:s}\n"
        if self.ff_opt:
            s += "forcefield: On-the-fly\n"
        else:
            s += "forcefield: Predefined\n"

        if self.parameters is not None:
            s += f"ff_style  : {self.ff_style:s}\n"
            if type(self.ff_parameters) == list:
                for para in self.ff_parameters:
                    s += f"ff_params : {para:s}\n"
            else:
                s += f"ff_params : {self.ff_parameters:s}\n"
            s += f"references: {self.reference_file:s}\n"
            s += str(self.parameters)

        return s

    def __repr__(self):
        return str(self)

    def new_struc(self, xtal, xtals):
        return new_struc(xtal, xtals)

    def select_xtals(self, ref_xtals, ids, N_max):
        """
        Select only unique structures
        """
        xtals = []
        for id in ids:
            xtal = ref_xtals[id]
            if xtal.energy <= self.E_max and self.new_struc(xtal, xtals):
                xtals.append(xtal)  # .to_ase(resort=False))
            if len(xtals) == N_max:
                break
        # xtals = [xtal.to_ase(resort=False) for xtal in xtals]
        return xtals

    def ff_optimization(self, xtals, N_added, N_min=50, dE=2.5, FMSE=2.5):
        """
        Optimize the current FF based on newly explored data

        Args:
            xtals: list of xtals from the current gen
            N_added (int): the number of structures that have been added
            N_min (int): minimum number of configurations to trigger the FF training
            dE (float): the cutoff energy value
            FMSE (float): the cutoff Force MSE value
        """
        from ost.utils import reset_lammps_cell

        numMols = [xtal.numMols for xtal in xtals]
        xtals = [xtal.to_ase(resort=False) for xtal in xtals]

        gen = self.generation
        # Initialize ff parameters and references
        params_opt, err_dict = self.parameters.load_parameters(self.ff_parameters)
        if os.path.exists(self.reference_file):
            ref_dics = self.parameters.load_references(self.reference_file)
            if self.ref_criteria is not None:
                ref_dics = self.parameters.clean_ref_dics(ref_dics, self.ref_criteria)
                ref_dics = self.parameters.cut_references_by_error(ref_dics, params_opt, dE=dE, FMSE=FMSE)
            # self.parameters.generate_report(ref_dics, params_opt)
        else:
            ref_dics = []

        # Add references
        print("Current number of reference structures", len(ref_dics))
        t0 = time()
        if len(ref_dics) > 100:  # no fit if ref_dics is large
            # Here we find the lowest engs and select only low-E struc
            ref_engs = [ref_dic["energy"] / ref_dic["replicate"] for ref_dic in ref_dics]
            ref_e2 = np.array(ref_engs).min()
            print("Min Reference Energy", ref_e2)

            if len(err_dict) == 0:
                # update the offset if necessary
                _, params_opt = self.parameters.optimize_offset(ref_dics, params_opt)
                results = self.parameters.evaluate_multi_references(ref_dics, params_opt, 1000, 1000)
                (ff_values, ref_values, rmse_values, r2_values) = results
                err_dict = {"rmse_values": rmse_values}

            _ref_dics = []
            rmse_values = err_dict["rmse_values"]
            lmp_in = self.parameters.ff.get_lammps_in()
            self.parameters.ase_templates = {}
            self.lmp_dat = {}
            ff_engs, ff_fors, ff_strs = [], [], []
            rf_engs, rf_fors, rf_strs = [], [], []
            for numMol, xtal in zip(numMols, xtals):
                struc = reset_lammps_cell(xtal)
                lmp_struc, lmp_dat = self.parameters.get_lmp_input_from_structure(struc, numMol, set_template=False)
                replicate = len(lmp_struc.atoms) / self.parameters.natoms_per_unit

                try:
                    e1, f1, s1 = get_lmp_efs(lmp_struc, lmp_in, lmp_dat)  # ; print('Debug KONTIQ', struc, e1)
                except:
                    e1 = self.E_max

                # filter very high energy structures
                if e1 < 1000:
                    struc.set_calculator(self.parameters.calculator)
                    e2 = struc.get_potential_energy()
                    e2 /= replicate
                    # Ignore very high energy structures
                    if e2 < ref_e2 + self.eng_cutoff:
                        e1 /= replicate
                        f2 = struc.get_forces()
                        s2 = struc.get_stress()
                        struc.set_calculator()
                        e_err = abs(e1 - e2 + params_opt[-1])
                        f_err = np.sqrt(((f1.flatten() - f2.flatten()) ** 2).mean())
                        s_err = np.sqrt(((s1 - s2) ** 2).mean())
                        e_check = e_err < 0.5 * rmse_values[0]
                        f_check = f_err < 1.0 * rmse_values[1]
                        s_check = s_err < 1.0 * rmse_values[2]
                        strs = f"Errors of csp structure in gen{gen:3d} "
                        strs += f"{e_err:.4f} {f_err:.4f} {s_err:.4f} "
                        strs += f"{e1:8.4f} {e2:8.4f}"
                        print(strs, e_check, f_check, s_check)

                        # avoid very unphysical structures
                        if e_err < 4.0 and f_err < 4.0:
                            ff_engs.append(e1 + params_opt[-1])
                            rf_engs.append(e2)
                            ff_fors.extend(f1.flatten())
                            rf_fors.extend(f2.flatten())
                            ff_strs.extend(s1.flatten())
                            rf_strs.extend(s2.flatten())

                            if False in [e_check, f_check, s_check]:
                                _ref_dic = {
                                    "structure": struc,
                                    "energy": e2 * replicate,
                                    "forces": f2,
                                    "stress": s2,
                                    "replicate": replicate,
                                    "options": [True, not f_check, True],
                                    "tag": "CSP",
                                    "numMols": numMol,
                                }
                                _ref_dics.append(_ref_dic)
                    else:
                        print("Ignore the structure due to high energy", e2)

            # QZ: Output FF performances MSE, R2 for the selected structures
            if len(_ref_dics) == 0:
                print("There is a serious problem in depositing high energy")
                raise ValueError("The program needs to stop here")

            if self.ref_criteria is not None:
                _ref_dics = self.parameters.clean_ref_dics(_ref_dics, self.ref_criteria)

            # print("Added {:d} new reference structures into training".format(len(_ref_dics)))
            print("FF performances")
            ff_engs = np.array(ff_engs)
            rf_engs = np.array(rf_engs)
            ff_fors = np.array(ff_fors)
            rf_fors = np.array(rf_fors)
            ff_strs = np.array(ff_strs)
            rf_strs = np.array(rf_strs)
            r2_engs = compute_r2(ff_engs, rf_engs)
            r2_fors = compute_r2(ff_fors, rf_fors)
            r2_strs = compute_r2(ff_strs, rf_strs)
            mse_engs = np.sqrt(np.mean((ff_engs - rf_engs) ** 2))
            mse_fors = np.sqrt(np.mean((ff_fors - rf_fors) ** 2))
            mse_strs = np.sqrt(np.mean((ff_strs - rf_strs) ** 2))

            print(f"R2   {r2_engs:8.4f} {r2_fors:8.4f} {r2_strs:8.4f}")
            print(f"RMSE {mse_engs:8.4f} {mse_fors:8.4f} {mse_strs:8.4f}")
            # self.parameters.generate_report(_ref_dics, params_opt)
            # import sys; sys.exit()
        else:
            # reduce the number of structures to save some time
            N_selected = min([N_min, self.ncpu])
            print("Create the reference data by augmentation", N_selected)
            if len(xtals) >= N_selected:
                ids = sample(list(range(len(xtals))), N_selected)
                xtals = [xtals[id] for id in ids]
                numMols = [numMols[id] for id in ids]

            _ref_dics = self.parameters.add_multi_references(
                xtals,
                numMols,
                augment=True,
                steps=20,  # 50,
                N_vibs=1,
                logfile="ase.log",
            )
            if len(ref_dics) == 0:
                _, params_opt = self.parameters.optimize_offset(_ref_dics, params_opt)
                self.parameters.generate_report(ref_dics, params_opt)
                # import sys; sys.exit()

        ref_dics.extend(_ref_dics)
        print(f"Add {len(_ref_dics):d} references in {(time() - t0) / 60:.2f} min")
        self.parameters.export_references(ref_dics, self.reference_file)

        # Optimize ff parameters if we get enough number of configurations
        N_added += len(_ref_dics)
        if N_added < N_min:
            print("Do not update ff, the current number of configurations is", N_added)
        else:
            t0 = time()
            _, params_opt = self.parameters.optimize_offset(ref_dics, params_opt)

            for data in [
                (["bond", "angle", "proper"], 50),
                (["proper", "vdW", "charge"], 50),
                (["bond", "angle", "proper", "vdW", "charge"], 50),
            ]:
                (terms, steps) = data

                # Actual optimization
                opt_dict = self.parameters.get_opt_dict(terms, None, params_opt)
                x, fun, values, it = self.parameters.optimize_global(
                    ref_dics, opt_dict, params_opt, obj="R2", t0=0.1, steps=25
                )

                params_opt = self.parameters.set_sub_parameters(values, terms, params_opt)
                opt_dict = self.parameters.get_opt_dict(terms, None, params_opt)
                x, fun, values, it = self.parameters.optimize_local(
                    ref_dics, opt_dict, params_opt, obj="R2", steps=steps
                )

                params_opt = self.parameters.set_sub_parameters(values, terms, params_opt)
                _, params_opt = self.parameters.optimize_offset(ref_dics, params_opt)

                # To add Early termination
            t = (time() - t0) / 60
            print(f"FF optimization {t:.2f} min ", fun)
            # Reset N_added to 0
            N_added = 0

        # Export FF performances
        if gen < 10:
            gen_prefix = "gen_00" + str(gen)
        elif gen < 100:
            gen_prefix = "gen_0" + str(gen)
        else:
            gen_prefix = "gen_" + str(gen)

        performance_fig = f"FF_performance_{gen_prefix:s}.png"
        errs = self.parameters.plot_ff_results(performance_fig, ref_dics, [params_opt], labels=gen_prefix)

        param_fig = f"parameters_{gen_prefix:s}.png"
        self.parameters.plot_ff_parameters(param_fig, [params_opt])

        # Save parameters
        self.parameters.export_parameters(self.ff_parameters, params_opt, errs[0])
        self.prepare_chm_info(params_opt)
        # self.parameters.generate_report(ref_dics, params_opt)

        return N_added

    def prepare_chm_info(self, params0, params1=None, suffix="calc/pyxtal"):
        """
        TODO: A base classs for optimization
        prepar_chm_info with the updated params.

        Args:
            params0 (array or list): FF parameters array
            params1 (array or list): FF parameters array
            suffix (str): suffix of the temporary file

        Returns:
            update the self.atom_info for next calculations
        """
        pwd = os.getcwd()
        os.chdir(self.workdir)

        # To remove the old pyxtal1 files
        if os.path.exists(suffix + ".rtf"):
            os.remove(suffix + ".rtf")
        if os.path.exists(suffix + ".prm"):
            os.remove(suffix + ".prm")

        ase_with_ff = self.parameters.get_ase_charmm(params0)
        ase_with_ff.write_charmmfiles(base=suffix)
        if params1 is not None:
            ase_with_ff = self.parameters.get_ase_charmm(params1)
            ase_with_ff.write_charmmfiles(base=suffix)
        os.chdir(pwd)

        # Info
        self.atom_info = ase_with_ff.get_atom_info()

    def early_termination(self, xtals, matches, engs, tags, ref_pmg, ref_eng):
        """
        Exit if a match is found
        """

        e, d1, d2 = 10000, 0, 0
        match_id = None
        if ref_pmg is not None:
            for id, match in enumerate(matches):
                if match and engs[id] < e:
                    xtal = xtals[id]
                    res = self._print_match(xtal, ref_pmg)
                    if res[0] is not None:
                        e, match_id, d1, d2 = engs[id], id, res[0], res[1]

            if match_id is not None:
                all_engs = np.sort(np.array(self.engs))
                rank = len(all_engs[all_engs < (e - 0.001)]) + 1
                tag = tags[match_id][0]
                done = False
                if rank / self.N_struc < 0.5:
                    done = True
                else:
                    if ref_eng is not None:
                        if e < ref_eng + 0.2:
                            done = True
                    else:
                        done = True
                if done:
                    return {
                        "energy": e,
                        "tag": tag,
                        "l_rms": d1,
                        "a_rms": d2,
                        "rank": rank,
                    }
                return None
            return None

        else:
            return None

    def get_label(self, i):
        if i < 10:
            folder = f"cpu00{i}"
        elif i < 100:
            folder = f"cpu0{i}"
        else:
            folder = f"cpu0{i}"
        return folder

    def _print_match(self, xtal, ref_pmg):
        """
        print the matched structure

        Args:
            rep: 1d rep
            eng: energy
            ref_pmg: reference pmg structure
        """

        rep0 = xtal.get_1D_representation()
        pmg_s1 = xtal.to_pymatgen()
        pmg_s1.remove_species("H")
        strs = rep0.to_string(eng=xtal.energy / sum(xtal.numMols))
        rmsd = self.matcher.get_rms_dist(ref_pmg, pmg_s1)
        if rmsd is not None:
            strs += f"{rmsd[0]:6.3f}{rmsd[1]:6.3f} Match Ref"
            print(strs)
            return rmsd[0], rmsd[1]
        else:
            return None, None

    def _apply_gaussian(self, reps, engs, h1=0.1, h2=0.1, w1=0.2, w2=3):
        """
        Apply Gaussian to discourage the sampling of already visited configs.
        Consider both lattice abc and torsion
        """
        from copy import deepcopy

        # check torsion
        N_id = 8
        engs_gau = deepcopy(engs)
        for i, rep in enumerate(reps):
            gau = 0
            if rep is not None and engs[i] < 9999:
                sg1, abc1 = rep[0][0], np.array(rep[0][1:])
                tor1 = np.zeros(self.N_torsion)
                count = 0
                for j in range(1, len(rep)):
                    if len(rep[j]) > N_id:  # for Cl-
                        tor1[count : count + len(rep[j]) - N_id - 1] = rep[j][N_id:-1]
                        count += len(rep[j]) - N_id

                for ref in self.best_reps:
                    sg2, abc2 = ref[0][0], np.array(ref[0][1:])
                    # Cell
                    g1 = 0
                    if sg1 == sg2:
                        diff1 = np.sum((abc1 - abc2) ** 2) / w1**2
                        g1 = h1 * np.exp(-0.5 * diff1)  # cell
                    # Torsion
                    g2 = 0
                    if len(tor1) > 0:
                        tor2 = np.zeros(self.N_torsion)
                        count = 0
                        for j in range(1, len(rep)):
                            if len(rep[j]) > N_id:  # for Cl-
                                tor2[count : count + len(ref[j]) - N_id - 1] = ref[j][N_id:-1]
                                count += len(ref[j]) - N_id

                        diff2 = np.sum((tor1 - tor2) ** 2) / w2**2
                        g2 = h2 * np.exp(-0.5 * diff2)  # torsion
                    gau += g1 + g2
                # if gau > 1e-2: print(sg1, diag1, abc1, tor1)
            # print("Gaussian", i, "eng", engs[i], "gau", gau)
            # import sys; sys.exit()
            engs_gau[i] += gau
        return np.array(engs_gau)

    def check_ref(self, reps=None, reference=None, filename="pyxtal.cif"):
        """
        check if ground state structure is found

        Args:
            reps: list of representations
            refernce: [pmg, eng]
            filename: filename
        """
        if os.path.exists(filename):
            os.remove(filename)

        if reference is not None:
            [pmg0, eng] = reference
            pmg0.remove_species("H")
            print("check if ground state structure is found")

            if reps is None:
                reps = np.array(self.reps)

            if eng is None:
                eng = np.min(reps[:, -1]) + 0.25

            reps = reps[reps[:, -1] < (eng + 0.1)]
            ids = np.argsort(reps[:, -1])
            reps = reps[ids]
            new_reps = []
            for rep in reps:
                eng1 = rep[-1]
                rep0 = representation(rep[:-1], self.smiles)
                xtal = rep0.to_pyxtal()
                pmg_s1 = xtal.to_pymatgen()
                pmg_s1.remove_species("H")
                new = True
                for ref in new_reps:
                    eng2 = ref[-1]
                    pmg_s2 = representation(rep[:-1], self.smiles).to_pyxtal().to_pymatgen()
                    pmg_s2.remove_species("H")
                    if abs(eng1 - eng2) < 1e-2 and sm.StructureMatcher().fit(pmg_s1, pmg_s2):
                        new = False
                        break
                if new:
                    new_reps.append(rep)
                    header = f"{len(new_reps):d}: {eng1:12.4f}"
                    xtal.to_file(filename, header=header, permission="a+")
                    strs = rep0.to_string(eng=eng1)
                    rmsd = self.matcher.get_rms_dist(pmg0, pmg_s1)
                    if rmsd is not None:
                        strs += f"{rmsd[0]:6.3f}{rmsd[1]:6.3f} True"
                        print(strs)
                        return True
                    else:
                        print(strs)
        return False


if __name__ == "__main__":
    print("test")
