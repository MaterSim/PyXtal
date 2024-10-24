"""
mof_builder aims to generate crystal structure from the defined
building blocks (e.g., SiO2 with 4-coordined Si and 2-coordinated O)
1. Generate all possible wp combinations
2. For each wp combination,

    2.1 generate structure randomly
    2.2 optimize the geomtry
    2.3 parse the coordination
    2.4 save the qualified structure to ase.db

Todo:
1. Symmetry support for dSdX
2. add parallel for gulp optimiza
"""

# Standard Libraries
import os
import numpy as np
import random
from scipy.optimize import minimize
from pyxtal.lego.basinhopping import basinhopping
from pyxtal.lego.SO3 import SO3
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import Pool
from collections import deque
import gc

# Material science Libraries
from ase.io import write
from ase.db import connect
from pyxtal import pyxtal
from pyxtal.util import generate_wp_lib
from pyxtal.symmetry import Group
from pyxtal.lattice import Lattice
from pyxtal.db import database_topology

# Logging and debugging
import logging
from time import time
# np.set_printoptions(precision=3, suppress=True)

VECTORS = np.array([[x1, y1, z1] for x1 in range(-1, 2) for y1 in range(-1, 2) for z1 in range(-1, 2)])

def generate_wp_lib_par(spgs, composition, num_wp, num_fu, num_dof):
    """
    A wrapper to generate the wp list in parallel
    """

    my_spgs, wp_libs = [], []
    for spg in spgs:
        wp_lib = generate_wp_lib([spg], composition, num_wp, num_fu, num_dof),
        if len(wp_lib) > 0:
            wp_libs.append(wp_lib)
            my_spgs.append(spg)
    return (my_spgs, wp_libs)

def generate_xtal_par(wp_libs, niter, dim, elements, calculator, ref_environments,
                      criteria, T, N_max, early_quit):
    """
    A wrapper to call generate_xtal function in parallel
    """
    xtals, sims = [], []
    for wp_lib in wp_libs:
        (number, spg, wps, dof) = wp_lib
        xtal, sim = generate_xtal(dim, spg, wps, niter*dof, elements, calculator,
                                  ref_environments, criteria, T, N_max, early_quit)
        if xtal is not None:
            xtals.append(xtal)
            sims.append(sim)

    return (xtals, sims)

def minimize_from_x_par(*args):
    """
    A wrapper to call minimize_from_x function in parallel
    """
    dim, wp_libs, elements, calculator, ref_environments, opt_type, T, niter, early_quit, minimizers = args[0]
    xtals = []
    xs = []
    for wp_lib in wp_libs:
        (x, spg, wps) = wp_lib
        res = minimize_from_x(x, dim, spg, wps, elements, calculator,
                              ref_environments, T, niter,
                              early_quit, opt_type, minimizers,
                              filename=None)
        if res is not None:
            xtals.append(res[0])
            xs.append(res[1])
    return xtals, xs


def generate_xtal(dim, spg, wps, niter, elements, calculator,
                  ref_environments, criteria, T, N_max, early_quit,
                  dump=False, random_state=None, verbose=False):
    """
    Generate a crystal with the desired local environment

    Args:
        dim (int): 1, 2, 3
        spg: pyxtal.symmetry.Group object
        wps: list of wps for the disired crystal (e.g., [wp1, wp2])
        ref_env: reference enviroment
        f: callable function to compute env
        n_iter (int):
        T (float): for basinhopping

    Returns:
        xtal and its similarity
    """

    # Here we start to optimize the xtal based on similarity

    print("\n", dim, spg, wps, T, N_max, early_quit)
    count = 0
    while True:

        filename = 'opt.txt' if dump else None
        result = minimize_from_x(None, dim, spg, wps, elements, calculator,
                                 ref_environments, T, niter, early_quit,
                                 'global', filename=filename,
                                 random_state=random_state)

        if result is not None:
            (xtal, xs) = result
            if xtal.check_validity(criteria, verbose=verbose):
                x = xtal.get_1d_rep_x()
                sim1 = calculate_S(x, xtal, ref_environments, calculator)
                print(xtal.get_xtal_string({'sim': sim1}))
            else:
                xtal = None
                sim1 = None
            return xtal, sim1

        count += 1
        if count == N_max:
            break
    return None, None


def minimize_from_x(x, dim, spg, wps, elements, calculator, ref_environments,
                    T=0.2, niter=20, early_quit=0.02, opt_type='local',
                    minimizers=[('Nelder-Mead', 100), ('L-BFGS-B', 100)],
                    filename='local_opt_data.txt', random_state=None,
                    #derivative=True):
                    derivative=False):
    """
    Generate xtal from the 1d representation

    Args:
        x: list of 1D array
        spg (int): space group number (1-230)
        wps (string): e.g. [['4a', '8b']]
        elements (string): e.g., ['Si', 'O']
    """
    if derivative:
        jac = calculate_dSdx
    else:
        jac = None

    g, wps, dof = get_input_from_letters(spg, wps, dim)
    l_type = g.lattice_type
    sites_wp = []
    sites = []
    numIons = []
    ref_envs = None
    for i, wp in enumerate(wps):
        site = []
        numIon = 0
        for w in wp:
            sites.append((elements[i], w))  # .get_label())
            site.append(w.get_label())
            numIon += w.multiplicity
            if ref_envs is None:
                ref_envs = ref_environments[i]
            else:
                ref_envs = np.vstack((ref_envs, ref_environments[i]))
        sites_wp.append(site)
        numIons.append(numIon)

    if len(ref_envs.shape) == 1:
        ref_envs = ref_envs.reshape((1, len(ref_envs)))

    xtal = pyxtal()
    if x is None:
        count = 0
        while True:
            count += 1
            try:
                xtal.from_random(dim, g, elements, numIons,
                                 sites=sites_wp, factor=1.0,
                                 random_state=random_state)
            except RuntimeError:
                print(g.number, numIons, sites)
                print("Trouble in generating random xtals from pyxtal")
            if xtal.valid:
                atoms = xtal.to_ase(resort=False, add_vaccum=False)
                try:
                    des = calculator.calculate(atoms)['x']
                except:
                    if filename is not None:
                        print('Not a good structure, skip')
                    continue
                x = xtal.get_1d_rep_x()
                break
            elif count == 5:
                return None
    else:
        sites = []
        for ele, _wps in zip(elements, wps):
            for wp in _wps:
                sites.append((ele, wp))
        try:
            xtal.from_1d_rep(x, sites, dim=dim)
        except:
            return None

    x0 = np.array(x.copy())
    # Extract variables, call from Pyxtal
    [N_abc, N_ang] = Lattice.get_dofs(xtal.lattice.ltype)
    rep = xtal.get_1D_representation()
    xyzs = rep.x[1:]

    # Set constraints and minimization
    bounds = [(1.5, 50)] * (N_abc) + [(30, 150)] * (N_ang)

    # Special treatment in case the random lattice is small
    for i in range(N_abc):
        if x[i] < 1.5:
            x[i] = 1.5
        if x[i] > 50.0:
            x[i] = 50.0

    for i in range(N_abc, N_abc + N_ang):
        if x[i] < 30.0:
            x[i] = 30.0
        if x[i] > 150.0:
            x[i] = 150.0

    for xyz in xyzs:
        if len(xyz) > 2:
            bounds += [(0.0, 1.0)] * len(xyz[2:])

    if len(x) != len(bounds):
        print('debug before min', xtal, x, bounds, len(x), len(bounds))

    sim0 = calculate_S(x, xtal, ref_envs, calculator)
    if filename is not None:
        with open(filename, 'a+') as f0:
            f0.write('\nSpace Group: {:d}\n'.format(xtal.group.number))
            for element, numIon, site in zip(elements, numIons, sites_wp):
                strs = 'Element: {:2s} {:4d} '.format(element, numIon)
                for s in site:
                    strs += '{:s} '.format(s)
                strs += '\n'
                f0.write(strs)
            # Initial value
            strs = 'Init: {:9.3f} '.format(sim0)
            for x0 in x:
                strs += '{:8.4f} '.format(x0)
            strs += '\n'
            print(strs)
            f0.write(strs)

    # Run local minimization
    if opt_type == 'local':
        # set call back function for debugging
        def print_local_fun(x):
            f = calculate_S(x, xtal, ref_envs, calculator)
            print("{:.4f} ".format(f), x)
            if filename is not None:
                with open(filename, 'a+') as f0:
                    strs = 'Iter: {:9.3f} '.format(f)
                    for x0 in x[:3]:
                        strs += '{:8.4f} '.format(x0)
                    strs += '\n'
                    f0.write(strs)
        callback = print_local_fun if filename is not None else None

        for minimizer in minimizers:
            (method, step) = minimizer
            #print("Starting", xtal.lattice, method, step)
            if len(x) != len(bounds):
                print('debug min', xtal, x, bounds, len(x), len(bounds))
            res = minimize(calculate_S, x,
                           method=method,
                           args=(xtal, ref_envs, calculator),
                           jac=None if method=='Nelder-Mead' else jac,
                           bounds=bounds,
                           options={'maxiter': step}, #'disp': True},
                           callback=callback)
            x = res.x
            if xtal.lattice is None:
                return None
            #import sys; sys.exit()

        if filename is not None:
            with open(filename, 'a+') as f0:
                f0.write('END\n')
    else:
        # set call back function for debugging
        def print_fun_local(x):
            f = calculate_S(x, xtal, ref_envs, calculator)
            # print("{:.4f} ".format(f), x)
            if f < early_quit:
                return True
            else:
                return None
        callback = print_fun_local  # if verbose else None

        minimizer_kwargs = {'method': ['Nelder-Mead', 'l-bfgs-b',
                                       'Nelder-Mead', 'l-bfgs-b'],
                            'args': (xtal, ref_envs, calculator),
                            'bounds': bounds,
                            'callback': callback,
                            'options': {'maxiter': 100,
                                        'fatol': 1e-6,
                                        'ftol': 1e-6}}

        bounded_step = RandomDispBounds(np.array([b[0] for b in bounds]),
                                        np.array([b[1] for b in bounds]),
                                        id1=N_abc + N_ang,
                                        id2=N_abc)

        # set call back function for debugging
        def print_fun(x, f, accepted):
            if filename is not None:
                print("minimum {:.4f}[{:.4f}] accepted {:d} ".format(
                    f, early_quit, int(accepted)), x[:N_abc])
            if f < early_quit:
                # print("Return True", True is not None)
                return True
            else:
                return None
        callback = print_fun  # if verbose else None

        # Run BH optimization
        res = basinhopping(calculate_S, x, T=T,
                           minimizer_kwargs=minimizer_kwargs,
                           niter=niter,
                           take_step=bounded_step,
                           callback=callback)
        if xtal.lattice is None:
            return None

    # Extract the optimized xtal
    xtal = pyxtal()
    try:
        xtal.from_1d_rep(res.x, sites, dim=dim)
        return xtal, (x0, res.x)
    except:
        return None


def calculate_dSdx(x, xtal, des_ref, f, eps=1e-4, symmetry=True, verbose=False):
    """
    Compute dSdx via the chain rule of dSdr*drdx.
    This routine will serve as the jacobian for scipy.minimize

    Args:
        x (float): input variable array to describe a crystal structure
        xtal (instance): pyxtal object
        des_ref (array): reference environment
        f (callable): function to compute the environment
        symmetry (bool): whether or not turn on symmetry
        verbose (bool): output more information
    """
    xtal.update_from_1d_rep(x)
    ids = [0]
    weights = []
    for site in xtal.atom_sites:
        ids.append(site.wp.multiplicity + ids[-1])
        weights.append(site.wp.multiplicity)
    ids = ids[:-1]
    weights = np.array(weights, dtype=float)

    atoms = xtal.to_ase()
    dPdr, P = f.compute_dpdr_5d(atoms, ids)
    dPdr = np.einsum('i, ijklm -> ijklm', weights, dPdr)

    # Compute dSdr [N, M] [N, N, M, 3, 27] => [N, 3, 27]
    # print(P.shape, des_ref.shape, dPdr.shape)
    dSdr = np.einsum("ik, ijklm -> jlm", 2*(P - des_ref), dPdr)

    # Get supercell positions
    ref_pos = np.repeat(atoms.positions[:, :, np.newaxis], 27, axis=2)
    cell_shifts = np.dot(VECTORS, atoms.cell)
    ref_pos += cell_shifts.T[np.newaxis, :, :]

    # Compute drdx via numerical func
    drdx = np.zeros([len(atoms), 3, 27, len(x)])

    xtal0 = xtal.copy()
    for i in range(len(x)):
        x0 = x.copy()
        x0[i] += eps
        # Maybe expensive Reduce calls, just need positions
        xtal0.update_from_1d_rep(x0)
        atoms = xtal0.to_ase()

        # Get supercell positions
        pos = np.repeat(atoms.positions[:, :, np.newaxis], 27, axis=2)
        cell_shifts = np.dot(VECTORS, atoms.cell)
        pos += cell_shifts.T[np.newaxis, :, :]
        drdx[:, :, :, i] += (pos - ref_pos)/eps

    # [N, 3, 27] [N, 3, 27, H] => H
    dSdx = np.einsum("ijk, ijkl -> l", dSdr, drdx)
    return dSdx

def calculate_S(x, xtal, des_ref, f, ref_CN=None, verbose=False):
    """
    Optimization function used for structure generation

    Args:
        x (float): input variable array to describe a crystal structure
        xtal (instance): pyxtal object
        des_ref: reference environment
        f (callable): function to compute the enivronment
        verbose (boolean): output more information
    """
    if xtal.lattice is None:
        des = np.zeros(des_ref.shape)
        weights = 1
    else:
        xtal.update_from_1d_rep(x)  # ; print(xtal)
        if xtal.lattice is not None:
            ids = [0]
            weights = []
            for site in xtal.atom_sites:
                ids.append(site.wp.multiplicity + ids[-1])
                weights.append(site.wp.multiplicity)
            ids = ids[:-1]
            weights = np.array(weights, dtype=float)
            atoms = xtal.to_ase(resort=False)
            #des = f.calculate(atoms)['x'][ids]
            des, CNs = f.compute_p(atoms, ids, return_CN=True)
        else:
            des = np.zeros(des_ref.shape)
            weights = 1

    sim = np.sum((des-des_ref)**2, axis=1)  # /des.shape[1]
    if ref_CN is not None:
        coefs = CNs[:, 1] - ref_CN
        coefs[coefs < 0] = 0
        penalty = coefs * (f.rcut + 0.01 - CNs[:, 0])
        if penalty.sum() > 0: print('debug', penalty.sum(), x[:3])
        sim += 5 * coefs * penalty
    obj = np.sum(sim * weights)
    #print(des-des_ref); print(sim); print(obj)#; import sys; sys.exit()
    return obj


def get_input_from_letters(spg, sites, dim=3):
    """
    Prepare the input from letters

    Args:
        spg (int): 1-230
        sites (list): [['4a'], ['4b']] or [[0], [1]]
        dim (int): 0, 1, 2, 3

    Returns:
        space group instance + wp_list + dof
    """
    g = Group(spg, dim=dim)
    wp_combo = []
    dof = g.get_lattice_dof()
    for site in sites:
        wp_specie = []
        for s in site:
            if type(s) in [int, np.int64]:
                wp = g[s]
            else:
                wp = g.get_wyckoff_position(s)
            wp_specie.append(wp)
            dof += wp.get_dof()
        wp_combo.append(wp_specie)
    return g, wp_combo, dof


def create_trajectory(dumpfile, trjfile, modes=['Init', 'Iter'], dim=3):
    """
    create the ase trajectory from the dump file.
    This is used to analyze the performance of structure optimization
    """
    keyword_start = 'Space Group'  # \n'
    keyword_end = 'END'  # \n'

    with open(dumpfile, 'r') as f:
        lines = f.readlines()
        starts, ends = [], []
        for i in range(len(lines)):
            l = lines[i].lstrip()
            if l.startswith(keyword_start):
                starts.append(i)
            elif l.startswith(keyword_end):
                ends.append(i)
    # print(len(starts), len(ends))
    assert (len(starts) == len(ends))

    atoms = []
    for i, j in zip(starts, ends):
        # parse each run
        xtal_strs = lines[i:j]
        spg = int(xtal_strs[0].split(':')[1])
        l_type = Group(spg).lattice_type
        elements = []
        numIons = []
        wps = []
        for count, tmp in enumerate(xtal_strs[1:]):
            if tmp.startswith('Element'):
                tmp1 = tmp.split(':')
                tmp2 = tmp1[1].split()
                elements.append(tmp2[0])
                numIons.append(int(tmp2[1]))
                wps.append(tmp2[2:])
            else:
                line_struc = count + 1
                break

        # print(spg, elements, numIons, wps)
        g, wps, dof = get_input_from_letters(spg, wps, dim)
        for line_number in range(line_struc, len(xtal_strs)):
            tmp0 = xtal_strs[line_number].split(':')
            tag = tmp0[0]
            if tag in modes:
                tmp = tmp0[1].split()
                sim = float(tmp[0])
                x = [float(t) for t in tmp[1:]]

                # QZ: to fix
                xtal = pyxtal()
                xtal.from_1d_rep(x, wps, numIons, l_type, elements, dim)

                struc = xtal.to_ase()
                struc.info = {'time': line_number-line_struc, 'fun': sim}
                atoms.append(struc)
                # print(xtal.lattice)
        write(filename=trjfile, images=atoms, format='extxyz')


class mof_builder(object):
    """
    Class for generating MOF like structures

    Args:
        elements (str): e.g. ['Si', 'O']
        composition (int): e.g. [1, 2]
        dim (int): e.g., 0, 1, 2, 3
        db_file (str): default is 'mof.db'
        log_file (str): default is 'mof.log'

    Examples
    --------

    To create a new structure instance

    >>> from mof_builder import mof_builder
    >>> builder = mof_builder(['P', 'O', 'N'], [1, 1, 1], db_file='PON.db')
    """

    def __init__(self, elements, composition, dim=3, prefix='mof',
                 db_file=None, log_file=None, rank=0, verbose=False):

        self.rank = rank
        self.prefix = f"{prefix}-{rank}"
        # Define the chemical system
        self.dim = dim
        self.elements = elements
        self.composition = composition

        # Initialize neccessary functions and attributes
        self.calculator = None       # will be a callable function
        self.ref_environments = None  # will be a numpy array
        self.criteria = {}           # will be a dictionary
        self.verbose = verbose

        # Define the I/O
        logging.getLogger().handlers.clear()
        if log_file is not None:
            self.log_file = log_file
        else:
            self.log_file = self.prefix + '.log'
        logging.basicConfig(format="%(asctime)s| %(message)s",
                            filename=self.log_file,
                            level=logging.INFO)

        self.logging = logging
        if db_file is not None:
            self.db_file = db_file
        else:
            self.db_file = self.prefix + '.db'
        self.db = database_topology(self.db_file, log_file=self.log_file)

    def __str__(self):

        s = "\n------MOF Builder------"
        s += "\nSystem: "
        for element, comp in zip(self.elements, self.composition):
            s += "{:s}{:d} ".format(element, comp)
        s += "\nDatabase: {}".format(self.db_file)
        s += "\nLog_file: {}".format(self.log_file)
        if self.calculator is not None:
            s += "\nDescriptor: {}".format(self.calculator)
        if self.ref_environments is not None:
            (d1, d2) = self.ref_environments.shape
            s += "Reference enviroments ({:}, {:})".format(d1, d2)
        if len(self.criteria.keys()) > 0:
            for key in self.criteria.keys():
                s += '\nCriterion_{:}: {:}'.format(key, self.criteria[key])
        return s

    def __repr__(self):
        return str(self)

    def print_memory_usage(self):
        import psutil
        process = psutil.Process(os.getpid())
        mem = process.memory_info().rss / 1024 ** 2
        self.logging.info(f"Rank {self.rank} memory: {mem:.1f} MB")
        print(f"Rank {self.rank} memory: {mem:.1f} MB")

    def set_descriptor_calculator(self, dtype='SO3', mykwargs={}):
        """
        Set up the calculator for descriptor computation.
        Here we mostly use the pyxtal_ff module

        Arg:
            dytpe (str): only SO3 is suppoted now
        """
        if dtype == 'SO3':
            kwargs = {'lmax': 4,
                      'nmax': 2,
                      'rcut': 2.2,
                      'alpha': 1.5,
                      'weight_on': True,
                      }
            kwargs.update(mykwargs)

            self.calculator = SO3(**kwargs)

    def set_reference_enviroments(self, cif_file, substitute=None):
        """
        Get the reference enviroments

        Args:
            cif_file (str): cif structure
            substitute (dict): substitution directory
        """

        if self.calculator is None:
            raise RuntimeError(
                "Must call set_descriptor_calculator in advance")

        xtal = pyxtal()
        xtal.from_seed(cif_file)
        if substitute is not None:
            xtal.substitute(substitute)  # ; print(xtal)
        xtal.resort_species(self.elements)

        ids = [0] * len(self.elements)
        count = 0
        for site in xtal.atom_sites:
            for i, element in enumerate(self.elements):
                if element == site.specie:
                    ids[i] = count
                    break
            count += site.multiplicity
        if self.verbose:
            print("ids from Reference xtal", ids)
        atoms = xtal.to_ase(resort=False)
        self.ref_environments = self.calculator.compute_p(atoms, ids)
        if self.verbose:
            print(self.ref_environments)
        self.ref_xtal = xtal

    def set_criteria(self, CN=None, dimension=None, min_density=None, exclude_ii=False):
        """
        define the criteria to check if a structure is good

        Args:
            CN (int): coordination number
            dimension (int): target dimensionality (e.g. we want the 3D structure)
            min_density (float): minimum density
            exclude_ii (bool): allow the bond between same element
        """

        if CN is not None:
            self.criteria["CN"] = CN
            if 'cutoff' in CN.keys():
                self.criteria['cutoff'] = CN['cutoff']
            else:
                self.criteria['cutoff'] = None
        if dimension is not None:
            self.criteria["Dimension"] = dimension
        if min_density is not None:
            self.criteria["MIN_Density"] = min_density
        self.criteria['exclude_ii'] = exclude_ii

    def get_input_from_letters(self, spg, wps):
        """
        A short cut to get (spg, wps, dof) from the get_input functions

        Args:
            spg (int): space group number 1-230
            wps (list): e.g. [['4a', '4b']]
        """
        return get_input_from_letters(spg, wps, self.dim)

    def get_input_from_ref_xtal(self, xtal, substitute=None):
        """
        Generate the input from a given pyxtal

        Args:
            xtal: pyxtal object
            substitute: a dictionary to describe chemical substitution
        """

        # Make sure it is a pyxtal object
        if type(xtal) == str:
            c = pyxtal()
            c.from_seed(xtal)
            xtal = c
        if substitute is not None:
            xtal.substitute(substitute)

        g = xtal.group
        sites = [[] for _ in self.elements]
        dof = xtal.lattice.dof
        for s in xtal.atom_sites:
            for i, specie in enumerate(self.elements):
                if s.specie == specie:
                    # wp_combo[i].append(s.wp)
                    sites[i].append(s.wp.get_label())
                    break
            dof += s.wp.get_dof()

        return g.number, sites, dof

    def get_similarity(self, xtal):
        """
        Compute the similrity for a given xtal
        QZ: call calculate_S()

        Args:
            xtal: pyxtal object

        Returns:
            the similarity value .w.r.t the reference environment
        """
        x = xtal.get_1d_rep_x()
        return calculate_S(x, xtal, self.ref_environments,
                           self.calculator)

    def process_xtals(self, xtals, xs, add_db, symmetrize):
        # Now process each of the results
        valid_xtals = []
        count = 0
        for xtal, _xs in zip(xtals, xs):
            status = xtal.check_validity(self.criteria, verbose=self.verbose)
            if status:
                valid_xtals.append(xtal)
                sim1 = self.get_similarity(xtal)
                if symmetrize:
                    pmg = xtal.to_pymatgen()
                    xtal = pyxtal()
                    xtal.from_seed(pmg)
                if add_db:
                    self.process_xtal(xtal, [0, sim1], count, xs=_xs)
                    count += 1
                else:
                    dicts = {'sim': "{:6.3f}".format(sim1)}
                    print(xtal.get_xtal_string(dicts))
        return valid_xtals

    def optimize_xtals(self, xtals, ncpu=1, opt_type='local',
                       T=0.2, niter=20, early_quit=0.02,
                       add_db=True, symmetrize=False,
                       minimizers=[('Nelder-Mead', 100), ('L-BFGS-B', 100)],
                       ):
        """
        Perform optimization for each structure

        Args:
            xtals: list of xtals
            ncpu (int):

        """
        args = (opt_type, T, niter, early_quit, add_db, symmetrize, minimizers)
        if ncpu == 1:
            valid_xtals = self.optimize_xtals_serial(xtals, args)
        else:
            valid_xtals = self.optimize_xtals_mproc(xtals, ncpu, args)
        return valid_xtals

    def optimize_xtals_serial(self, xtals, args):
        """
        Optimization in serial mode.

        Args:
            xtals: list of xtals
            args: (opt_type, T, n_iter, early_quit, add_db, symmetrize, minimizers)
        """
        # (opt_type, T, n_iter, early_quit, add_db, symmetrize, minimizers) = args
        xtals_opt = []
        for i, xtal in enumerate(xtals):
            xtal, sim, _xs = self.optimize_xtal(xtal, i, *args)
            if xtal is not None:
                xtals_opt.append(xtal)
        return xtals_opt

    def optimize_xtals_mproc(self, xtals, ncpu, args):
        """
        Optimization in multiprocess mode.

        Args:
            xtals: list of xtals
            ncpu (int): number of parallel python processes
            args: (opt_type, T, n_iter, early_quit, add_db, symmetrize, minimizers)
        """

        pool = Pool(processes=ncpu)
        (opt_type, T, niter, early_quit, add_db, symmetrize, minimizers) = args
        xtals_opt = deque()

        # Split the input structures to minibatches
        N_batches = 50 * ncpu
        for _i, i in enumerate(range(0, len(xtals), N_batches)):
            start, end = i, min([i+N_batches, len(xtals)])
            ids = list(range(start, end))
            print(f"Rank {self.rank} minibatch {start} {end}")
            self.print_memory_usage()

            def generate_args():
                """
                A generator to yield argument lists for minimize_from_x_par.
                """
                for j in range(ncpu):
                    _ids = ids[j::ncpu]
                    wp_libs = []
                    for id in _ids:
                        xtal = xtals[id]
                        x = xtal.get_1d_rep_x()
                        spg, wps, _ = self.get_input_from_ref_xtal(xtal)
                        wp_libs.append((x, xtal.group.number, wps))
                    yield (self.dim, wp_libs, self.elements, self.calculator,
                           self.ref_environments, opt_type, T, niter,
                           early_quit, minimizers)

            # Use the generator to pass args to reduce memory usage
            _xtal, _xs = None, None
            for result in pool.imap_unordered(minimize_from_x_par,
                                              generate_args(),
                                              chunksize=1):
                if result is not None:
                    (_xtals, _xs) = result
                    valid_xtals = self.process_xtals(
                        _xtals, _xs, add_db, symmetrize)
                    xtals_opt.extend(valid_xtals)  # Use deque to reduce memory

            # Remove the duplicate structures
            self.db.update_row_topology(overwrite=False, prefix=self.prefix)
            self.db.clean_structures_spg_topology(dim=self.dim)

            # After each minibatch, delete the local variables and run garbage collection
            del ids, _xtals, _xs
            gc.collect()  # Explicitly call garbage collector to free memory

        xtals_opt = list(xtals_opt)
        print(f"Rank {self.rank} finish optimize_xtals_mproc {len(xtals_opt)}")
        return xtals_opt

    def optimize_reps(self, reps, ncpu=1, opt_type='local',
                      T=0.2, niter=20, early_quit=0.02,
                      add_db=True, symmetrize=False,
                      minimizers=[('Nelder-Mead', 100), ('L-BFGS-B', 100)],
                      N_grids=None):
        """
        Perform optimization for each structure

        Args:
            reps: list of reps
            ncpu (int):
        """
        args = (opt_type, T, niter, early_quit, add_db, symmetrize, minimizers)
        if ncpu == 1:
            valid_xtals = self.optimize_reps_serial(reps, args, N_grids)
        else:
            valid_xtals = self.optimize_reps_mproc(reps, ncpu, args, N_grids)
        return valid_xtals

    def optimize_reps_serial(self, reps, args, N_grids):
        """
        Optimization in multiprocess mode.

        Args:
            reps: list of reps
            ncpu (int): number of parallel python processes
            args: (opt_type, T, n_iter, early_quit, add_db, symmetrize, minimizers)
        """
        xtals_opt = []
        for i, rep in enumerate(reps):
            #print('start', i, rep, len(rep))
            xtal = pyxtal()
            xtal.from_tabular_representation(rep,
                                             normalize=False,
                                             discrete=N_grids)
                                             #verbose=True)
            xtal, sim, _xs = self.optimize_xtal(xtal, i, *args)
            if xtal is not None:
                xtals_opt.append(xtal)
            #else:
            #    print("Debug===="); import sys; sys.exit()
        return xtals_opt

    def optimize_reps_mproc(self, reps, ncpu, args, N_grids):
        """
        Optimization in multiprocess mode.

        Args:
            reps: list of reps
            ncpu (int): number of parallel python processes
            args: (opt_type, T, n_iter, early_quit, add_db, symmetrize, minimizers)
        """

        pool = Pool(processes=ncpu)
        (opt_type, T, niter, early_quit, add_db, symmetrize, minimizers) = args
        xtals_opt = deque()

        # Split the input structures to minibatches
        N_batches = 50 * ncpu
        for _i, i in enumerate(range(0, len(reps), N_batches)):
            start, end = i, min([i+N_batches, len(reps)])
            ids = list(range(start, end))
            print(f"Rank {self.rank} minibatch {start} {end}")
            self.logging.info(f"Rank {self.rank} minibatch {start} {end}")
            self.print_memory_usage()

            def generate_args():
                """
                A generator to yield argument lists for minimize_from_x_par.
                """
                for j in range(ncpu):
                    _ids = ids[j::ncpu]
                    wp_libs = []
                    for id in _ids:
                        rep = reps[id]
                        xtal = pyxtal()
                        discrete = False if N_grids is None else True
                        xtal.from_tabular_representation(rep,
                                                         normalize=False,
                                                         discrete=discrete,
                                                         N_grids=N_grids)
                        x = xtal.get_1d_rep_x()
                        spg, wps, _ = self.get_input_from_ref_xtal(xtal)
                        wp_libs.append((x, spg, wps))
                    yield (self.dim, wp_libs, self.elements, self.calculator,
                           self.ref_environments, opt_type, T, niter,
                           early_quit, minimizers)

            # Use the generator to pass args to reduce memory usage
            _xtal, _xs = None, None
            for result in pool.imap_unordered(minimize_from_x_par,
                                              generate_args(),
                                              chunksize=1):
                if result is not None:
                    (_xtals, _xs) = result
                    valid_xtals = self.process_xtals(
                        _xtals, _xs, add_db, symmetrize)
                    xtals_opt.extend(valid_xtals)  # Use deque to reduce memory

            # Remove the duplicate structures
            self.db.update_row_topology(overwrite=False, prefix=self.prefix)
            self.db.clean_structures_spg_topology(dim=self.dim)

            # After each minibatch, delete the local variables and run garbage collection
            del ids, _xtals, _xs
            gc.collect()  # Explicitly call garbage collector to free memory

        xtals_opt = list(xtals_opt)
        print(f"Rank {self.rank} finish optimize_reps_mproc {len(xtals_opt)}")
        return xtals_opt

    def optimize_xtal(self, xtal, count=0, opt_type='local',
                      T=0.2, niter=20, early_quit=0.02,
                      add_db=True, symmetrize=False,
                      minimizers=[('Nelder-Mead', 100), ('L-BFGS-B', 100)],
                      filename=None):
        """
        Further optimize the input xtal w.r.t reference environment

        Args:
            xtal (instance): pyxtal
        """
        # Change the angle to a better rep
        if xtal.dim == 3 and xtal.lattice is not None and xtal.lattice.ltype in ['triclinic', 'monoclinic']:
            xtal.optimize_lattice(standard=True)
        #xtal.to_file(f'init_{count}.cif')#; print(xtal)
        x = xtal.get_1d_rep_x()
        _, wps, _ = self.get_input_from_ref_xtal(xtal)

        sim0 = self.get_similarity(xtal)
        if xtal.lattice is not None:
            result = minimize_from_x(x, xtal.dim, xtal.group.number, wps,
                                     self.elements, self.calculator,
                                     self.ref_environments,
                                     opt_type=opt_type,
                                     T=T,
                                     niter=niter,
                                     early_quit=early_quit,
                                     minimizers=minimizers,
                                     filename=filename)
            xtal, xs = result
            if result is not None:
                status = xtal.check_validity(self.criteria, verbose=self.verbose)
                sim1 = self.get_similarity(xtal)
                #print("after optim", sim1, status)
            else:
                xtal, xs, status, sim1 = None, None, False, None
        else:
            print("Lattice is None")#, xtal.get_xtal_string())
            xtal, xs, status, sim1 = None, None, False, None
            #import sys; sys.exit()

        if status:
            if symmetrize:
                pmg = xtal.to_pymatgen()
                xtal = pyxtal()
                xtal.from_seed(pmg)
            if add_db:
                self.process_xtal(xtal, [sim0, sim1], count, xs)
            else:
                dicts = {'sim': "{:6.3f} => {:6.3f}".format(sim0, sim1)}
                print(xtal.get_xtal_string(dicts))
        else:
            if self.verbose:
                print('invalid relaxation', count)
                print(xtal.get_xtal_string())
            #import sys; sys.exit()
            #xtal.to_file(f'{count}.cif')
            xtal = None

        return xtal, sim1, xs

    def generate_xtal(self, spg, wps, niter, T=0.2, N_max=5,
                      early_quit=0.03, dump=False, verbose=None,
                      add_db=True, random_state=None):
        """
        Generate a crystal with the desired local environment

        Args:
            spg (int): group number
            wps (list): list of wps for the disired crystal (e.g., [spg, wp1, wp2])
            n_iter (int): number of iterations for basin hopping
            T (float): for basinhopping
            N_max (int): number of maximum
            early_quit (float): threshhold for early termination
            dump (bool): whether or not dump the trajectory
        """
        if verbose is None:
            verbose = self.verbose
        xtal, sim = generate_xtal(self.dim, spg, wps, niter,
                                  self.elements,
                                  self.calculator,
                                  self.ref_environments,
                                  self.criteria,
                                  T, N_max, early_quit, dump,
                                  random_state,
                                  verbose=verbose)

        if xtal is not None and xtal.check_validity(self.criteria):
            if add_db:
                self.process_xtal(xtal, [0, sim], 0)

        return xtal, sim

    def generate_xtals_from_wp_libs(self, wp_libs, N_max=5, ncpu=1,
                                    T=0.2, factor=5, early_quit=0.02):
        """
        Run multiple crystal generation from the given wp_libs
        This is the core part for structure generation.

        Args:
            wp_libs (tuple): (number, spg, wps, dof)
            N_max (int): Number of maximum runs
            ncpu (int): Num of parallel processes
            T (float): basinhopping temperature
            factor (int): the number of Basinhopping iterations = factor * dof
            early_quit (float): early termination for basinhopping
        """

        # Generate xtals
        _args = (self.dim,
                 self.elements,
                 self.calculator,
                 self.ref_environments,
                 self.criteria,
                 T, N_max, early_quit)
        count = 0
        xtals, sims = [], []
        if ncpu == 1:
            for wp_lib in wp_libs:
                (number, spg, wps, dof) = wp_lib
                xtal, sim = generate_xtal(spg, wps, factor*dof, *_args)
                if xtal is not None:
                    xtals.append(xtal)
                    sims.append(sim)

        else:
            N_cycle = int(np.ceil(len(wp_libs)/ncpu))
            print(
                "\n# Parallel Calculation in generate_xtals_from_wp_libs", ncpu, N_cycle)

            args_list = []
            for i in range(ncpu):
                id1 = i * N_cycle
                id2 = min([id1 + N_cycle, len(wp_libs)])
                args_list.append((wp_libs[id1:id2], factor) + _args)

            with ProcessPoolExecutor(max_workers=ncpu) as executor:
                results = [executor.submit(generate_xtal_par, *p)
                           for p in args_list]
                for result in results:
                    (_xtals, _sims) = result.result()
                    xtals.extend(_xtals)
                    sims.extend(_sims)

        for i, xtal in enumerate(xtals):
            self.process_xtal(xtal, [0, sims[i]], i)

        return xtals

    def get_wp_libs_from_spglist(self, spg_list,
                                 num_wp=(None, None),
                                 num_fu=(None, None),
                                 num_dof=(1, 10),
                                 per_spg=30,
                                 ncpu=1):
        """
        Generate wp choices from the list of space groups

        Args:
            spglist (list): list of space group numbers
            num_wp (int): a tuple of (min_wp, max_wp)
            num_fu (int): a tuple of (min_fu, max_fu)
            num_dof (int): a tuple of (min_dof, max_dof)
            per_spg (int): maximum number of wp combinations
            ncpu (int): number of processors
        """

        print('\nGet wp_libs from the given spglist')
        composition = self.composition
        (min_wp, max_wp) = num_wp
        if min_wp is None:
            min_wp = len(composition)
        if max_wp is None:
            max_wp = max([min_wp, len(composition)])
        num_wp = (min_wp, max_wp)

        def process_wp_lib(spg, wp_lib):
            strs = "{:d} wp combos in space group {:d}".format(
                len(wp_lib), spg)
            print(strs)
            self.logging.info(strs)

            if len(wp_lib) > per_spg:
                ids = np.random.choice(
                    range(len(wp_lib)), per_spg, replace=False)
                strs = "Randomly choose {:} wp combinations".format(per_spg)
                wp_lib = [wp_lib[x] for x in ids]
                print(strs)
                self.logging.info(strs)
            return wp_lib

        # Get wp_libs
        wp_libs_total = []
        if ncpu == 1:
            for spg in spg_list:
                wp_lib = generate_wp_lib(
                    [spg], composition, num_wp, num_fu, num_dof)
                wp_lib = process_wp_lib(spg, wp_lib)
                if len(wp_lib) > 0:
                    wp_libs_total.extend(wp_lib)

        else:
            N_cycle = int(np.ceil(len(spg_list)/ncpu))
            args_list = []
            print("\n# Parallel Calculation in get_wp_libs_from_spglist", ncpu, N_cycle)
            for i in range(ncpu):
                id1 = i*N_cycle
                id2 = min([id1+N_cycle, len(spg_list)])
                args_list.append((spg_list[id1:id2],
                                 composition,
                                 num_wp,
                                 num_fu,
                                 num_dof))

            # collect the results
            with ProcessPoolExecutor(max_workers=ncpu) as executor:
                results = [executor.submit(generate_wp_lib_par, *p)
                           for p in args_list]
                for result in results:
                    (spgs, wp_libs) = result.result()
                    for spg, wp_lib in zip(spgs, wp_libs[0]):
                        wp_lib = process_wp_lib(spg, wp_lib)
                        wp_libs_total.extend(wp_lib)

        return sorted(wp_libs_total)

    def get_wp_libs_from_xtals(self, db_file=None,
                               num_wp=(None, None),
                               num_dof=(1, 20),
                               num_atoms=[1, 500]):
        """
        For each struc in the database, get the (spg, wp) info

        Args:
            df_file (str): database file
            num_wp (int): tuple of (min_wp, max_wp)
            num_dof (int): tuple of (min_dof, max_dof)
            num_atoms (int): tuple of (min_at, max_at)
        """

        (min_dof, max_dof) = num_dof
        (min_wp, max_wp) = num_wp
        (min_at, max_at) = num_atoms
        if min_wp is None:
            min_wp = len(self.composition)
        if max_wp is None:
            max_wp = max([min_wp, len(self.composition)])

        wp_libs_total = []
        if db_file is not None:
            strs = "====Loading the structures from {:}".format(db_file)
            print("\n", strs)
            self.logging.info(strs)

            with connect(db_file) as db:
                for row in db.select():
                    atoms = db.get_atoms(row.id)
                    if min_at <= len(atoms) <= max_at:
                        xtal = pyxtal()
                        xtal.from_seed(atoms, tol=0.1)
                        if min_wp <= len(xtal.atom_sites) <= max_wp and \
                                min_dof <= xtal.get_dof() <= max_dof and \
                                xtal.check_validity(self.criteria):
                            t0 = time()
                            sim0 = self.get_similarity(xtal)
                            dicts = {'sim': "0 => {:6.3f}".format(sim0)}
                            strs = xtal.get_xtal_string(dicts)
                            print(strs, "{:.6f}".format(time()-t0))
                            spg, wps, dof = self.get_input_from_ref_xtal(xtal)
                            wp_libs_total.append(
                                (sum(xtal.numIons), spg, wps, dof))

        return sorted(wp_libs_total)

    def import_structures(self, db_file, ids=(None, None),
                          check=True, same_group=True,
                          spglist=range(1, 231), bounds=[1, 500],
                          relax=True,
                          ):
        """
        Import the structures from the external ase database

        Args:
            db_file (str): ase database
            check (boolean): whether or not
            spglist (int): list of spg numbers
            bounds: number of atoms
            energy (boolean): add energy or not
        """
        [lb, ub] = bounds
        count = 0
        strs = "====Loading the structures from {:}".format(db_file)
        print("\n", strs)
        self.logging.info(strs)
        with connect(db_file) as db:
            (min_id, max_id) = ids
            if min_id is None:
                min_id = 1
            if max_id is None:
                max_id = db.count() + 100000
            for row in db.select():
                if min_id <= row.id <= max_id:
                    atoms = row.toatoms()
                    xtal = pyxtal()
                    try:
                        xtal.from_seed(atoms, tol=0.1)
                        if lb <= len(atoms) <= ub and xtal.group.number in spglist:
                            status = True
                        else:
                            status = False
                    except:
                        print("Error in loading structure")
                        status = False

                    # Check structures
                    if status:
                        # Relax the structure
                        if relax:
                            xtal, sim, _ = self.optimize_xtal(
                                xtal, add_db=False)

                        if xtal is not None:
                            energy = row.ff_energy if hasattr(
                                row, 'ff_energy') else None
                            topology = row.topology if hasattr(
                                row, 'topology') else None
                            self.process_xtal(xtal, [0, sim], count,
                                              energy=energy,
                                              topology=topology,
                                              same_group=same_group)
                            count += 1

    def process_xtal(self, xtal, sim, count=0, xs=None, energy=None,
                     topology=None, same_group=True, db=None, check=False):
        """
        Check, print and add xtal to the database

        Args:
            xtal: pyxtal object
            sim: list of two similarity numbers
            count (int): id of the structure in the database
            xs (tuple): list of reps before and after optimization
            energy (float): the energy value
            same_group (bool): keep the same group or not
            db (str): db path
            check (bool): whether or not check if the structure is a duplicate
        """
        if db is None:
            db = self.db
        if check:
            status = db.check_new_structure(xtal, same_group)
        else:
            status = True
        header = "{:4d}".format(count)
        dicts = {'energy': energy,
                 'status': status,
                 'sim': "{:12.3f} => {:6.3f}".format(sim[0], sim[1])
                 }
        strs = xtal.get_xtal_string(dicts, header)
        print(strs)
        self.logging.info(strs)
        kvp = {
            'similarity0': sim[0],
            'similarity': sim[1],
        }
        if xs is not None:
            kvp['x_init'] = np.array2string(xs[0])
            kvp['x_opt'] = np.array2string(xs[1])
        if energy is not None:
            kvp['ff_energy'] = energy
        if topology is not None:
            kvp['topology'] = topology
        if status:
            db.add_xtal(xtal, kvp)


"""
Custom step-function for basin hopping optimization
"""


class RandomDispBounds(object):
    """
    random displacement with bounds:
    see: https://stackoverflow.com/a/21967888/2320035
    """

    def __init__(self, xmin, xmax, id1, id2, stepsize=0.1, dumpfile=None):
        self.xmin = xmin  # bound_min
        self.xmax = xmax  # bound_max
        self.id1 = id1
        self.id2 = id2
        self.stepsize = stepsize
        self.dumpfile = dumpfile

    def __call__(self, x):
        """
        Move the step proportionally within the bound
        """
        # To strongly rebounce the values hitting the wall
        for i in range(self.id2):
            if abs(x[i] - self.xmax[i]) < 0.1:
                # print("To strongly rebounce max values", x[i], self.xmax[i])
                x[i] *= 0.5
            elif abs(x[i] - self.xmin[i]) < 0.1:
                # print("To strongly rebounce min values", x[i], self.xmin[i])
                x[i] *= 2.0

        random_step = np.random.uniform(low=-self.stepsize,
                                        high=self.stepsize,
                                        size=x.shape)
        xnew = x + random_step
        # Cell
        coefs = 1+0.2*(np.random.sample(self.id2)-0.5)
        xnew[:self.id2] *= coefs

        # xyz
        xnew[self.id1:] -= np.floor(xnew[self.id1:])
        # xnew = np.maximum(self.xmin, xnew)
        # xnew = np.minimum(self.xmax, xnew)

        # Randomly introduce compression to prevent non-pd xtal
        if np.random.random() < 0.5:
            id = np.argmax(xnew[:self.id2])
            xnew[id] *= 0.8
        xnew = np.maximum(self.xmin, xnew)
        xnew = np.minimum(self.xmax, xnew)

        # min_step = np.maximum(self.xmin - x, -self.stepsize)
        # max_step = np.minimum(self.xmax - x, self.stepsize)
        # random_step = np.random.uniform(low=min_step, high=max_step, size=x.shape)
        # xnew = x + random_step
        if self.dumpfile is not None:
            with open(self.dumpfile, 'a+') as f0:
                # Initial value
                strs = 'Init: {:9.3f} '.format(10.0)
                for x0 in xnew:
                    strs += '{:8.4f} '.format(x0)
                strs += '\n'
                f0.write(strs)

        return xnew


if __name__ == "__main__":

    xtal = pyxtal()
    xtal.from_spg_wps_rep(194, ['2c', '2b'], [2.46, 6.70])
    cif_file = xtal.to_pymatgen()
    builder = mof_builder(['C'], [1], db_file='reaxff.db',
                          verbose=False)  # True)
    builder.set_descriptor_calculator(mykwargs={'rcut': 1.9})
    builder.set_reference_enviroments(cif_file)
    builder.set_criteria(CN={'C': [3]})
    print(builder)
    print(builder.ref_xtal)
    if False:
        wp_libs = builder.get_wp_libs_from_spglist([191, 179], ncpu=1)
        for wp_lib in wp_libs[:4]:
            print(wp_lib)
        builder.generate_xtals_from_wp_libs(
            wp_libs[4:8], ncpu=2, N_max=4, early_quit=0.05)

    if False:
        spg, wps = 179, ['6a', '6a', '6a', '6a']
        xtals = []
        for x in [
            # [ 9.6244, 2.5459, 0.1749, 0.7701, 0.4501, 0.6114],
            # [15.0223, 1.5013, 0.8951, 0.6298, 0.4530, 0.1876],
            # [10.0129, 2.6424, 0.3331, 0.7246, 0.4719, 0.8628],
            # [10.2520, 3.1457, 0.2367, 0.6994, 0.2522, 0.6533],
            # [ 9.3994,   2.5525,   0.3072,   0.8414,   0.3480,   0.9638],
            [11.1120,   2.6428,   0.2973,   0.7513,   0.4236,   0.8777],
            [7.9522,   2.6057,   0.5922,   0.9268,   0.6081,   0.3077],
        ]:
            xtal = pyxtal()
            xtal.from_spg_wps_rep(spg, wps, x, ['C']*len(wps))
            xtals.append(xtal)
            builder.optimize_xtal(xtal)
        builder.optimize_xtals(xtals)
