import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from ase.db import connect
from pyxtal import pyxtal
from pyxtal.XRD_indexer import get_cell_from_thetas, CellManager, WPManager, XtalManager
from pyxtal.XRD import Similarity
from pyxtal.interface.ase_opt import ASE_relax
from math import gcd

database = connect('demo.db')

for row in database.select():
    xtal0 = pyxtal()
    xtal0.from_seed(row.toatoms())
    spg = xtal0.group.number
    if spg < 75: continue
    print(xtal0)
    found = False

    # Draw PXRD
    my_pxrd = xtal0.get_XRD(thetas=[10, 80], SCALED_INTENSITY_TOL=0.01)
    x1, y1 = my_pxrd.get_plot(bg_ratio=0)
    y1 = np.array(y1)
    peaks, _ = find_peaks(y1, height=1, prominence=2.5)

    composition = {}
    divsor = gcd(*xtal0.numIons)
    for i, el in enumerate(xtal0.species):
        composition[el] = xtal0.numIons[i] // divsor

    pxrd = np.zeros([len(x1), 2])
    pxrd[:, 0], pxrd[:, 1] = x1, y1
    peak_dicts = {'locations': x1[peaks], 'intensities': y1[peaks]}
    print(peak_dicts)
    #print(my_pxrd.by_hkl(N_max=50))

    solutions = get_cell_from_thetas(spg, peak_dicts['locations'], max_mismatch=20,
                                     hkl_max=(2, 3, 10), max_square=20,
                                     theta_tol=0.5, verbose=False)

    sols = [(spg, sol['cell'], len(sol['mis_matched_peaks'])) for sol in solutions]
    cells = CellManager.consolidate(sols, merge_tol=0.15)

    # Store all results
    eng_min = 1e10
    for cell in cells[:2]:
        print(f"Trying cell: {cell.dims}, missing peaks: {cell.missing}")
        wp_manager = WPManager(spg, cell.dims, composition)
        sols = wp_manager.get_wyckoff_positions()
        for sol in sols:
            (spg, comp, lattice, wp_ids, num_wps, dof) = sol
            #print(f"{comp}, {lattice}, {spg}, WPs: {wp_ids}")
            xm = XtalManager(spg, composition.keys(), comp, lattice, wp_ids)
            for i in range(2*xm.dof + 1):
                xtal = xm.generate_structure()
                if not xtal.valid: continue
                if xm.dof > 0:
                    atoms = ASE_relax(xtal.to_ase(), opt_lat=False, fmax=0.05, step=50, logfile='ase.log')
                    if atoms is None: continue
                    xtal.from_seed(atoms)
                else:
                    atoms = ASE_relax(xtal.to_ase(), opt_lat=False, fmax=0.05, step=10, logfile='ase.log')
                    if atoms is None: continue
                xrd = xtal.get_XRD(thetas=[10, 80], SCALED_INTENSITY_TOL=0.01)
                x2, y2 = xrd.get_plot(bg_ratio=0)
                sim = Similarity((x1, y1), (x2, y2)).value
                eng = atoms.get_potential_energy()/len(atoms)
                stress = atoms.get_stress()[:3].mean()
                print(xtal.get_xtal_string(), f"{sim:.3f}, {eng:.3f}, {stress:.3f}")
                if eng < eng_min and abs(stress) < 0.5:
                    eng_min = eng
                    print('+++++++++++++++++++++Update energy', eng_min)

                if sim > 0.9:
                    plt.figure(figsize=(8, 3))
                    strs = f'{xtal.get_xtal_string()}, Energy: {eng:.3f} eV/atom'
                    plt.title(strs)
                    plt.plot(x1, y1, color='green', alpha=0.5, label='Observed')
                    plt.plot(x2, y2, color='blue', alpha=0.5, label=f'Match ({sim:.3f})')
                    plt.plot(x1[peaks], y1[peaks], "x", c='orange')
                    plt.legend()
                    plt.xlabel('2θ (degrees)')
                    plt.ylabel('Intensity (a.u.)')
                    plt.savefig(f'Match_{row.id}.png', dpi=300)
                    plt.close()
                    xtal.to_file(f'Match_{row.id}.cif')
                    xtal0.to_file(f'Original_{row.id}.cif')
                    if eng - eng_min < 1e-2:
                        found = True
                        break
                        #import sys; sys.exit()
            if found:
                break
        if found:
            break
