"""
Module for crystal packing descriptor from energy decomposition
"""

import numpy as np
from scipy.optimize import minimize
from scipy.spatial.transform import Rotation
from scipy.special import sph_harm
from scipy.stats import qmc


def _qlm(dists, l=4):
    """
    Calculates the vector associated with an atomic site and
    one of its neighbors

    Args:
        distss: a list of distance vectors
        l:  free integer quantum number
    Returns:
        q: numpy array(complex128), the complex vector qlm normalized
            by the number of nearest neighbors
    """
    # initiate variable as a complex number
    q = np.zeros(2 * l + 1, dtype=np.complex128)
    neighbors_count = len(dists)

    for i, m in enumerate(range(-l, l + 1)):
        for _j, r_vec in enumerate(dists):
            # find the position vector of the site/neighbor pair
            r_mag = np.linalg.norm(r_vec)
            theta = np.arccos(r_vec[2] / r_mag)
            if abs((r_vec[2] / r_mag) - 1.0) < 10.0 ** (-8.0):
                theta = 0.0
            elif abs((r_vec[2] / r_mag) + 1.0) < 10.0 ** (-8.0):
                theta = np.pi

            # phi
            if r_vec[0] < 0.0:
                phi = np.pi + np.arctan(r_vec[1] / r_vec[0])
            elif r_vec[0] > 0.0 and r_vec[1] < 0.0:
                phi = 2 * np.pi + np.arctan(r_vec[1] / r_vec[0])
            elif r_vec[0] > 0.0 and r_vec[1] >= 0.0:
                phi = np.arctan(r_vec[1] / r_vec[0])
            elif r_vec[0] == 0.0 and r_vec[1] > 0.0:
                phi = 0.5 * np.pi
            elif r_vec[0] == 0.0 and r_vec[1] < 0.0:
                phi = 1.5 * np.pi
            else:
                phi = 0.0

            q[i] += sph_harm(m, l, phi, theta)
    # normalize by number of neighbors
    return q / neighbors_count


def correlation(coef1, coef2, angle=None, s=0):
    """
    Compute the correlation between to sph coefs

    Args:
        coef1: sph coefficients 1
        coef2: sph coefficients 2
        angle: [alpha, beta, gamma]
        s: starting index of coefs

    Return:
        distance scaled in [0, 1]
    """
    if angle is not None:
        coef2 = coef2.rotate(angle[0], angle[1], angle[2], degrees=True)
    power = np.sqrt(coef1.spectrum()[s:].sum() * coef2.spectrum()[s:].sum())
    cross = coef1.cross_spectrum(coef2)[s:]
    return cross.sum() / power


def correlation_opt(coef1, coef2, angle, s=0):
    """
    Compute the correlation between two sph coefs

    Args:
        coef1: sph coefficients 1
        coef2: sph coefficients 2
        angle: [alpha, beta, gamma]
        s: starting index of coefs

    Return:
        distance scaled in [0, 1]
    """

    def fun(x0, coef1, coef2, s):
        coef = coef2.rotate(x0[0], x0[1], x0[2], degrees=True)
        return -correlation(coef1, coef, s=s)

    res = minimize(
        fun,
        angle,
        args=(coef1, coef2, s),
        method="Nelder-Mead",
        options={"maxiter": 20},
    )
    return -res.fun, res.x


def correlation_go(coef1, coef2, M=6, s=0, d_cut=0.92):
    """
    global optimization of two coefs based on quasi random sampling

    Args:
        coef1: sph coefficients 1
        coef2: sph coefficients 2
        M: 2^M sampling points
        s: starting index of coefs

    Return:
        distance scaled in [0, 1]

    """
    sampler = qmc.Sobol(d=3, scramble=False)
    sample = sampler.random_base2(m=M)
    sample = qmc.scale(sample, [-180, -90, -180], [180, 90, 180])

    ds, angles = [], []
    for angle in sample:
        d, angle = correlation_opt(coef1, coef2, angle, s=0)
        ds.append(d)
        angles.append(angle)
        if d > d_cut * 1.1:
            break
    ds = np.array(ds)
    id = np.argmax(ds)
    return ds[id], angles[id]


def fibonacci_sphere(N=1000):
    """
    Sampling the sphere grids

    Args:
        N: number of pts to generate

    Returns:
        3D points array in Cartesian coordinates
    """
    points = []
    phi = np.pi * (3.0 - np.sqrt(5.0))  # golden angle in radians
    for i in range(N):
        y = 1 - (i / float(N - 1)) * 2  # y goes from 1 to -1
        radius = np.sqrt(1 - y * y)  # radius at y
        theta = phi * i  # golden angle increment
        x = np.cos(theta) * radius
        z = np.sin(theta) * radius
        points.append((x, y, z))

    return np.array(points)


def cart2sph(x, y, z):
    """
    convert the x, y, z to spherical coordinates (phi, theta, r)
    phi: [-pi, pi]
    theta: [-pi/2, pi/2]
    """
    hxy = np.hypot(x, y)
    r = np.hypot(hxy, z)
    theta = np.arctan2(z, hxy)
    phi = np.arctan2(y, x)
    return phi, theta, r


def sph2cart(phi, theta, r):
    """
    convert spherical coordinates (phi, theta, r) to Cartesian (x, y, z)
    """
    rcos_theta = r * np.cos(theta)
    x = rcos_theta * np.cos(phi)
    y = rcos_theta * np.sin(phi)
    z = r * np.sin(theta)
    return x, y, z


def xyz2sph(xyzs, radian=True):
    """
    convert the vectors (x, y, z) to the sphere representation (theta, phi)

    Args:
        xyzs: 3D xyz coordinates
        radian: return in radian (otherwise degree)
    """
    pts = np.zeros([len(xyzs), 2])
    for i, r_vec in enumerate(xyzs):
        r_mag = np.linalg.norm(r_vec)
        theta0 = np.arccos(r_vec[2] / r_mag)
        if abs((r_vec[2] / r_mag) - 1.0) < 10.0 ** (-8.0):
            theta0 = 0.0
        elif abs((r_vec[2] / r_mag) + 1.0) < 10.0 ** (-8.0):
            theta0 = np.pi

        if r_vec[0] < 0.0:
            phi0 = np.pi + np.arctan(r_vec[1] / r_vec[0])
        elif r_vec[0] > 0.0 and r_vec[1] < 0.0:
            phi0 = 2 * np.pi + np.arctan(r_vec[1] / r_vec[0])
        elif r_vec[0] > 0.0 and r_vec[1] >= 0.0:
            phi0 = np.arctan(r_vec[1] / r_vec[0])
        elif r_vec[0] == 0.0 and r_vec[1] > 0.0:
            phi0 = 0.5 * np.pi
        elif r_vec[0] == 0.0 and r_vec[1] < 0.0:
            phi0 = 1.5 * np.pi
        else:
            phi0 = 0.0
        pts[i, :] = [theta0, phi0]
    if not radian:
        pts = np.degree(pts)

    return pts


def expand_sph(pts, l_max, norm=4, csphase=-1):
    from pyshtools.expand import SHExpandLSQ

    """
    Transform the grid points to spherical harmonics

    Args:
        pts: 3D array
            (thetas, phis, vals) with a length of N

        lmax: Integer
            The maximum degree of spherical harmonic.

        coeff_norm: integer
            The normalization of SHExpandLSQ().
            1 (default) = Geodesy 4-pi normalized harmonics;
            2 = Schmidt semi-normalized harmonics;
            3 = unnormalized harmonics;
            4 = orthonormal harmonics.

        csphase: Integer
            whether (-1) or not (1) apply the Condon-Shortley phase factor.

    Return:
        cilm: float, dimension (2, lmax+1, lmax+1)
            Coefficients of the spherican harmonic

        chi2: float
            The residual sum of squares misfit for an overdetermined inversion.
    """
    thetas, phis, vals = pts[:, 0], pts[:, 1], pts[:, 2]
    phis = np.degrees(phis)
    # shift from (0, 180) to (-90, 90)
    thetas = np.degrees(thetas)
    thetas -= 90

    # if thetas is within [0, 180]
    # print('check thetas', thetas.min(), thetas.max())
    # if abs(thetas.min()) < 1e-3: thetas -= 90
    cilm, chi2 = SHExpandLSQ(vals, thetas, phis, l_max, norm=norm, csphase=csphase)

    return cilm, chi2


def get_alignment(pts, degrees=True):
    """
    Here we define the equator is the plane with three most important neighbors.
    Get the required rotation angles to get that representation.

    Args:
        pts; important points (3D array)

    Returns:
        angles: [alpha, beta, gamma]
    """
    # get the three most importants ids
    tps = pts[:3]

    xyz0 = np.array(
        [
            np.sin(tps[0, 0]) * np.cos(tps[0, 1]),
            np.sin(tps[0, 1]) * np.sin(tps[0, 0]),
            np.cos(tps[0, 0]),
        ]
    )
    xyz1 = np.array(
        [
            np.sin(tps[1, 0]) * np.cos(tps[1, 1]),
            np.sin(tps[1, 1]) * np.sin(tps[1, 0]),
            np.cos(tps[1, 0]),
        ]
    )
    xyz2 = np.array(
        [
            np.sin(tps[2, 0]) * np.cos(tps[2, 1]),
            np.sin(tps[2, 1]) * np.sin(tps[2, 0]),
            np.cos(tps[2, 0]),
        ]
    )
    line1, line2 = xyz1 - xyz0, xyz2 - xyz0
    ax1 = np.cross(line1, line2)
    ax1 /= np.linalg.norm(ax1)
    ax0 = np.array([0.0, 0.0, 1.0])
    #
    vec = np.cross(ax0, ax1)
    angle = np.arccos(np.dot(ax0, ax1))
    r = Rotation.from_rotvec(angle * vec)
    angles = r.as_euler("zyz")
    if degrees:
        angles = np.degrees(angles)
    return angles


class spherical_image:
    import pyshtools as pysh

    """
    A class to handle the crystal packing descriptor from spherical image

    Args:
        xtal: pyxtal structure
        model: 'molecule' or 'contact'
        max_d: maximum intermolecular distances
        lmax: maximum bandwidth for spherical harmonic expansion
        sigma: Gaussian width to project into the unit sphere
        N: number of grid points on the unit sphere
    """

    def __init__(self, xtal, model="molecule", max_d=10, factor=2.2, lmax=13, sigma=0.1, N=10000):
        for i in range(len(xtal.mol_sites)):
            try:
                numbers = xtal.mol_sites[i].molecule.mol.atomic_numbers
                if numbers.count(7) > 0 or numbers.count(8) > 0:
                    xtal.mol_sites[i].molecule.set_labels()
            except:
                print("Warning! Needs the smiles information!")
        self.xtal = xtal
        self.max_d = max_d
        self.factor = factor
        self.lmax = lmax
        self.sigma = sigma
        self.N = N

        xyzs = fibonacci_sphere(N)
        grids = np.zeros([N, 3])
        grids[:, :2] = xyz2sph(xyzs)

        if model == "molecule":
            self.pts = self.get_molecules()
        else:
            self.pts = self.get_contacts()
        self.coefs = []
        for pt in self.pts:
            grids[:, 2] = self.calculate_density(pt, xyzs)
            cilm, chi2 = expand_sph(grids, lmax)
            self.coefs.append(pysh.SHCoeffs.from_array(cilm))
        self.ds = np.ones(len(self.coefs))

    def calculate_density(self, pt, xyzs):
        """
        calculate the projected density on the unit sphere
        """
        vals = np.zeros(len(xyzs))
        for _pt in pt:
            t0, p0, h = _pt
            x0, y0, z0 = np.sin(t0) * np.cos(p0), np.sin(t0) * np.sin(p0), np.cos(t0)
            dst = np.linalg.norm(xyzs - np.array([x0, y0, z0]), axis=1)
            vals += h * np.exp(-(dst**2 / (2.0 * self.sigma**2)))
        return vals

    def get_molecules(self):
        """
        compute the spherical images from neighboring molecules

        Returns:
            pts: [N, 3] array, (theta, phi, eng)
        """
        pts = []
        for i, site in enumerate(self.xtal.mol_sites):
            _, neighs, comps, _, engs = self.xtal.get_neighboring_molecules(
                i, factor=self.factor, max_d=self.max_d, ignore_E=False
            )
            xyz, _ = site._get_coords_and_species(absolute=True, first=True)
            center = site.molecule.get_center(xyz)
            coords = np.zeros([len(neighs), 3])
            for _i, xyz in enumerate(neighs):
                # coords[_i, :] = site.molecule.get_center(xyz) - center
                coords[_i, :] = self.xtal.molecules[comps[_i]].get_center(xyz) - center
            pt = np.zeros([len(coords), 3])
            pt[:, :2] = xyz2sph(coords)
            pt[:, 2] = engs / np.sum(engs)
            pts.append(pt)
        return pts

    def get_contacts(self):
        """
        Compute the spherical images from the neighboring distances

        Returns:
            pts: [N, 3] array, (theta, phi, eng)
        """
        pts = []
        for i, _site in enumerate(self.xtal.mol_sites):
            engs, pairs, dists = self.xtal.get_neighboring_dists(i, factor=self.factor, max_d=self.max_d)
            pt = np.zeros([len(pairs), 3])
            pt[:, :2] = xyz2sph(pairs)
            pt[:, 2] = engs / np.sum(engs)
            pts.append(pt)
        return pts

    def plot_sph_images(self, lmax=None, figname=None, molecule=False):
        """
        Plot the spherical images in both 3d and 2d

        Args:
            lmax: maximum truncation
            figname: name of figure file
            molecule: draw 2D molecule diagram or not
        """
        import matplotlib.gridspec as gridspec
        import matplotlib.pyplot as plt

        if molecule:
            nrows = len(self.coefs) + 1
            shift = 1
        else:
            nrows = len(self.coefs)
            shift = 0

        fig = plt.figure(figsize=(9, 4 * nrows))
        gs = gridspec.GridSpec(nrows=nrows, ncols=2, wspace=0.15, width_ratios=[0.7, 1])
        if molecule:
            from rdkit import Chem
            from rdkit.Chem import Draw

            smi = ""
            for i, m in enumerate(self.xtal.molecules):
                smi += m.smile
                if i + 1 < len(self.xtal.molecules):
                    smi += "."

            m = Chem.MolFromSmiles(smi)
            im = Draw.MolToImage(m)
            ax0 = fig.add_subplot(gs[0, :])
            plt.imshow(im)
            ax0.axis("off")

        if lmax is None:
            lmax = self.lmax
        elif lmax > self.lmax:
            print("Warning: Cannot set lmax greater than the ", self.lmax)
            lmax = self.lmax

        for i in range(len(self.coefs)):
            ax1 = fig.add_subplot(gs[i + shift, 0])
            ax2 = fig.add_subplot(gs[i + shift, 1])
            coef = self.coefs[i]
            grid = coef.expand(lmax=lmax)
            grid.plot3d(0, 0, title=f"{self.ds[i]:6.3f}", show=False, ax=ax1)
            grid.plot(
                show=False,
                ax=ax2,
                tick_interval=[120, 90],
                tick_labelsize=14,
                axes_labelsize=16,
            )
            ax2.set_xlim([1, 359])

        if figname is None:
            plt.show()
        else:
            plt.savefig(figname)

    def plot_real_image(self, id=0):
        """
        Plot the real molecular contacts in the crystal
        """
        return self.xtal.show_mol_cluster(id, factor=self.factor, max_d=self.max_d, ignore_E=False, plot=False)

    def align(self, M=6):
        """
        Align spherical image in a way that three most important contributions
        are parallel to the equatorial plane. Experimental stage now!

        Args:
            M: number of power in quasi random sampling
        """
        self.coefs[0]
        angles = get_alignment(self.pts[0])
        self.coefs[0] = self.coefs[0].rotate(angles[0], angles[1], angles[2])

        for i in range(1, len(self.coefs)):
            coef1 = self.coefs[i]
            d, angles = correlation_go(self.coefs[0], coef1, M=M)
            self.coefs[i] = coef1.rotate(angles[0], angles[1], angles[2])
            self.ds[i] = d

    def rotate(self, alpha=0, beta=0, gamma=0):
        """
        uniformly rotate the coefs

        Args:
            alpha: rotation in degrees
            beta: rotation in degress
            gamma: rotation in degress
        """
        for i in range(len(self.coefs)):
            self.coefs[i] = self.coefs[i].rotate(alpha, beta, gamma)

    def get_similarity(self, sph2, M=6, cutoff=0.95):
        """
        Compute the similarity matrix between two sphs

        Args:
            sph2: the 2nd sph class
            M: number of power in quasi random sampling
            cutoff: cutoff similarity to terminate search early
        """
        S = np.zeros([len(self.coefs), len(sph2.coefs)])
        for i in range(len(self.coefs)):
            coef1 = self.coefs[i]
            for j in range(len(sph2.coefs)):
                coef2 = sph2.coefs[j]
                d, _ = correlation_go(coef1, coef2, M=M, d_cut=cutoff)
                S[i, j] = d
        return S


class orientation_order:
    """
    Computes the Steinhardt orientation order parameters

    Args:
        xtal: pyxtal structure
        max_d: maximum intermolecular distances
        lmax: maximum bandwidth for spherical harmonic expansion
    """

    def __init__(self, xtal, max_CN=14):
        self.xtal = xtal
        self.max_CN = max_CN
        self.dists = self.get_neighbors()

    def get_neighbors(self):
        """
        get neighboring molecules

        Returns:
            pts: [N, 3] array, (theta, phi, eng)
        """
        pts = []
        for i, site in enumerate(self.xtal.mol_sites):
            _, neighs, comps, _, engs = self.xtal.get_neighboring_molecules(i)
            xyz, _ = site._get_coords_and_species(absolute=True, first=True)
            center = site.molecule.get_center(xyz)
            coords = np.zeros([len(neighs), 3])
            # print(len(neighs))
            if len(neighs) > self.max_CN:
                neighs = neighs[: self.max_CN]
            for _i, xyz in enumerate(neighs):
                coords[_i, :] = self.xtal.molecules[comps[_i]].get_center(xyz) - center
            pts.append(coords)
        return pts

    def get_parameters(self, ls=None):
        """
        Computes
        Args:
            center: center xyz coordinate
            neighbors: a list of neighboring xyz coordinates
            weights: a list of weights for each neighbor
        Returns:
             q: numpy array(complex128), the complex vector qlm normalized
                by the number of nearest neighbors
        """
        if ls is None:
            ls = [4, 6]
        qs = []
        for dist in self.dists:
            for l in ls:
                (4 * np.pi) / (2 * l + 1)
                qlms = _qlm(dist, l)
                dot = float(np.sum(qlms * np.conjugate(qlms)))
                qs.append(np.sqrt((4 * np.pi) / (2 * l + 1) * dot))

        return qs


if __name__ == "__main__":
    import importlib.resources

    from pyxtal import pyxtal

    with importlib.resources.as_file(importlib.resources.files("pyxtal") / "database" / "cifs") as path:
        cif_path = path
    c1 = pyxtal(molecular=True)
    for name in ["benzene", "resorcinol", "aspirin", "naphthalene"]:
        c1.from_seed(seed=str(cif_path / f"{name}.cif"), molecules=[name])
        for model in ["contact", "molecule"]:
            print(name, model)
            sph = spherical_image(c1, model=model, lmax=18)
            sph.align()
            sph.plot_sph_images(figname=name + "-" + model + ".png")

    # for name in ['BENZEN', 'ACSALA', 'RESORA']:
    #    c1.from_CSD(name)
    #    for model in ['contact', 'molecule']:
    #        sph = spherical_image(c1, model=model, lmax=18)
    #        sph.align()
    #        sph.plot_sph_images(figname=name+'-'+model+'.png', molecule=True)
