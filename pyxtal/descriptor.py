"""
Module for crystal packing descriptor from energy decomposition
"""
import numpy as np
from pyshtools.expand import SHExpandLSQ
import pyshtools as pysh
from scipy.stats import qmc
from scipy.spatial.transform import Rotation
from scipy.optimize import minimize

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
    power = np.sqrt(coef1.spectrum()[s:].sum()*coef2.spectrum()[s:].sum())
    cross = coef1.cross_spectrum(coef2)[s:]
    return cross.sum()/power

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

    res = minimize(fun, angle, args=(coef1, coef2, s), 
            method='Nelder-Mead', options={'maxiter': 20})
    return -res.fun, res.x

def correlation_go(coef1, coef2, M=6, s=0, d_cut=0.9):
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
    phi = np.pi * (3. - np.sqrt(5.))  # golden angle in radians
    for i in range(N):
        y = 1 - (i / float(N - 1)) * 2  # y goes from 1 to -1
        radius = np.sqrt(1 - y * y)  # radius at y
        theta = phi * i  # golden angle increment
        x = np.cos(theta) * radius
        z = np.sin(theta) * radius
        points.append((x, y, z))

    return np.array(points)

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
        theta0 = np.arccos(r_vec[2]/r_mag)
        if abs((r_vec[2] / r_mag) - 1.0) < 10.**(-8.):
            theta0 = 0.0
        elif abs((r_vec[2] / r_mag) + 1.0) < 10.**(-8.):
            theta0 = np.pi
       
        if r_vec[0] < 0.:
            phi0 = np.pi + np.arctan(r_vec[1] / r_vec[0])
        elif 0. < r_vec[0] and r_vec[1] < 0.:
            phi0 = 2 * np.pi + np.arctan(r_vec[1] / r_vec[0])
        elif 0. < r_vec[0] and 0. <= r_vec[1]:
            phi0 = np.arctan(r_vec[1] / r_vec[0])
        elif r_vec[0] == 0. and 0. < r_vec[1]:
            phi0 = 0.5 * np.pi
        elif r_vec[0] == 0. and r_vec[1] < 0.:
            phi0 = 1.5 * np.pi
        else:
            phi0 = 0.
        pts[i, :] = [theta0, phi0]
    if not radian:
        pts = np.degree(pts)

    return pts

def expand_sph(pts, l_max, norm=4, csphase=-1):
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
    thetas, phis, vals = pts[:,0], pts[:,1], pts[:,2]
    thetas = np.degrees(thetas) - 90
    phis = np.degrees(phis)
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
    #get the three most importants ids
    tps = pts[:3]

    xyz0 = np.array([np.sin(tps[0,0])*np.cos(tps[0,1]), np.sin(tps[0,1])*np.sin(tps[0,0]), np.cos(tps[0,0])])
    xyz1 = np.array([np.sin(tps[1,0])*np.cos(tps[1,1]), np.sin(tps[1,1])*np.sin(tps[1,0]), np.cos(tps[1,0])])
    xyz2 = np.array([np.sin(tps[2,0])*np.cos(tps[2,1]), np.sin(tps[2,1])*np.sin(tps[2,0]), np.cos(tps[2,0])])
    line1, line2 = xyz1-xyz0, xyz2-xyz0
    ax1 = np.cross(line1, line2)
    ax1 /= np.linalg.norm(ax1)
    ax0 = np.array([0.0, 0.0, 1.0])
    #
    vec = np.cross(ax0, ax1)
    angle = np.arccos(np.dot(ax0, ax1))
    r = Rotation.from_rotvec(angle * vec)
    angles = r.as_euler('zyz')
    if degrees: angles = np.degrees(angles)
    return angles



class spherical_image():
    """
    A class to handle the crystal packing descriptor from spherical image

    Args:
        xtal: pyxtal structure
        model: 'molecule' or 'contact'
        max_d: maximum intermolecular distances
        max_CN: maximum number of neighboring molecules
        lmax: maximum bandwidth for spherical harmonic expansion
        sigma: Gaussian width to project into the unit sphere
        N: number of grid points on the unit sphere
    """
    def __init__(self, xtal, model='molecule', max_d=10, max_CN=36, 
            lmax=13, sigma=0.1, N=10000):

        self.xtal = xtal
        self.max_d = max_d
        self.max_CN = max_CN
        self.lmax = lmax
        self.sigma = sigma
        self.N = N

        xyzs = fibonacci_sphere(N)
        grids = np.zeros([N, 3])
        grids[:, :2] = xyz2sph(xyzs)

        if model == 'molecule':
            self.pts = self.get_molecules()
        else:
            self.pts = self.get_contacts()

        self.coefs = []
        for pt in self.pts:
            grids[:, 2] = self.calculate_density(pt, xyzs)
            cilm, chi2 = expand_sph(grids, lmax)
            self.coefs.append(pysh.SHCoeffs.from_array(cilm))

    def calculate_density(self, pt, xyzs):
        """
        calculate the projected density on the unit sphere
        """
        vals = np.zeros(len(xyzs))
        for _pt in pt:
            t0, p0, h = _pt
            x0, y0, z0 = np.sin(t0)*np.cos(p0), np.sin(t0)*np.sin(p0), np.cos(t0)
            dst = np.linalg.norm(xyzs - np.array([x0, y0, z0]), axis=1)
            vals += h*np.exp(-(dst**2/(2.0*self.sigma**2))) 
        return vals
            
    def get_molecules(self):
        '''
        compute the spherical images from neighboring molecules
    
        Args:
            scale: scale the coordinates by priciple axis?
    
        Returns
            pts: [N, 3] array, (theta, phi, eng)
        '''
        pts = []
        for i, site in enumerate(self.xtal.mol_sites):
            _, neighs, _, _, engs = self.xtal.get_neighboring_molecules(i, 
                                                    factor=2.2, 
                                                    max_d=self.max_d, 
                                                    CN=self.max_CN)
            xyz, _ = site._get_coords_and_species(absolute=True, first=True) 
            center = site.molecule.get_center(xyz)
            coords = np.zeros([len(neighs), 3])
            for _i, xyz in enumerate(neighs):
                coords[_i, :] = site.molecule.get_center(xyz) - center
            pt = np.zeros([len(coords), 3])
            pt[:, :2] = xyz2sph(coords)
            pt[:, 2] = engs/np.sum(engs)
            pts.append(pt)
        return pts
    
    def get_contacts(self):
        '''
        compute the spherical images from the neighboring distances
    
        Args:
            scale: scale the coordinates by priciple axis?
    
        Returns
            pts: [N, 3] array, (theta, phi, eng)
        '''
        pts = []
        for i, site in enumerate(self.xtal.mol_sites):
            engs, pairs, dists = self.xtal.get_neighboring_dists(i, 
                                            factor=2.2, 
                                            max_d=self.max_d)
            pt = np.zeros([len(pairs), 3])
            pt[:, :2] = xyz2sph(pairs)
            pt[:, 2] = engs/np.sum(engs)
            pts.append(pt)
        return pts
    
    def plot(self, figname=None):
        """
        plot the descriptors
        """
        import matplotlib.pyplot as plt
        import matplotlib.gridspec as gridspec
        
        fig = plt.figure(figsize=(9, 4*len(self.coefs)))
        gs = gridspec.GridSpec(nrows=len(self.coefs), ncols=2, 
                               wspace=0.15, width_ratios=[0.7, 1])
        for i in range(len(self.coefs)):
            ax1 = fig.add_subplot(gs[i, 0])
            ax2 = fig.add_subplot(gs[i, 1])
            coef = self.coefs[i]
            grid = coef.expand()
            grid.plot3d(0, 0, show=False, ax=ax1)
            grid.plot(show=False, ax=ax2, tick_interval=[120, 90], 
                    tick_labelsize=14, axes_labelsize=16)
            ax2.set_xlim([1, 359])
        if figname is None:
            plt.show()
        else:
            plt.savefig(figname)


    def align(self):
        for i in range(len(self.coefs)):
            angles = get_alignment(self.pts[i])
            self.coefs[i] = self.coefs[i].rotate(angles[0], angles[1], angles[2]) 


if __name__ == '__main__':
    
    from pkg_resources import resource_filename
    from pyxtal import pyxtal

    cif_path = resource_filename("pyxtal", "database/cifs/")
    c1 = pyxtal(molecular=True)
    for name in ['benzene', 'aspirin', 'naphthalene']:
        c1.from_seed(seed=cif_path+name+".cif", molecules=[name])
        for model in ['contact', 'molecule']:
            sph = spherical_image(c1, model=model)
            sph.align()
            sph.plot(name+'-'+model+'.png')
