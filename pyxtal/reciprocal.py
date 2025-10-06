"""
Module for Reciprocal Representation Simulation
"""
import importlib.resources
import numpy as np
from monty.serialization import loadfn
from pyxtal.database.element import Element
from pyxtal.XRD import create_index
import torch
from e3nn.o3 import spherical_harmonics

with importlib.resources.as_file(
    importlib.resources.files("pyxtal") / "database" / "atomic_scattering_params.json"
) as path:
    ATOMIC_SCATTERING_PARAMS = loadfn(path)

def bessel_basis(r, nmax=6, r_cut=0.24):
    """
    Spherical Bessel basis functions - excellent for reciprocal space
    """
    from scipy.special import spherical_jn

    r_numpy = r.cpu().numpy()
    r_scaled = r_numpy / r_cut * (nmax + 1) * np.pi

    basis = torch.zeros((r.shape[0], nmax), dtype=r.dtype, device=r.device)
    for n in range(nmax):
        # Add small offset to avoid singularity at r=0
        jn = spherical_jn(n, r_scaled + 1e-10)
        #basis[:, n] = torch.tensor(jn, dtype=r.dtype, device=r.device)
        basis[:, n] = torch.tensor(jn, dtype=r.dtype, device=r.device).view(-1)

    return basis

def gto_basis(r, nmax=6, r_cut=0.24):
    """
    Create a set of radial basis functions

    Args:
        r: radial distances, shape (N, 1)
        nmax: maximum radial quantum number
        r_cut: cutoff radius for normalization (defaults to max radius)

    Returns:
        Tensor of shape (N, nmax)
    """
    # Scale r to [0, 1] range
    r_scaled = r / r_cut

    # Create empty tensor to hold basis functions
    if isinstance(r, torch.Tensor):
        basis = torch.zeros((r.shape[0], nmax), dtype=r.dtype, device=r.device)
    else:
        basis = np.zeros((r.shape[0], nmax))

    # Fill with basis functions (Gaussian-type orbitals)
    for n in range(nmax):
        # GTO-like function: r^n * exp(-alpha * r^2)
        alpha = 1.0 / (n + 1)  # Different width for each basis function
        basis[:, n] = ((r_scaled ** n) * torch.exp(-alpha * r_scaled ** 2)).view(-1)

    return basis  # Shape (N, nmax)

def chebyshev_basis(r, nmax=6, r_cut=0.24):
    """
    Chebyshev polynomial basis - excellent for oscillatory features
    """

    # Scale r to [-1, 1] range for Chebyshev polynomials
    r_scaled = 2 * (r / r_cut) - 1
    basis = torch.zeros((r.shape[0], nmax), dtype=r.dtype, device=r.device)

    # T0(x) = 1, T1(x) = x
    basis[:, 0] = 1
    if nmax > 1:
        basis[:, 1] = r_scaled.view(-1)

    # Recurrence relation: Tn+1(x) = 2x*Tn(x) - Tn-1(x)
    for n in range(2, nmax):
        basis[:, n] = 2 * r_scaled.view(-1) * basis[:, n-1] - basis[:, n-2]

    return basis

class RECP:
    """
    A class to compute the crystal in the reciprocal space.

    Args:
        d_max (float): maximum d-spacing to consider in the reciprocal space
        nmax: int, degree of radial expansion
        lmax: int, degree of spherical harmonic expansion
        alpha: float, gaussian width parameter
    """

    def __init__(self, dmax=6.0, nmax=4, lmax=4, rbasis='chebyshev', res=0.1, sigma=None):
        self.dmax = dmax
        self.nmax = nmax
        self.lmax = lmax
        self.num_bins = int(np.ceil(self.dmax / res))
        self.res = res
        if sigma is None: sigma = 5*self.res
        self.sigma = sigma
        self.rbasis = rbasis
        self.rcut = 0.5 * self.dmax / np.pi

    def __str__(self):
        s = "Reciprocal space expansion with Cutoff: {self.dmax:6.3f} per Ang\n"
        s += "lmax: {self.lmax}, nmax: {self.nmax}, alpha: {self.alpha:.3f}\n"
        return s

    def __repr__(self):
        return str(self)

    def build_reciprocal(self, atoms):
        """
        3x3 representation -> 1x6 (a, b, c, alpha, beta, gamma)
        """
        #print(atoms)
        rec_matrix = atoms.cell.reciprocal()
        hkl_index = create_index()
        hkl_max = np.array([1, 1, 1])
        hkl_min = np.array([-1, -1, -1])

        for index in hkl_index:
            d = np.linalg.norm(np.dot(index, rec_matrix)) * (2 * np.pi)
            multiple = int(np.ceil(self.dmax / d))
            index *= multiple
            for i in range(len(hkl_max)):
                if hkl_max[i] < index[i]:
                    hkl_max[i] = index[i]
                if hkl_min[i] > index[i]:
                    hkl_min[i] = index[i]
        h0, k0, l0 = hkl_min
        h1, k1, l1 = hkl_max
        h = np.arange(h0, h1 + 1)
        k = np.arange(k0, k1 + 1)
        l = np.arange(l0, l1 + 1)

        hkl = np.array(np.meshgrid(h, k, l)).transpose()
        hkl = np.reshape(hkl, [len(h) * len(k) * len(l), 3])
        hkl = hkl[np.where(hkl.any(axis=1))[0]]
        d_hkl = np.linalg.norm(hkl@rec_matrix, axis=1) * (2 * np.pi)

        mask = np.where(d_hkl < self.dmax)[0]
        hkl, d_hkl = hkl[mask], d_hkl[mask]
        #print("hkl_max", h0, h1, k0, k1, l0, l1, d_hkl.shape, "d_hkl", d_hkl.min(), d_hkl.max())

        N_atoms = len(atoms)
        s2 = d_hkl ** 2 / (16 * np.pi ** 2)
        # Compute the atomic scattering factors
        coeffs = np.zeros([N_atoms, 4, 2])
        zs = np.zeros([N_atoms, 1], dtype=int)
        for i, elem in enumerate(atoms.get_chemical_symbols()):
            coeffs[i, :, :] = ATOMIC_SCATTERING_PARAMS[elem]
            zs[i] = Element(elem).z

        tmp1 = np.exp(np.einsum("ij,k->ijk", -coeffs[:, :, 1], s2))  # N*4, M
        tmp2 = np.einsum("ij,ijk->ik", coeffs[:, :, 0], tmp1)  # N*4, N*M
        sfs = np.add(-41.78214 * np.einsum("ij,j->ij", tmp2, s2), zs)  # N*M, M -> N*M
        # to add dampling factor to ensure the decay to 0?

        # Compute the structure factors
        const = -2j * np.pi
        positions = atoms.get_scaled_positions()#; print("positions", positions)
        g_dot_rs = np.dot(positions, hkl.T)  # N_atoms * M_hkl
        exps = np.exp(const * g_dot_rs)
        fs = np.sum(sfs * exps, axis=0)
        intensities = (fs * fs.conjugate()).real  # M
        masks = np.where(intensities > 1e-4)[0]
        hkl, intensities, d_hkl, fs = hkl[masks], intensities[masks], d_hkl[masks], fs[masks]
        sfs = sfs[:, masks]
        I0 = np.sum(zs)**2
        intensities /= I0  # Normalize the intensities
        intensities *= np.cos(0.5 * np.pi * d_hkl / self.dmax)  # Apply the Gaussian factor
        #print("intensities", intensities.shape, intensities.min(), intensities.max())
        #max_idx = np.argmax(intensities); print("max", intensities[max_idx], hkl[max_idx], d_hkl[max_idx], fs[max_idx])
        #min_idx = np.argmin(intensities); print("min", intensities[min_idx], hkl[min_idx], d_hkl[min_idx], fs[min_idx])
        # find the id of hkl=[1,1,1]
        #id_111 = np.where(np.all(hkl == [1, 1, 1], axis=1))[0]
        #if len(id_111) > 0:
        #    print("hkl=[1,1,1]", intensities[id_111], d_hkl[id_111], sfs[0, id_111], fs[id_111])
        return hkl@rec_matrix, intensities, d_hkl

    def compute(self, atoms, norm=False):
        """
        d for any give abitray [h,k,l] index
        """
        coords, vals, ds = self.build_reciprocal(atoms)
        p = self.compute_sph_torch(coords, vals, norm=norm)
        #print("\np shape:", p.shape, "p min:", p.min(), "p max:", p.max())
        rdf = self.compute_rdf(ds, vals)
        #D = np.concatenate([p, rdf], axis=0)
        #print("D shape:", D.shape, "D min:", D.min(), "D max:", D.max())

        return p, rdf

    def compute_rdf(self, ds, vals):
        """
        Get the radial distribution function (RDF) from the d-spacing and values.
        """
        from scipy.ndimage import gaussian_filter1d

        bins = np.linspace(0, self.dmax, self.num_bins)
        rdf, _ = np.histogram(ds, bins=bins, weights=vals)
        rdf = gaussian_filter1d(rdf, sigma=self.sigma)
        #ids = np.where(rdf > 0.01)[0]; print(len(ds)); print("loc", bins[ids][:5]); print("pek", rdf[ids][:5])
        #print("rdf shape:", rdf.shape, "rdf min:", rdf.min(), "rdf max:", rdf.max())
        return rdf

    def compute_sph_torch(self, xyz, v, norm=False):
        """
        Compute a descriptor using spherical harmonics and radial basis functions.

        Args:
            xyz: Tensor of shape (N, 3) representing 3D coordinates.
            v: Tensor of shape (N,) representing scalar values at each point.
            norm: Whether to normalize the final descriptor.

        Returns:
            Tensor of shape (sum(2l+1) * nmax,) representing the descriptor.
        """
        # Convert NumPy arrays to PyTorch tensors
        xyz = torch.tensor(xyz, dtype=torch.float32)
        v = torch.tensor(v, dtype=torch.float32).view(-1, 1)  # Add a dimension for repeat
        r_hat = xyz / (torch.norm(xyz, dim=1, keepdim=True) + 1e-12)  # unit direction
        r = torch.norm(xyz, dim=1, keepdim=True)  # radial distances
        #print("Debug xyz", xyz.shape, "r_hat", r_hat, "r", r)

        # Compute spherical harmonics up to lmax
        if self.rbasis == 'chebyshev':
            R = chebyshev_basis(r, self.nmax, self.rcut)
        elif self.rbasis == 'gto':
            R = gto_basis(r, self.nmax, self.rcut)
        else:
            R = bessel_basis(r, self.nmax, self.rcut)
        #print("Debug Radial", R)

        degrees = list(range(0, self.lmax + 1))
        #degrees = list(range(0, self.lmax + 1, 2))
        Y = spherical_harmonics(degrees, r_hat, normalize=False, normalization='norm')#component')
        #Y = spherical_harmonics(degrees, r_hat, normalize=False, normalization='component')
        #print("Debug Spherical", Y.shape, "Y min:", Y.min(), "Y max:", Y.max())

        descriptor = []
        # For each radial basis function
        for n in range(self.nmax):
            r_basis = R[:, n:n+1]  # Shape (N, 1)

            # Weight the spherical harmonics by this radial basis
            f_n = v * r_basis * Y  # Shape (N, sum(2l+1))

            # Process by angular momentum
            offset = 0
            for l in degrees:
                #print(f"\nDebug l={l}, n={n}, offset={offset}")
                dim = 2 * l + 1
                Y_l = f_n[:, offset:offset + dim]
                c_nl = torch.sum(Y_l, dim=0)  # Weighted sum over all points
                power_l = torch.mean(c_nl ** 2)
                #print(f"\nDebug Y_l at basis {n}, l={l}")
                #print("point 0", Y_l[0, :], (Y_l[0, :]**2).sum())
                #print("point 1", Y_l[1, :], (Y_l[1, :]**2).sum())
                #print("point 2", Y_l[2, :], (Y_l[2, :]**2).sum())
                #print("Total power:", power_l.item())

                descriptor.append(power_l)
                offset += dim
        descriptor = torch.stack(descriptor)
        if norm:
            norm = torch.linalg.norm(descriptor)
            descriptor /= (norm + 1e-9)  # Add epsilon for stability
        return descriptor

    def plot(self, data, filename='reciprocal.png'):
        """
        Plot the computed reciprocal space representation and RDF.

        Args:
            data: Tuple of (rdf, p, label) where rdf is the radial distribution function
                  and p is the expansion coefficients, label is a string for the plot title.
            filename: Name of the file to save the plot.
        """
        import matplotlib.pyplot as plt

        plt.figure(figsize=(12, 4))
        plt.subplot(1, 2, 1)
        for (p, rdf, label) in data:
            #plt.plot(p, label=label, alpha=0.5, lw=0.9, marker='o', markersize=3, linestyle='-')
            plt.plot(p, label=label, alpha=0.5, lw=0.9)

        plt.xlabel(f'Index ({len(p)})')
        plt.ylabel('Expansion Coefficients')
        plt.legend(loc=1)

        plt.subplot(1, 2, 2)
        x = np.linspace(0, self.rcut, len(rdf))
        for (p, rdf, label) in data:
            plt.plot(x, rdf, label=label, alpha=0.5, lw=0.9)
        plt.xlabel(f'd-spacing (per Angstrom)')
        plt.ylabel('RDF')
        plt.legend(loc=1)

        plt.tight_layout()
        plt.savefig(filename, dpi=300)
        plt.close()

    def reconstruction(self, spg, wps, elements, rep0, P_ref, rdf_ref, verbose=False):
        """
        Generate a crystal with the desired local P_ref

        Args:
            spg (int): pyxtal.symmetry.Group object
            wps: list of wps for the disired crystal (e.g., [wp1, wp2])
            P_ref: reference enviroment

        Returns:
            xtal and its mse loss
        """
        torch.autograd.set_detect_anomaly(True)

        def apply_bounds(tensor):
            """Clamps tensor values between 0 and 1."""
            with torch.no_grad():
                tensor.clamp_(0.0, 1.0)

        # Clone and enable gradients for `reps`
        rep = rep.clone().detach().requires_grad_(True)
        generators = generators.clone().detach()

        # Choose optimizer
        optimizer = torch.optim.Adam([rep_batch], lr=lr)
        scheduler = StepLR(optimizer, step_size=50, gamma=0.1)

        # Optimization loop
        for step in range(num_steps):
            optimizer.zero_grad()

            # Compute losses per sample (B,)
            loss = self.loss(spg, wps, elements, P_ref, RDF_ref)
            loss.backward(torch.ones_like(loss))
            torch.nn.utils.clip_grad_norm_(rep_batch, max_norm=10.0)  # Gradient clipping

            # Step the scheduler
            optimizer.step()
            if step > 100: scheduler.step()

            if verbose and step % 1 == 0:
                print(f"Step {step}, {loss_sum:.6f}, LR={scheduler.get_last_lr()[0]:.6f}")
            if step + 1 == num_steps:
                print(f"stopping at last iteration")
        xtal =
        return rep.detach(), losses.detach()

    def loss(self, spg, wps, elements, P_ref, RDF_ref):
        res =  WP.get()
        p, xrd, rdf = self.compute()
        loss1 = torch.sum()
        return loss





if __name__ == "__main__":
    from pyxtal import pyxtal

    xtal1 = pyxtal(); xtal1.from_prototype('diamond') #; xtal.to_file('dia.cif'); print(xtal)
    xtal_sub = xtal1.subgroup_once(H=141, eps=0.1); print(xtal_sub)
    xtal2 = pyxtal(); xtal2.from_prototype('graphite')
    xtal3 = xtal1.copy(); xtal3.substitute({'C': 'Si'});
    xtal3.lattice = xtal3.lattice.scale(1.52); print(xtal3)

    recp = RECP(dmax=8.0, nmax=5, lmax=3)
    p1, rdf1 = recp.compute(xtal1.to_ase())
    p2, rdf2 = recp.compute(xtal1.to_ase()*2)
    p3, rdf3 = recp.compute(xtal_sub.to_ase())
    p4, rdf4 = recp.compute(xtal2.to_ase())
    p5, rdf5 = recp.compute(xtal3.to_ase())
