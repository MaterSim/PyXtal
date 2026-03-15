"""
Module for XRD simulation
"""

import collections
import importlib.resources
import os

import numpy as np
from monty.serialization import loadfn
from scipy.interpolate import interp1d
from scipy.special import erf, wofz
from pyxtal.database.element import Element

with importlib.resources.as_file(
    importlib.resources.files("pyxtal") / "database" / "atomic_scattering_params.json"
) as path:
    ATOMIC_SCATTERING_PARAMS = loadfn(path)


class XRD:
    """
    A class to compute the powder XRD.

    Args:
        crystal: ase atoms object
        wavelength (float): wavelength of the X-ray in Angstrom (default: 1.54184)
        thetas (list): list of 2theta angles in degrees (default: [0, 180])
        res (float): resolution of the XRD in degrees (default: 0.01)
        per_N (int): number of hkl per process (default: 30000)
        ncpu: int, number of cpu to use (default: 1)
        preferred_orientation: boolean, whether to use preferred orientation
        march_parameter: float, the march parameter for preferred orientation
        TWO_THETA_TOL: tolerance to find repeating angles
        SCALED_INTENSITY_TOL: threshold for intensities
    """

    def __init__(
        self,
        crystal,
        wavelength=1.54184,
        thetas=[0, 180],
        res=0.01,
        per_N=30000,
        ncpu=1,
        filename=None,
        preferred_orientation=False,
        march_parameter=None,
        TWO_THETA_TOL=1e-5,
        SCALED_INTENSITY_TOL=1e-5,
    ):

        self.res = np.radians(res)
        if filename is None:
            self.wavelength = wavelength
            self.min2theta = np.radians(thetas[0])
            self.max2theta = np.radians(thetas[1])
            self.per_N = per_N
            self.ncpu = ncpu
            self.name = crystal.get_chemical_formula()
            self.preferred_orientation = preferred_orientation
            self.debye_waller_factor = 1.0  # default no debye waller factor
            self.march_parameter = march_parameter
            self.SCALED_INTENSITY_TOL = SCALED_INTENSITY_TOL
            self.TWO_THETA_TOL = TWO_THETA_TOL
            self.all_dhkl(crystal)
            self.skip_hkl = self.intensity(crystal)
            self.pxrdf()
        else:
            self.load(filename)

    def save(self, filename):
        """
        savetxt file
        """
        header = f"wavelength/thetas {self.wavelength:12.6f} {np.degrees(self.min2theta):6.2f} {np.degrees(self.max2theta):6.2f}"
        np.savetxt(filename, self.pxrd, header=header)

    def load(self, filename):
        """
        Load the pxrd from txt file
        """
        fp = open(filename)
        tmp = fp.readline()
        res = tmp.split()[2:]
        self.wavelength = float(res[0])
        self.min2theta = np.radians(float(res[1]))
        self.max2theta = np.radians(float(res[2]))

        pxrd = np.loadtxt(filename)
        self.theta2 = pxrd[:, 0]
        self.d_hkls = pxrd[:, 1]
        self.xrd_intensity = pxrd[:, -1]
        hkl_labels = []
        for i in range(len(pxrd)):
            h, k, l = int(pxrd[i, 2]), int(pxrd[i, 3]), int(pxrd[i, 4])
            hkl_labels.append([{"hkl": (h, k, l), "multiplicity": 1}])
        self.hkl_labels = hkl_labels
        self.pxrd = pxrd
        self.name = filename

    def __str__(self):
        return self.by_hkl()

    def __repr__(self):
        return str(self)

    def by_hkl(self, hkl=None, N_max=None):
        """
        d for any give abitray [h,k,l] index
        """
        s = ""
        if hkl is None:
            id1 = self.hkl_labels
            seqs = range(len(id1))
            if N_max is not None:
                seqs = range(min(N_max, len(id1)))
        else:
            seqs = None
            for id, label in enumerate(self.hkl_labels):
                hkl0 = list(label[0]["hkl"])  # label['multiplicity']
                if hkl == hkl0:
                    seqs = [id]

        if seqs is not None:
            s += "  2theta     d_hkl     hkl       Intensity  Multi\n"
            for i in seqs:
                s += f"{self.theta2[i]:8.3f}  {self.d_hkls[i]:8.3f}   "
                s += "[{:2d} {:2d} {:2d}]".format(*self.hkl_labels[i][0]["hkl"])
                s += f" {100 * self.xrd_intensity[i] / max(self.xrd_intensity):8.2f} "
                s += "{:8d}\n".format(self.hkl_labels[i][0]["multiplicity"])
        else:
            s += "This hkl is not in the given 2theta range"

        return s

    def all_dhkl(self, crystal):
        """
        3x3 representation -> 1x6 (a, b, c, alpha, beta, gamma)
        """
        rec_matrix = crystal.cell.reciprocal()
        eps = 1e-8  # small value to avoid division by zero
        d_max = self.wavelength / (np.sin(self.min2theta / 2) + eps) / 2
        d_min = self.wavelength / (np.sin(self.max2theta / 2) + eps) / 2

        # This block is to find the shortest d_hkl
        hkl_index = create_index()  # 2, 2, 2)
        hkl_max = np.array([1, 1, 1])

        for index in hkl_index:
            d = np.linalg.norm(np.dot(index, rec_matrix))
            multiple = int(np.ceil(1 / d / d_min))
            index *= multiple
            for i in range(len(hkl_max)):
                if hkl_max[i] < index[i]:
                    hkl_max[i] = index[i]

        h1, k1, l1 = hkl_max
        h = np.arange(-h1, h1 + 1)
        k = np.arange(-k1, k1 + 1)
        l = np.arange(-l1, l1 + 1)

        hkl = np.array(np.meshgrid(h, k, l)).transpose()
        hkl_list = np.reshape(hkl, [len(h) * len(k) * len(l), 3])
        hkl_list = hkl_list[np.where(hkl_list.any(axis=1))[0]]
        d_hkl = 1 / np.linalg.norm(np.dot(hkl_list, rec_matrix), axis=1)

        shortlist = np.where((d_hkl >= d_min) & (d_hkl < d_max))[0]
        d_hkl = d_hkl[shortlist]
        hkl_list = hkl_list[shortlist]
        sintheta = self.wavelength / 2 / d_hkl

        self.theta = np.arcsin(sintheta)#; print(self.theta[0:5000:])#; import sys; sys.exit()
        self.hkl_list = np.array(hkl_list, dtype=int)#; print(self.hkl_list[0:5000:]); import sys; sys.exit()
        self.d_hkl = d_hkl

    def intensity(self, crystal):
        """
        This function calculates all that is necessary to find the intensities.
        This scheme is similar to pymatgen
        If the number of hkl is significanly large, will automtically switch to
        the fast mode in which we only calculate the intensity and do not care
        the exact hkl families

        Args:

        """
        # obtain scattering parameters, atomic numbers, and occus
        # print("total number of hkl lists", len(self.hkl_list))
        # print("total number of coords:", len(crystal.get_scaled_positions()))
        # from time import time
        # t0 = time()

        N_atom, N_hkls = len(crystal), len(self.hkl_list)
        # Make sure the code don't split into too many cycles
        if self.per_N < N_atom:
            self.per_N = N_atom

        coeffs = np.zeros([N_atom, 4, 2])
        zs = np.zeros([N_atom, 1], dtype=int)
        for i, elem in enumerate(crystal.get_chemical_symbols()):
            if elem == "D":
                elem = "H"
            coeffs[i, :, :] = ATOMIC_SCATTERING_PARAMS[elem]
            zs[i] = Element(elem).z

        # A heavy calculation, Partition it to prevent the memory issue
        #s2s = self.d_hkl**2 #(np.sin(self.theta) / self.wavelength) ** 2  # M
        s2s = 1 / (4 * self.d_hkl ** 2)  # M
        hkl_per_proc = int(self.per_N / N_atom)
        N_cycle = int(np.ceil(N_hkls / hkl_per_proc))
        positions = crystal.get_scaled_positions()

        if self.ncpu == 1:
            N_cycles = range(N_cycle)
            Is = get_all_intensity(N_cycles, N_atom, self.per_N, positions, self.hkl_list, s2s, coeffs, zs)
        else:
            import multiprocessing as mp
            queue = mp.Queue()
            cycle_per_cpu = int(np.ceil(N_cycle / self.ncpu))
            processes = []
            for cpu in range(self.ncpu):
                N1 = cpu * cycle_per_cpu
                Start = N1 * hkl_per_proc
                if cpu + 1 == self.ncpu:
                    N2 = N_cycle
                    End = N_hkls
                else:
                    N2 = (cpu + 1) * cycle_per_cpu
                    End = N2 * hkl_per_proc

                cycles = range(N1, N2)

                p = mp.Process(
                    target=get_all_intensity_par,
                    args=(cpu, queue, cycles, Start, End,
                        hkl_per_proc, positions,
                        self.hkl_list[Start:End],
                        s2s[Start:End],
                        coeffs, zs))
                p.start()
                processes.append(p)

            unsorted_result = [queue.get() for p in processes]
            for p in processes:
                p.join()

            # collect results
            Is = np.zeros([N_hkls])
            for t in unsorted_result:
                Is[t[1] : t[2]] += t[3]

        # Lorentz polarization factor
        lfs = (1 + np.cos(2 * self.theta) ** 2) / (np.sin(self.theta) ** 2 * np.cos(self.theta))

        # Preferred orientation factor
        if self.preferred_orientation is not False:
            G = self.march_parameter
            pos = ((G * np.cos(self.theta)) ** 2 + 1 / G * np.sin(self.theta) ** 2) ** (-3 / 2)
        else:
            pos = np.ones(N_hkls)

        # Group the peaks by theta values
        _two_thetas = np.degrees(2 * self.theta)
        self.peaks = {}

        N = int((self.max2theta - self.min2theta) / self.res)
        if len(self.hkl_list) > N:
            skip_hkl = True
            refs = np.degrees(np.linspace(self.min2theta, self.max2theta, N + 1))
            dtol = np.degrees(self.res / 2)
            for ref_theta in refs:
                ids = np.where(np.abs(_two_thetas - ref_theta) < dtol)[0]
                if len(ids) > 0:
                    intensity = np.sum(Is[ids] * lfs[ids] * pos[ids])
                    self.peaks[ref_theta] = [
                        intensity,
                        self.hkl_list[ids],
                        self.d_hkl[ids[0]],
                    ]
        else:
            skip_hkl = False
            two_thetas = []
            for id in range(len(self.hkl_list)):
                hkl, d_hkl = self.hkl_list[id], self.d_hkl[id]
                # find where the scattered angles are equal
                ind = np.where(np.abs(np.subtract(two_thetas, _two_thetas[id])) < self.TWO_THETA_TOL)
                if len(ind[0]) > 0:
                    # append intensity, hkl plane, and thetas to lists
                    self.peaks[two_thetas[ind[0][0]]][0] += Is[id] * lfs[id] * pos[id]
                    self.peaks[two_thetas[ind[0][0]]][1].append(tuple(hkl))
                else:
                    self.peaks[_two_thetas[id]] = [
                        Is[id] * lfs[id] * pos[id],
                        [tuple(hkl)],
                        d_hkl,
                    ]
                    two_thetas.append(_two_thetas[id])

        # obtain important intensities (defined by SCALED_INTENSITY_TOL)
        # and corresponding 2*theta, hkl plane + multiplicity, and d_hkl
        max_intensity = max([v[0] for v in self.peaks.values()])
        x = []
        y = []
        hkls = []
        d_hkls = []
        count = 0
        for k in sorted(self.peaks.keys()):
            count += 1
            v = self.peaks[k]
            #if skip_hkl:
            #    fam = {}
            #    fam[tuple(v[1][0])] = len(v[1])
            #else:
            fam = self.get_unique_families(v[1])#; print(v[1], fam)
            if v[0] / max_intensity * 100 > self.SCALED_INTENSITY_TOL:
                # print(k, v[0]/max_intensity)
                x.append(k)
                y.append(v[0])

                hkls.append([{"hkl": hkl, "multiplicity": mult} for hkl, mult in fam.items()])
                d_hkls.append(v[2])

        self.theta2 = x
        self.xrd_intensity = y
        self.hkl_labels = hkls
        self.d_hkls = d_hkls

        return skip_hkl

    def pxrdf(self):
        """
        Group the equivalent hkl planes together by 2\theta angle
        N*6 arrays, Angle, d_hkl, h, k, l, intensity
        """
        rank = range(len(self.theta2))  # np.argsort(self.theta2)
        PL = []
        last = 0
        for i in rank:
            if self.xrd_intensity[i] > 0.01:
                angle = self.theta2[i]
                if abs(angle - last) < 1e-4:
                    PL[-1][-1] += self.xrd_intensity[i]
                else:
                    PL.append(
                        [
                            angle,
                            self.d_hkls[i],
                            self.hkl_labels[i][0]["hkl"][0],
                            self.hkl_labels[i][0]["hkl"][1],
                            self.hkl_labels[i][0]["hkl"][2],
                            self.xrd_intensity[i],
                        ]
                    )
                last = angle

        PL = np.array(PL)
        PL[:, -1] = PL[:, -1] / max(PL[:, -1])
        self.pxrd = PL

    def get_unique_families(self, hkls, verbose=False):
        """
        Returns unique families of Miller indices. Families must be permutations
        of each other.

        Args:
            hkls ([h, k, l]): List of Miller indices.
            verbose (bool): Whether or not to print out information on families.

        Returns:
            {hkl: multiplicity}: A dict with unique hkl and multiplicity.
        """

        # TODO: Definitely speed it up.
        def is_perm(hkl1, hkl2):
            h1 = np.abs(hkl1)
            h2 = np.abs(hkl2)
            return all(i == j for i, j in zip(sorted(h1), sorted(h2)))

        unique = collections.defaultdict(list)
        for hkl1 in hkls:
            found = False
            hkl1_tuple = tuple(hkl1)
            for hkl2 in unique:
                if is_perm(hkl1, hkl2):
                    found = True
                    unique[hkl2].append(hkl1_tuple)
                    break
            if not found:
                unique[hkl1_tuple].append(hkl1_tuple)

        pretty_unique = {}
        for v in unique.values():
            pretty_unique[sorted(v)[-1]] = len(v)

        return pretty_unique

    @staticmethod
    def draw_hkl(hkl):
        """
        turn negative numbers in hkl to overbar
        """

        hkl_str = []
        for i in hkl:
            if i < 0:
                label = str(int(-i))
                label = r"$\bar{" + label + "}$"
                hkl_str.append(str(label))
            else:
                hkl_str.append(str(int(i)))

        return hkl_str

    def plot_pxrd(
        self,
        filename=None,
        profile=None,
        minimum_I=0.01,
        show_hkl=True,
        fontsize=None,
        figsize=(20, 10),
        res=0.02,
        fwhm=0.1,
        ax=None,
        xlim=None,
        width=1.0,
        legend=None,
        show=False,
    ):
        """
        plot PXRD

        Args:
            filename (None): name of the xrd plot. If None, show the plot
            profile: type of peak profile
            minimum_I (0.01): the minimum intensity to include in the plot
            show_hkl (True): whether or not show hkl labels
            fontsize (None): fontsize of text in the plot
            figsize ((20, 10)): figsize
            xlim (None): the 2theta range [x_min, x_max]
        """
        import matplotlib as mpl
        import matplotlib.pyplot as plt

        if fontsize is not None:
            mpl.rcParams.update({"font.size": fontsize})

        if xlim is None:
            x_min, x_max = 0, np.degrees(self.max2theta)
        else:
            x_min, x_max = xlim[0], xlim[1]

        if ax is None:
            fig, axes = plt.subplots(1, 1, figsize=figsize)  # plt.figure(figsize=figsize)
            axes.set_title("PXRD of " + self.name)
        else:
            axes = ax

        if profile is None:
            dx = x_max - x_min
            for i in self.pxrd:
                axes.bar(i[0], i[-1], color="b", width=width * dx / 180)
                if i[-1] > minimum_I and x_min <= i[0] <= x_max and show_hkl:
                    label = self.draw_hkl(i[2:5])
                    axes.text(i[0] - dx / 40, i[-1], label[0] + label[1] + label[2])
        else:
            spectra = self.get_profile(method=profile, res=res, user_kwargs={"FWHM": fwhm})
            label = "Profile: " + profile if legend is None else legend
            axes.plot(spectra[0], spectra[1], label=label)
            axes.legend()

        axes.set_xlim([x_min, x_max])
        axes.set_xlabel(r"2$\Theta$ ($\lambda$=" + str(self.wavelength) + r" $\AA$)")
        axes.set_ylabel("Intensity")

        if ax is None:
            axes.grid()
            if filename is None:
                if show:
                    fig.show()
            else:
                fig.savefig(filename)
                # fig.close()
            return fig, axes
        return None

    def plotly_pxrd(
        self,
        profile="gaussian",
        minimum_I=0.01,
        res=0.02,
        FWHM=0.1,
        height=450,
        html=None,
    ):
        import plotly.graph_objects as go

        """
        interactive plot for pxrd powered by plotly
        Args:
            xrd: xrd object
            html: html filename (str)
        """

        x, y, labels = [], [], []
        for i in range(len(self.pxrd)):
            theta2, d, h, k, l, I = self.pxrd[i]
            h, k, l = int(h), int(k), int(l)
            if minimum_I < I:
                label = f"<br>2&#952;: {theta2:6.2f}<br>d: {d:6.4f}<br>"
                label += f"I: {I:6.4f}</br>hkl: ({h:d}{k:d}{l:d})"
                x.append(theta2)
                y.append(-0.1)
                labels.append(label)

        trace1 = go.Bar(
            x=x,
            y=y,
            text=labels,
            hovertemplate="%{text}",
            width=0.5,
            name="hkl indices",
        )
        if profile is None:
            fig = go.Figure(data=[trace1])
        else:
            spectra = self.get_profile(method=profile, res=res, user_kwargs={"FWHM": FWHM})
            trace2 = go.Scatter(x=spectra[0], y=spectra[1], name="Profile: " + profile)
            fig = go.Figure(data=[trace2, trace1])

        fig.update_layout(
            height=height,
            xaxis_title=f"2&#952; ({self.wavelength:.4f} &#8491;)",
            yaxis_title="Intensity",
            title="PXRD of " + self.name,
        )

        if os.environ.get("_", "").find("jupyter") == -1:
            if html is None:
                return fig.to_html()
            else:
                fig.write_html(html)
                return None
        else:
            print("This is running on Jupyter Notebook")
            return fig

    def get_profile(self, method="gaussian", res=0.01, user_kwargs=None):
        """
        return the profile detail
        """

        return Profile(method, res, user_kwargs).get_profile(
            self.theta2,
            self.xrd_intensity,
            np.degrees(self.min2theta),
            np.degrees(self.max2theta),
        )

    def get_plot(self, grainsize=20, orientation=0.1, thermo=0.1,
                 L=500, H=50, S=25, bg_order=6, bg_ratio=0.05,
                 mix_ratio=0.02, dx=0.02):
        """
        Generate a simulated XRD plot with various parameters.
        Inspired from Pysimxrd at PyPI.
        Needs to double check the parameters.

        Args:
            grainsize (float): Grain size in micrometers.
            orientation (float): Preferred orientation factor.
            thermo (float): Thermal vibration factor.
            L (float): Axial divergence length.
            H (float): Axial divergence height.
            S (float): Slit width.
            bg_order (int): Order of the polynomial background.
            bg_ratio (float): Ratio of background intensity.
            mix_ratio (float): Ratio of random noise intensity.
            dx (float): Step size for the simulated XRD.

        Returns:
            tuple: Simulated 2-theta values and corresponding intensities.
        """

        # Marked locations and intensities
        x, y = self.pxrd[:, 0], self.pxrd[:, -1] * 100
        thetas = np.radians(x/2)

        # Calculate FWHM using Scherrer equation
        # Standard Scherrer: FWHM_L = K*λ/(L*cosθ) where K≈0.9
        # For Lorentzian HWHM: γ = FWHM/2
        K = 0.9  # Scherrer constant (shape factor)
        fwhm = K * self.wavelength / (grainsize * np.cos(thetas) + 1e-10)  # FWHM in radians
        gamma = fwhm / 2  # Lorentzian HWHM

        # Convert HWHM to Gaussian variance for Voigt profile
        # For pure Gaussian: HWHM = sqrt(2*ln2) * σ
        # Therefore: σ² = HWHM² / (2*ln2)
        sigma2 = gamma ** 2 / (2 * np.log(2))  # Gaussian variance

        # Apply preferred orientation and Debye-Waller factor
        ori_m, ori_p = 1 - orientation, 1 + orientation
        ori = np.clip(np.random.normal(loc=1, scale=0.2), ori_m, ori_p)

        # Apply Debye-Waller factor
        deb = np.exp(-16/3 * np.pi**2 * thermo**2 * (np.sin(thetas) / self.wavelength)**2)
        y *= ori * deb
        #print(x, y, gamma, sigma2)

        # Get profiles
        theta_min, theta_max = np.degrees(self.min2theta), min(90.0, np.degrees(self.max2theta))
        x_sim = np.arange(theta_min, theta_max, dx)
        y_sim = np.zeros_like(x_sim)

        # Add each peak contribution
        for k in range(len(x)):
            if x[k] < theta_max:
                #print("Adding peak at 2theta =", x[k])
                y_sim += add_peak(x_sim, x[k], gamma[k], sigma2[k], L, H, S, dx) * y[k]

        # Add each peak contribution
        area = np.trapz(y_sim, x_sim)
        y_sim /= area#; print(area, y_sim.max())

        # Add background
        bg_coeffs = np.abs(np.random.randn(bg_order + 1))
        bg_coeffs[0] = -bg_coeffs[0]  # Ensure decreasing trend
        bg_fun = np.poly1d(bg_coeffs)
        #bg_fun = np.poly1d(np.random.randn(bg_order + 1))
        bg = bg_fun(x_sim)
        bg -= bg.min()
        bg_y = bg / bg.max() * y_sim.max() * bg_ratio
        mixture = np.random.uniform(0, y_sim.max() * mix_ratio, size=len(x_sim))
        y_sim +=  bg_y + mixture

        # Scale to (0, 100)
        y_min, y_max = y_sim.min(), y_sim.max()
        if y_max > y_min:  # Avoid division by zero
            y_sim = (y_sim - y_min) / (y_max - y_min) * 100
        else:
            y_sim = np.zeros_like(y_sim)
        #import matplotlib.pyplot as plt
        #plt.plot(x_sim, y_sim)
        #plt.show()
        return x_sim, y_sim

def add_peak(twotheta, mu, gamma, sigma2, L, H, S, step=0.02, width=0.1, sigma2_distor=0.001):
    """
    Add a single peak to the XRD pattern using Voigt profile,
    axial divergence, slit function, and lattice distortion.

    Args:
        twotheta (array-like): Array of 2-theta
        mu (float): Peak center (2-theta) in degrees.
        gamma (float): Lorentzian HWHM parameter.
        sigma2 (float): Gaussian variance parameter.
        L (float): Axial divergence length.
        H (float): Axial divergence height.
        S (float): Slit half-width.
        step (float): Step size for the 2-theta array.
        width (float): Width of the slit function in degrees.
        sigma2_distor (float): Variance for lattice distortion Gaussian.

    Returns:
        ndarray: Array of same shape as twotheta with the peak intensity.
    """
    # Determine l_gap based on mu value
    if mu <= 10:
        l_gap = 7.8
    elif 10 < mu <= 15:
        l_gap = 10
    elif 15 < mu <= 20:
        l_gap = 15
    elif 20 < mu <= 30:
        l_gap = 20
    else:
        l_gap = 30

    # Ensure mu-l_gap and mu+l_gap are recorded in twotheta or its extension
    x = np.arange(np.round(mu - l_gap, 2), np.round(mu + l_gap, 2), step)

    # Voigt profile calculation
    sigma = np.sqrt(sigma2)
    z = ((x - mu) + 1j * gamma) / (sigma * np.sqrt(2))
    voigt = np.real(wofz(z) / (sigma * np.sqrt(2 * np.pi)))

    # Axial divergence calculation
    axial = axial_div(x, mu, L, H, S)

    # Slit function calculation
    height = 1.0 / width
    slit = np.where((x >= mu - width / 2) & (x <= mu + width / 2), height, 0)

    #slit = np.zeros_like(x)
    #mask = np.abs(x - mu) <= width / 2
    #slit[mask] = 1.0 / width  # Normalized rectangular function

    # Lattice distortion calculation
    sigma = np.sqrt(sigma2_distor)
    distor = (1 / (sigma * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((x - mu) / sigma)**2)

    # Convolve the peaks
    combined = np.convolve(voigt, axial, mode='same')
    combined = np.convolve(combined, slit, mode='same')
    combined = np.convolve(combined, distor, mode='same')
    if np.sum(combined) > 0:
        combined /= np.sum(combined) * step  # Normalize peak and apply weight
        return map_intensity(combined, x, twotheta)
    else:
        return np.zeros_like(twotheta)

def axial_div_bak(x, mu, L, H, S):
    """
    Calculate the axial divergence peak contribution using the Van Laar model.

    Args:
        x (array-like): Array of 2-theta values in degrees.
        mu (float): Peak center (2-theta) in degrees.
        L (float): Axial divergence length (same units as H and S).
        H (float): Axial divergence height.
        S (float): Slit half-width (same units as H).

    Returns:
        ndarray: Array of same shape as x with the axial divergence shape (unnormalized).
    """
    axial_divergence = np.zeros_like(x)  # Initialize axial_divergence to zeros
    valid_indices = x <= mu  # Identify valid indices where x <= mu
    x_valid = np.radians(x[valid_indices])  # Get valid x values

    h = L * np.sqrt((np.cos(x_valid) / np.cos(np.radians(mu)))**2 - 1) + 1e-10 # Calculate h
    W = np.where((H - S <= h) & (h <= H + S), H + S - h, 0)  # Calculate W for valid h
    axial_divergence[valid_indices] = L / (2 * H * S * h * np.cos(x_valid)) * W
    #print('debug axial_div', mu, x_valid[-1], axial_divergence[valid_indices].max())
    axial_divergence /= (axial_divergence.max() + 1e-10)# in case numerical err
    cdf = np.zeros_like(x)
    mask = x < mu
    cdf[mask] = np.cumsum(axial_divergence[mask])
    return cdf

def axial_div(x, mu, L, H, S):
    """
    Van Laar axial divergence PDF (not CDF!)
    """
    x = np.asarray(x)
    f = np.zeros_like(x)

    mask = x < mu
    if not np.any(mask):
        return f

    x_m = np.radians(x[mask])
    mu_r = np.radians(mu)
    cos_mu = np.cos(mu_r)
    cos_x = np.cos(x_m)

    # Calculate h parameter with clipping to avoid negative square root
    cos_ratio_sq = (cos_x / cos_mu) ** 2
    h = L * np.sqrt(np.clip(cos_ratio_sq - 1, 0, None))

    # Window function: non-zero only when H - S <= h <= H + S
    W = np.clip(H + S - h, 0.0, 2 * S)

    # Van Laar axial divergence formula
    # Avoid division by zero
    denom = 2 * H * S * np.clip(h, 1e-10, None) * np.clip(cos_x, 1e-10, None)
    f[mask] = L * W / denom

    # Remove numerical noise and ensure non-negative
    f[~np.isfinite(f)] = 0.0
    f[f < 0] = 0.0

    # Normalize to unit area
    # Integrate only the non-zero part
    x_nonzero = x[mask]
    f_nonzero = f[mask]
    if len(x_nonzero) > 1 and np.sum(f_nonzero) > 0:
        area = np.trapz(f_nonzero, x_nonzero)
        if area > 0:
            f[mask] /= area
    return f

def map_intensity(peak, x, twotheta):
    """
    Map peak intensities from fine grid (x) to coarse grid (twotheta).
    Uses cubic spline interpolation to produce continuous, smooth profiles.

    Args:
        peak (array-like): Peak intensities on fine grid x.
        x (array-like): Fine grid positions (degrees).
        twotheta (array-like): Coarse grid positions (degrees).

    Returns:
        ndarray: Interpolated intensities on twotheta grid.
    """
    # Handle edge cases
    if len(peak) == 0 or len(x) == 0:
        return np.zeros_like(twotheta)

    # Check if x is monotonically increasing
    if not np.all(np.diff(x) > 0):
        # Sort by x if not already sorted
        sort_idx = np.argsort(x)
        x = x[sort_idx]
        peak = peak[sort_idx]

    # Use cubic spline interpolation for smooth results
    kind = 'linear' if len(x) < 4 else 'cubic'

    try:
        f_interp = interp1d(x, peak, kind=kind, bounds_error=False,
                           fill_value=0.0, assume_sorted=True)
        y_twotheta = f_interp(twotheta)
        # Ensure non-negative intensities
        y_twotheta[y_twotheta < 0] = 0.0
        return y_twotheta
    except (ValueError, RuntimeError) as e:
        print(f"Interpolation failed: {e}. Falling back to nearest-neighbor.")
        return _map_int_nearest_neighbor(peak, x, twotheta)


def _map_int_nearest_neighbor(peak, x, twotheta):
    """
    Fallback nearest-neighbor mapping when interpolation fails.

    Args:
        peak (array-like): Peak intensities on fine grid.
        x (array-like): Fine grid positions.
        twotheta (array-like): Coarse grid positions.

    Returns:
        ndarray: Nearest-neighbor intensities.
    """
    y_twotheta = np.zeros_like(twotheta, dtype=float)

    for x_val, peak_val in zip(x, peak):
        idx = np.argmin(np.abs(twotheta - x_val))
        if y_twotheta[idx] == 0:  # Only assign if not already set
            y_twotheta[idx] = peak_val
        else:
            y_twotheta[idx] += peak_val * 0.5  # Average with existing value

    return y_twotheta

# ----------------------------- Profile functions ------------------------------
class Profile:
    """
    This class applies a profiling function to simulated or
    experimentally obtained XRD spectra.

    Args:
        method (str): Type of function used to profile
        res (float): resolution of the profiling array in degree
        user_kwargs (dict): The parameters for the profiling method.
    """

    def __init__(self, method="mod_pseudo-voigt", res=0.02, user_kwargs=None):
        self.method = method
        self.user_kwargs = user_kwargs
        self.res = res
        kwargs = {}

        if method == "mod_pseudo-voigt":
            _kwargs = {
                "U": 5.776410e-03,
                "V": -1.673830e-03,
                "W": 5.668770e-03,
                "A": 1.03944,
                "eta_h": 0.504656,
                "eta_l": 0.611844,
            }
        elif method in ["gaussian", "lorentzian", "pseudo-voigt"]:
            _kwargs = {"FWHM": 0.1}

        else:
            msg = method + " isn't supported."
            raise NotImplementedError(msg)

        kwargs.update(_kwargs)

        if user_kwargs is not None:
            kwargs.update(user_kwargs)

        self.kwargs = kwargs

    def get_profile(self, two_thetas, intensities, min2theta, max2theta):
        """
        Performs profiling with selected function, resolution, and parameters

        Args:
            - two_thetas: 1d float array simulated/measured 2 theta values
            - intensities: simulated/measures peaks
        """
        N = int((max2theta - min2theta) / self.res)
        px = np.linspace(min2theta, max2theta, N)
        py = np.zeros(N)
        for two_theta, intensity in zip(two_thetas, intensities):
            # print(two_theta, intensity)
            if self.method == "gaussian":
                fwhm = self.kwargs["FWHM"]
                bin_edges = np.concatenate([px - self.res/2, [px[-1] + self.res/2]])
                tmp = np.zeros_like(px)
                for i in range(len(px)):
                    left, right = bin_edges[i], bin_edges[i+1]
                    tmp[i] = gaussian_integrated(left, right, two_theta, fwhm)
                #dtheta2 = ((px - two_theta) / fwhm) ** 2
                #tmp = np.exp(-4 * np.log(2) * dtheta2)
                # tmp = gaussian(two_theta, px, fwhm)

            elif self.method == "lorentzian":
                fwhm = self.kwargs["FWHM"]
                tmp = lorentzian(two_theta, px, fwhm)

            elif self.method == "pseudo-voigt":
                try:
                    fwhm_g = self.kwargs["FWHM-G"]
                    fwhm_l = self.kwargs["FWHM-L"]
                except:
                    fwhm_g = self.kwargs["FWHM"]
                    fwhm_l = self.kwargs["FWHM"]

                fwhm = (
                    fwhm_g**5
                    + 2.69269 * fwhm_g**4 * fwhm_l
                    + 2.42843 * fwhm_g**3 * fwhm_l**2
                    + 4.47163 * fwhm_g**2 * fwhm_l**3
                    + 0.07842 * fwhm_g * fwhm_l**4
                    + fwhm_l**5
                ) ** (1 / 5)
                eta = 1.36603 * fwhm_l / fwhm - 0.47719 * (fwhm_l / fwhm) ** 2 + 0.11116 * (fwhm_l / fwhm) ** 3
                tmp = pseudo_voigt(two_theta, px, fwhm, eta)

            elif self.method == "mod_pseudo-voigt":
                U = self.kwargs["U"]
                V = self.kwargs["V"]
                W = self.kwargs["W"]
                A = self.kwargs["A"]
                eta_h = self.kwargs["eta_h"]
                eta_l = self.kwargs["eta_l"]

                fwhm = np.sqrt(
                    U * np.tan(np.pi * two_theta / 2 / 180) ** 2 + V * np.tan(np.pi * two_theta / 2 / 180) + W
                )
                x = px - two_theta
                tmp = mod_pseudo_voigt(x, fwhm, A, eta_h, eta_l, N)

            py += intensity * tmp
            # print(intensity * tmp)

        py /= np.max(py)

        self.spectra = np.vstack((px, py))
        return self.spectra


# ------------------------------ Similarity between two XRDs ---------------------------------
class Similarity:
    def __init__(self, f, g, N=None, x_range=None, l=2.0, weight="cosine"):
        """
        Class to compute the similarity between two diffraction patterns
        Args:

            f: spectra1 (2D array)
            g: spectra2 (2D array)
            N: number of sampling points for the processed spectra
            x_range: the range of x values used to compute similarity ([x_min, x_max])
            l: cutoff value for shift (real)
            weight: weight function 'triangle' or 'cosine' (str)
        """

        fx, fy = f[0], f[1]
        gx, gy = g[0], g[1]
        self.l = abs(l)
        res1 = (fx[-1] - fx[0]) / len(fx)
        res2 = (gx[-1] - gx[0]) / len(gx)
        self.resolution = min([res1, res2]) / 3  # improve the resolution

        if N is None:
            self.N = int(2 * self.l / self.resolution)
        else:
            self.N = N
        self.r = np.linspace(-self.l, self.l, self.N)

        if x_range is None:
            x_min = max(np.min(fx), np.min(gx))
            x_max = min(np.max(fx), np.max(gx))
        else:
            x_min, x_max = x_range[0], x_range[1]

        self.x_range = [x_min, x_max]

        f_inter = interp1d(fx, fy, "cubic", fill_value="extrapolate")
        g_inter = interp1d(gx, gy, "cubic", fill_value="extrapolate")

        fgx_new = np.linspace(x_min, x_max, int((x_max - x_min) / self.resolution) + 1)
        fy_new = f_inter(fgx_new)
        gy_new = g_inter(fgx_new)

        self.fx, self.gx, self.fy, self.gy = fgx_new, fgx_new, fy_new, gy_new
        self.weight = weight
        if self.weight == "triangle":
            w = self.triangleFunction()
        elif self.weight == "cosine":
            w = self.cosineFunction()
        else:
            msg = self.weight + "is not supported"
            raise NotImplementedError(msg)

        Npts = len(self.fx)
        d = self.fx[1] - self.fx[0]

        self.value = similarity_calculate(self.r, w, d, Npts, self.fy, self.gy)

    def __str__(self):
        return f"The similarity between two PXRDs is {self.value:.4f}"

    def __repr__(self):
        return str(self)

    def triangleFunction(self):
        """
        Triangle function to weight correlations
        """
        w = 1 - np.abs(self.r / self.l)
        ids = np.abs(self.r) > self.l
        w[ids] = 0

        return w

    def cosineFunction(self):
        """
        cosine function to weight correlations
        """

        w = 0.5 * (np.cos(np.pi * self.r / self.l) + 1.0)
        ids = np.abs(self.r) > self.l
        w[ids] = 0

        return w

    def show(self, filename=None, fontsize=None, labels=None):
        """
        show the comparison plot

        Args:
            filename (None): name of the xrd plot. If None, show the plot
            labels [A, B]: labels of each plot
        """
        import matplotlib as mpl
        import matplotlib.pyplot as plt

        if labels is None:
            labels = ["profile 1", "profile 2"]
        if fontsize is not None:
            mpl.rcParams.update({"font.size": fontsize})

        fig1 = plt.figure(1, figsize=(15, 6))
        fig1.add_axes((0.1, 0.3, 0.8, 0.6))

        plt.plot(self.fx, self.fy, label=labels[0])
        plt.plot(self.fx, -self.gy, label=labels[1])
        plt.legend()

        # Residual plot
        residuals = self.gy - self.fy
        fig1.add_axes((0.1, 0.1, 0.8, 0.2))
        plt.plot(self.gx, residuals, ".r", markersize=0.5)
        plt.title(f"{self.value:6f}")

        if filename is None:
            plt.show()
        else:
            plt.savefig(filename)
            plt.close()


def mod_pseudo_voigt(x, fwhm, A, eta_h, eta_l, N):
    """
    A modified split-type pseudo-Voigt function for profiling peaks
    - Izumi, F., & Ikeda, T. (2000).
    """

    tmp = np.zeros(N)
    for xi, dx in enumerate(x):
        if dx < 0:
            A = A
            eta_l = eta_l
            eta_h = eta_h
        else:
            A = 1 / A
            eta_l = eta_h
            eta_h = eta_l

        tmp[xi] = (
            ((1 + A) * (eta_h + np.sqrt(np.pi * np.log(2)) * (1 - eta_h)))
            / (
                eta_l
                + np.sqrt(np.pi * np.log(2)) * (1 - eta_l)
                + A * (eta_h + np.sqrt(np.pi * np.log(2)) * (1 - eta_h))
            )
            * (
                eta_l * 2 / (np.pi * fwhm) * (1 + ((1 + A) / A) ** 2 * (dx / fwhm) ** 2) ** (-1)
                + (1 - eta_l)
                * np.sqrt(np.log(2) / np.pi)
                * 2
                / fwhm
                * np.exp(-np.log(2) * ((1 + A) / A) ** 2 * (dx / fwhm) ** 2)
            )
        )
    return tmp


def gaussian(theta2, alpha, fwhm):
    """
    Gaussian function for profiling peaks
    """

    tmp = ((alpha - theta2) / fwhm) ** 2
    return np.exp(-4 * np.log(2) * tmp)


def lorentzian(theta2, alpha, fwhm):
    """
    Lorentzian function for profiling peaks
    """

    tmp = 1 + 4 * ((alpha - theta2) / fwhm) ** 2
    return 1 / tmp


def pseudo_voigt(theta2, alpha, fwhm, eta):
    """
    Original Pseudo-Voigt function for profiling peaks
    - Thompson, D. E. Cox & J. B. Hastings (1986).
    """

    L = lorentzian(theta2, alpha, fwhm)
    G = gaussian(theta2, alpha, fwhm)
    return eta * L + (1 - eta) * G


def similarity_calculate(r, w, d, Npts, fy, gy):
    """
    Compute the similarity between the pair of spectra f, g
    """

    xCorrfg_w, aCorrff_w, aCorrgg_w = 0, 0, 0
    for r0, w0 in zip(r, w):
        Corrfg, Corrff, Corrgg = 0, 0, 0
        shift = int(r0 / d)
        for i in range(Npts):
            if 0 <= i + shift <= Npts - 1:
                Corrfg += fy[i] * gy[i + shift] * d
                Corrff += fy[i] * fy[i + shift] * d
                Corrgg += gy[i] * gy[i + shift] * d

        xCorrfg_w += w0 * Corrfg * d
        aCorrff_w += w0 * Corrff * d
        aCorrgg_w += w0 * Corrgg * d

    return np.abs(xCorrfg_w / np.sqrt(aCorrff_w * aCorrgg_w))


def create_index(imax=1, jmax=1, kmax=1):
    """
    shortcut to get the index
    """
    hkl_index = []
    for i in range(-imax, imax + 1):
        for j in range(-jmax, jmax + 1):
            for k in range(-kmax, kmax + 1):
                hkl = np.array([i, j, k])
                if sum(hkl * hkl) > 0:
                    hkl_index.append(hkl)
    return np.array(hkl_index).reshape([len(hkl_index), 3])


def get_intensity(positions, hkl, s2, coeffs, z):
    """
    Calculate the intensity for a given set of positions, hkl, s2, coefficients, and atomic numbers.
    Args:
        positions (np.ndarray): N*3 array of atomic positions in fractional coordinates.
        hkl (np.ndarray): 3*M array of Miller indices.
        s2 (np.ndarray): M array of sin^2(theta) values.
        coeffs (np.ndarray): N*4*2 array of coefficients for each atom.
        z (np.ndarray): N*1 array of atomic numbers.

    Returns:
        np.ndarray: M array of calculated intensities.
    """
    #N = len(positions); positions = np.random.rand(N, 3) #-= np.round(positions)  # ensure within [0,1)
    g_dot_rs = np.dot(positions, hkl)  # N*3 dot 3*M -> N*M
    exps = np.exp(-2j * np.pi * g_dot_rs)  # N*M
    #print(exps[0]); import sys; sys.exit()
    tmp1 = np.exp(np.einsum("ij,k->ijk", -coeffs[:, :, 1], s2))  # N*4, M
    tmp2 = np.einsum("ij,ijk->ik", coeffs[:, :, 0], tmp1)  # N*4, N*M
    sfs = np.add(-41.78214 * np.einsum("ij,j->ij", tmp2, s2), z)  # N*M, M -> N*M
    fs = np.sum(sfs * exps, axis=0)  # M

    # Final intensity values, M
    #for i, f in enumerate(fs):
    #    if f.real > 0.1 and s2[i] < 0.05:
    #        print("hkl", hkl[:, i], s2[i], '|F|', f)

    return (fs * fs.conjugate()).real


def get_all_intensity(N_cycles, N_atom, per_N, positions, hkls, s2s, coeffs, zs):
    Is = np.zeros(len(hkls))
    for i, cycle in enumerate(N_cycles):
        N1 = int(per_N * (cycle) / N_atom)
        if i + 1 == len(N_cycles):
            N2 = min([len(hkls), int(per_N * (cycle + 1) / N_atom)])
        else:
            N2 = int(per_N * (cycle + 1) / N_atom)
        hkl, s2 = hkls[N1:N2].T, s2s[N1:N2]
        Is[N1:N2] = get_intensity(positions, hkl, s2, coeffs, zs)
    return Is


def get_all_intensity_par(cpu, queue, cycles, Start, End, hkl_per_proc, positions, hkls, s2s, coeffs, zs):
    Is = np.zeros(End - Start)
    for i, cycle in enumerate(cycles):
        N1 = cycle * hkl_per_proc - Start
        N2 = End - Start if i + 1 == len(cycles) else N1 + hkl_per_proc
        hkl, s2 = hkls[N1:N2].T, s2s[N1:N2]
        Is[N1:N2] = get_intensity(positions, hkl, s2, coeffs, zs)
        # print('run', cpu, N1+Start, N2+Start, N1, N2)
    queue.put((cpu, Start, End, Is))


def gaussian_integrated(bin_left, bin_right, center, fwhm):
    sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))
    # Integrate the normalized Gaussian over [bin_left, bin_right]
    return 0.5 * (erf((bin_right - center) / (np.sqrt(2) * sigma)) -
                  erf((bin_left - center) / (np.sqrt(2) * sigma)))

def pxrd_refine(xtal, ref_pxrd, thetas, steps=50):
    """
    Improve the lattice w.r.t the reference PXRD

    Args:
        xtal: pyxtal object
        ref_pxrd: tuple of (thetas, intensities) for the reference PXRD
        thetas: list of angles to calculate the PXRD
        steps (int): number of steps for optimization

    Returns:
        xtal: refined pyxtal object
        x: parameters used for optimization
        sim: similarity value between the refined PXRD and the reference PXRD
    """
    from scipy.optimize import minimize

    def fun(x0, rep, ref_pxrd, thetas):
        rep.x[0][1:] = x0
        s = rep.to_pyxtal()
        xrd = s.get_XRD(thetas=thetas)
        pxrd = xrd.get_profile(res=0.15, user_kwargs={"FWHM": 0.25})
        sim = Similarity(ref_pxrd, pxrd, x_range=thetas).value
        return -sim

    rep = xtal.get_1D_representation()
    x0 = rep.x[0][1:]
    f0 = fun(x0, rep, ref_pxrd, thetas)
    if f0 < -0.8:
        res = minimize(fun, x0,
                       args=(rep, ref_pxrd, thetas),
                       method="Nelder-Mead",
                       options={"maxiter": steps})
        rep.x[0][1:] = res.x
        xtal = rep.to_pyxtal()
        return xtal, -res.fun
    else:
        #print("The initial PXRD is unlikely to match the reference PXRD well.")
        return xtal, -f0


def check_pxrd_match(xtal, ref_pxrd, s_tol=0.8, top_n=3, peak_tol=0.1, ang_tol=1.0,
                     wave_length=1.5406, verbose=False):
    """
    Check if there is a false match between the pyxtal structure and the reference PXRD.
    First, check the similarity between the two PXRDs. If the similarity is above s_tol,
    Second, for each of the top_n strongest peaks in the computed PXRD, check if the related
    peaks (within peak_tol) are present in the reference PXRD within a tolerance.

    Args:
        xtal: pyxtal object
        ref_pxrd: a 2D array of (thetas, intensities) for the reference PXRD
        s_tol: similarity tolerance, default is 0.8
        top_n: number of strongest peaks to consider, default is 3
        peak_tol: tolerance for peaks to be considered for a comparison, default is 0.05
        ang_tol: tolerance for matching peaks in degrees, default is 1.0
        wave_length: X-ray wavelength, default is Cu K-alpha
        verbose: whether or not print the information

    Returns:
        bool: True if there is a false match, False otherwise
    """
    xrd = xtal.get_XRD(thetas=[ref_pxrd[0][0], ref_pxrd[0][-1]], wavelength=wave_length)
    pxrd = xrd.get_profile(res=0.15, user_kwargs={"FWHM": 0.25})
    sim = Similarity(ref_pxrd, pxrd, x_range=[ref_pxrd[0][0], ref_pxrd[0][-1]]).value
    if verbose:
        print(f"Similarity between computed PXRD and reference PXRD: {sim:.4f}")
        print(xrd)
    if sim > s_tol:
        # get the strongest peaks from the computed PXRD from xtal
        peaks = xrd.pxrd[:, -1] # intensity
        hkls = xrd.pxrd[:, 2:5] # hkl
        thetas = xrd.pxrd[:, 0] # 2theta
        sorted_indices = np.argsort(peaks)[::-1]
        sorted_hkls = hkls[sorted_indices]
        sorted_peaks = peaks[sorted_indices]
        sorted_thetas = thetas[sorted_indices]
        for i in range(top_n):
            hkl = sorted_hkls[i]
            # get the hkls that are related to the current hkl
            for j, h in enumerate(hkls):
                if is_multiple(h, hkl) and not np.all(h == hkl):
                    theta = sorted_thetas[j]
                    peak = sorted_peaks[j]
                    if peak > peak_tol:
                        if verbose:
                            print(f"Checking {h}/{hkl} in top {i+1} peak  => {peak:.2f} at {theta:.2f}")
                        # check if there is a peak in the reference PXRD within ang_tol and peak_tol
                        close_peaks = ref_pxrd[0][(ref_pxrd[0] >= theta - ang_tol) & (ref_pxrd[0] <= theta + ang_tol)]
                        close_peaks = close_peaks[np.abs(close_peaks - theta) <= peak_tol]
                        if len(close_peaks) == 0:
                            if verbose:
                                print(f"False match at hkl {hkl}/{h} at {theta:.2f} not in ref. PXRD")
                            return 0
        return sim
    else:
        return sim  # Similarity too low to consider

def is_multiple(hkl, ref_hkl):
    # Avoid division by zero and require ref_hkl is not (0,0,0)
    if np.all(ref_hkl == 0):
        return False
    # Find the scaling factor for each component, ignore zeros in ref_hkl
    factors = []
    for h, r in zip(hkl, ref_hkl):
        if r == 0:
            if h != 0:
                return False
        else:
            factors.append(h / r)
    # All nonzero factors must be equal and positive integer
    if len(factors) == 0:
        return False
    first = factors[0]
    if not np.allclose(factors, first):
        return False
    # Check if the factor is a positive integer
    return first > 0 and np.isclose(first, int(round(first)))

def get_para_from_pxrd(ref_pxd, spg, wave_length=1.5406):
    """
    Estimate the lattice parameters from the reference PXRD using Bragg's law and cubic assumption.

    Args:
        ref_pxd: tuple of (thetas, intensities) for the reference PXRD
        spg: space group number
        wave_length: X-ray wavelength, default is Cu K-alpha

    Returns:
        a: estimated lattice parameter
    """
    # Get the first peak position
    #thetas, intensities = ref_pxd
    #peak_index = np.argmax(intensities)
    #theta = thetas[peak_index] / 2  # Convert 2theta to theta
    #a = wave_length / (2 * np.sin(np.radians(theta)))  # Bragg's law
    #if spg > 194:
    #    if spg in [196, 202, 203, 209, 210, 216, 219, 225, 226, 227, 228]: # F-cubic 111
    #        a /= np.sqrt(3)
    #    elif spg in [195, 198, 199, 200, 201, 205, 206, 207, 208, 211, 212, 213, 214, 215]: # I-cubic 200
    #        a /= np.sqrt(2)
    #    cell = [a, a, a, 90, 90, 90]
    #elif 143 <= spg <= 194: # (001) or (100)
    #    cell = [[a, a, a, 90, 90, 120], []]
    #elif 75 <= spg <= 142:
    #    if P: # (100) or (001)
    #    elif I: # (101) or (110)
    #        a /= np.sqrt(2)
    #elif 16 <= spg <= 74: # (100) or (001)
    #    if P: #(100), (010), (001)
    #    elif I: (101) or (110)
    #    elif C/F/I: (001)/(020)/(101) ???
    #elif 3 <= spg <= 15: # (001) or (100)
    #    if P: #(100), (010), (001)
    #    elif A/B/C: (001)/(010)/(100)
    #else: # (001), (100), (010)
    #    cell = [a, a, 10, 90, 90, 90]
    #return cell
    pass

if __name__ == "__main__":
    from optparse import OptionParser
    from matplotlib import pyplot as plt
    from pyxtal import pyxtal
    from pyxtal.util import parse_cif
    from pymatgen.core import Structure

    parser = OptionParser()
    parser.add_option("-f", dest="cif", help="cif file name")
    parser.add_option("-r", dest="ref", help="ref pxrd file", default="ref_pxrd.txt")
    parser.add_option("-s", dest="step", type=int, default=30, help="steps, optional")
    parser.add_option("-o", dest="out", default="PXRD-Matched.cif", help="output")
    parser.add_option("--smin1", dest="smin1", type=float, default=0.85,
                      help="minimum similarity to refine PXRD, default 0.85")
    parser.add_option("--smin2", dest="smin2", type=float, default=0.92,
                      help="minimum similarity to refine PXRD, default 0.92")
    parser.add_option("--thetas", dest="thetas", type=float, nargs=2, default=[5, 35],
                      help="thetas for PXRD calculation, default [5, 35]")

    (options, args) = parser.parse_args()
    thetas = options.thetas
    if options.ref is None:
        raise ValueError("Please provide a reference PXRD file using -r option.")
    if not os.path.exists(options.ref):
        raise FileNotFoundError(f"Reference PXRD '{options.ref}' does not exist.")

    smiles = None
    with open(options.cif) as f:
        lines = f.readlines()
        for l in lines:
            if 'smile' in l:
                smile_str = l.split(':')[1].strip()
                smiles = [s + '.smi' for s in smile_str.split('.')]
                break
    if smiles is None:
        raise ValueError("No smiles found in the CIF file. Please check the CIF format.")
    else:
        print("Smiles found:", smiles)
    with open(options.out, 'w') as f: f.write(f'smiles: {smile_str}\n')

    if options.ref.endswith(".cif"):
        s = pyxtal(molecular=True)
        s.from_seed(options.ref, molecules=smiles)
        ref_pxrd =s.get_XRD(thetas=thetas).get_profile(res=0.15, user_kwargs={"FWHM": 0.25})
    else:
        data = np.loadtxt(options.ref, skiprows=1)  # Load the reference PXRD data
        data[:, 1] /= np.max(data[:, 1])  # Normalize the intensity
        ref_pxrd = (data[:, 0], data[:, 1])

    #cifs, engs = parse_cif(options.cif)#, eng=True)
    cifs = parse_cif(options.cif)#, eng=True)
    pxrds = []
    s = pyxtal(molecular=True)
    for i, cif in enumerate(cifs):
        try:
            pmg = Structure.from_str(cif, fmt="cif")
            s.from_seed(pmg, molecules=smiles)
        except:
            print(f"Failed to parse CIF {i}. Skipping...")
            continue
        pxrd1 = s.get_XRD(thetas=thetas).get_profile(res=0.15, user_kwargs={"FWHM": 0.25})
        s1, val1 = pxrd_refine(s, ref_pxrd, thetas, steps=0)
        print(i, s1.lattice, val1)

        if val1 > options.smin1:
            s2, val2 = pxrd_refine(s1, ref_pxrd, thetas, steps=options.step)
            pxrd2 = s2.get_XRD(thetas=thetas).get_profile(res=0.15, user_kwargs={"FWHM": 0.25})
            print(i, s2.lattice, val2, '++++++')

            if val2 > options.smin2:
                str1 = s1.lattice.__str__(fmt="4.1f", ltype=False)
                str2 = s2.lattice.__str__(fmt="4.1f", ltype=False)
                label1 = f"{str1} - {val1:.3f} - {s1.group.number}"
                label2 = f"{str2} - {val2:.3f} - {s2.group.number}"
                pxrds.append((pxrd1, label1, pxrd2, label2))
                with open(options.out, 'a+') as f:
                    f.writelines(s1.to_file(header=label1))
                    f.writelines(s2.to_file(header=label2))

                fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(12, 4))

                item1 = label1.split('-')
                item2 = label2.split('-')
                cell1, sim1, spg1 = item1[0], item1[1], item1[2].strip()
                cell2, sim2, spg2 = item2[0], item2[1], item2[2].strip()

                l1 = f"Init. Similarity: {sim1} ({spg1})"
                l2 = f"Opt.  Similarity: {sim2} ({spg2})"

                ax1.plot(ref_pxrd[0], ref_pxrd[1], 'black', label=l1, lw=1.0, alpha=0.5)
                ax2.plot(ref_pxrd[0], ref_pxrd[1], 'black', label=l2, lw=1.0, alpha=0.5)
                ax1.plot(pxrd1[0], pxrd1[1], 'b:', label=cell1, lw=1.2, alpha=0.5)
                ax2.plot(pxrd2[0], pxrd2[1], 'b:', label=cell2, lw=1.2, alpha=0.5)

                ax1.set_xlabel('2θ (degrees)')
                ax1.set_ylabel('Intensity (a.u.)')
                ax2.set_xlabel('2θ (degrees)')
                ax2.set_ylabel('Intensity (a.u.)')

                ax1.legend()
                ax2.legend()

                fig.suptitle(f'PXRD Match {i+1}')
                fig.tight_layout()

                fig.savefig(f"pxrd_match_{i+1}.png", dpi=150)
                plt.close(fig)



    if len(pxrds) == 0:
        raise ValueError("No PXRDs found that match the reference PXRD well.")
    else:
        print("PXRDs found:", len(pxrds))

    #print("Saving individual PXRD match figures...")
    #for i, (pxrd1, label1, pxrd2, label2) in enumerate(pxrds):
    #print(f"Saved {len(pxrds)} individual PXRD match figures.")
