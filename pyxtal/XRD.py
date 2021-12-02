"""
Module for XRD simulation (experimental stage)
"""
import os
import collections
import numpy as np
import numba as nb
from scipy.interpolate import interp1d
from monty.serialization import loadfn
from pkg_resources import resource_filename
from pyxtal.database.element import Element

ATOMIC_SCATTERING_PARAMS = loadfn(resource_filename("pyxtal", "database/atomic_scattering_params.json"))

class XRD():
    """
    a class to compute the powder XRD.

    Args:
        crystal: ase atoms object
        wavelength: float
        max2theta: float
        preferred_orientation: boolean
        march_parameter: float
    """

    def __init__(self, crystal, wavelength=1.54184,
                 thetas = [0, 180],
                 preferred_orientation = False,
                 march_parameter = None):


        self.wavelength = wavelength
        self.min2theta = np.radians(thetas[0])
        self.max2theta = np.radians(thetas[1])
        self.name = crystal.get_chemical_formula()
        self.preferred_orientation = preferred_orientation
        self.march_parameter = march_parameter
        self.all_dhkl(crystal)
        self.intensity(crystal)
        self.pxrdf()

    def __str__(self):
        return self.by_hkl()

    def __repr__(self):
        return str(self)

    def by_hkl(self, hkl=None):
        """
        d for any give abitray [h,k,l] index
        """
        s = ""
        if hkl is None:
            id1 = self.hkl_labels
            seqs = range(len(id1))
        else:
            seqs = None
            for id, label in enumerate(self.hkl_labels):
                hkl0 = list(label[0]['hkl']) #label['multiplicity']
                if hkl == hkl0:
                    seqs = [id]

        if seqs is not None:
            s += '  2theta     d_hkl     hkl       Intensity  Multi\n'
            for i in seqs:
                s += "{:8.3f}  {:8.3f}   ".format(self.theta2[i], self.d_hkls[i])
                s += "[{:2d} {:2d} {:2d}]".format(*self.hkl_labels[i][0]["hkl"])
                s += " {:8.2f} ".format(100*self.xrd_intensity[i]/max(self.xrd_intensity))
                s += "{:8d}\n".format(self.hkl_labels[i][0]["multiplicity"])
        else:
            s += 'This hkl is not in the given 2theta range'

        return s

    def all_dhkl(self, crystal):
        """
        3x3 representation -> 1x6 (a, b, c, alpha, beta, gamma)
        """

        #rec_matrix = crystal.get_reciprocal_cell()
        rec_matrix = crystal.cell.reciprocal()
        d_min = self.wavelength/np.sin(self.max2theta/2)/2

        # This block is to find the shortest d_hkl
        # for all basic directions (1,0,0), (0,1,0), (1,1,0), (1,-1,0)

        hkl_index = create_index()
        hkl_max = np.array([1,1,1])

        for index in hkl_index:
            d = np.linalg.norm(np.dot(index, rec_matrix))
            multiple = int(np.ceil(1/d/d_min))
            index *= multiple
            for i in range(len(hkl_max)):
                if hkl_max[i] < index[i]:
                    hkl_max[i] = index[i]

        h1, k1, l1 = hkl_max
        h = np.arange(-h1,h1+1)
        k = np.arange(-k1,k1+1)
        l = np.arange(-l1,l1+1)

        hkl = np.array((np.meshgrid(h,k,l))).transpose()
        hkl_list = np.reshape(hkl, [len(h)*len(k)*len(l),3])
        hkl_list = hkl_list[np.where(hkl_list.any(axis=1))[0]]
        d_hkl = 1/np.linalg.norm( np.dot(hkl_list, rec_matrix), axis=1)

        shortlist = d_hkl > (d_min)
        d_hkl = d_hkl[shortlist]
        hkl_list = hkl_list[shortlist]
        sintheta = self.wavelength/2/d_hkl

        self.theta = np.arcsin(sintheta)
        self.hkl_list = np.array(hkl_list)
        self.d_hkl = d_hkl

    def intensity(self, crystal):

        """
        This function calculates all that is necessary to find the intensities.
        This scheme is based off of pymatgen
        Needs improvement from different correction factors.
        """

        d0 = (1/2/self.d_hkl)**2

        # obtiain scattering parameters, atomic numbers, and occus (need to look into occus)
        coeffs = []
        zs = []

        for elem in crystal.get_chemical_symbols():
            if elem == 'D':
                elem = 'H'
            c = ATOMIC_SCATTERING_PARAMS[elem]
            z = Element(elem).z
            coeffs.append(c)
            zs.append(z)

        coeffs = np.array(coeffs)
        self.peaks = {}
        two_thetas = []

        # self.march_parameter = 1
        TWO_THETA_TOL = 1e-5 # tolerance to find repeating angles
        SCALED_INTENSITY_TOL = 1e-5 # threshold for intensities

        for hkl, s2, theta, d_hkl in zip(self.hkl_list, d0, self.theta, self.d_hkl):

            # calculate the scattering factor sf
            g_dot_r = np.dot(crystal.get_scaled_positions(), np.transpose([hkl])).T[0]
            sf = zs - 41.78214 * s2 * np.sum(coeffs[:, :, 0] * np.exp(-coeffs[:, :, 1] * s2), axis=1)

            # calculate the structure factor f
            f = np.sum(sf * np.exp(2j * np.pi * g_dot_r))

            # calculate the lorentz polarization factor lf
            lf = (1 + np.cos(2 * theta) ** 2) / (np.sin(theta) ** 2 * np.cos(theta))

            # calculate the preferred orientation factor
            if self.preferred_orientation != False:
                G = self.march_parameter
                po = ((G * np.cos(theta))**2 + 1/G * np.sin(theta)**2)**(-3/2)
            else:
                po = 1

            # calculate the intensity I
            I = (f * f.conjugate()).real

            # calculate 2*theta
            two_theta = np.degrees(2 * theta)

            # find where the scattered angles are equal
            ind = np.where(np.abs(np.subtract(two_thetas, two_theta)) < TWO_THETA_TOL)

            # append intensity, hkl plane, and thetas to lists
            if len(ind[0]) > 0:
                self.peaks[two_thetas[ind[0][0]]][0] += I * lf * po
                self.peaks[two_thetas[ind[0][0]]][1].append(tuple(hkl))
            else:
                self.peaks[two_theta] = [I * lf * po, [tuple(hkl)],d_hkl]
                two_thetas.append(two_theta)

        # obtain important intensities (defined by SCALED_INTENSITY_TOL)
        # and corresponding 2*theta, hkl plane + multiplicity, and d_hkl

        max_intensity = max([v[0] for v in self.peaks.values()])
        x = []
        y = []
        hkls = []
        d_hkls = []
        count = 0
        for k in sorted(self.peaks.keys()):
            count +=1
            v = self.peaks[k]
            fam = self.get_unique_families(v[1])
            if v[0] / max_intensity * 100 > SCALED_INTENSITY_TOL:
                x.append(k)
                y.append(v[0])

                hkls.append([{"hkl": hkl, "multiplicity": mult}
                             for hkl, mult in fam.items()])
                d_hkls.append(v[2])

        self.theta2 = x
        self.xrd_intensity = y
        self.hkl_labels = hkls
        self.d_hkls = d_hkls

    def pxrdf(self):
        """
        Group the equivalent hkl planes together by 2\theta angle
        N*6 arrays, Angle, d_hkl, h, k, l, intensity
        """

        rank = range(len(self.theta2)) #np.argsort(self.theta2)
        PL = []
        last = 0
        for i in rank:
            if self.xrd_intensity[i] > 0.01:
                angle = self.theta2[i]
                if abs(angle-last) < 1e-4:
                    PL[-1][-1] += self.xrd_intensity[i]
                else:
                    PL.append([angle, self.d_hkls[i], \
                             self.hkl_labels[i][0]["hkl"][0], \
                             self.hkl_labels[i][0]["hkl"][1], \
                             self.hkl_labels[i][0]["hkl"][2], \
                             self.xrd_intensity[i]])
                last = angle

        PL = (np.array(PL))
        PL[:,-1] = PL[:,-1]/max(PL[:,-1])
        self.pxrd = PL

    def get_unique_families(self,hkls):
        """
        Returns unique families of Miller indices. Families must be permutations
        of each other.
        Args:
            hkls ([h, k, l]): List of Miller indices.
        Returns:
            {hkl: multiplicity}: A dict with unique hkl and multiplicity.
        """

       # TODO: Definitely can be sped up.
        def is_perm(hkl1, hkl2):
            h1 = np.abs(hkl1)
            h2 = np.abs(hkl2)
            return all([i == j for i, j in zip(sorted(h1), sorted(h2))])

        unique = collections.defaultdict(list)
        for hkl1 in hkls:
            found = False
            for hkl2 in unique.keys():
                if is_perm(hkl1, hkl2):
                    found = True
                    unique[hkl2].append(hkl1)
                    break
            if not found:
                unique[hkl1].append(hkl1)

        pretty_unique = {}
        for k, v in unique.items():
            pretty_unique[sorted(v)[-1]] = len(v)

        return pretty_unique

    @staticmethod
    def draw_hkl(hkl):
        """
        turn negative numbers in hkl to overbar
        """

        hkl_str= []
        for i in hkl:
            if i<0:
                label = str(int(-i))
                label = r"$\bar{" + label + '}$'
                hkl_str.append(str(label))
            else:
                hkl_str.append(str(int(i)))

        return hkl_str

    def plot_pxrd(self, filename=None, minimum_I=0.01, show_hkl=True,\
                fontsize=None, figsize=(20,10), xlim=None, width=1.0):
        """
        plot PXRD

        Args:
            filename (None): name of the xrd plot. If None, show the plot
            minimum_I (0.01): the minimum intensity to include in the plot
            show_hkl (True): whether or not show hkl labels
            fontsize (None): fontsize of text in the plot
            figsize ((20, 10)): figsize
            xlim (None): the 2theta range [x_min, x_max]
        """
        import matplotlib
        import matplotlib.pyplot as plt
        if fontsize is not None:
            matplotlib.rcParams.update({'font.size': fontsize})

        plt.figure(figsize=figsize)

        if xlim is None:
            x_min, x_max = 0, np.degrees(self.max2theta)
        else:
            x_min, x_max = xlim[0], xlim[1]
        dx = x_max-x_min
        for i in self.pxrd:
            plt.bar(i[0],i[-1], color='b', width=width*dx/180)
            if i[-1] > minimum_I and x_min <= i[0] <= x_max:
                if show_hkl:
                    label = self.draw_hkl(i[2:5])
                    plt.text(i[0]-dx/40, i[-1], label[0]+label[1]+label[2])

        #ax=plt.gca()
        plt.gca()
        plt.grid()
        plt.xlim([x_min, x_max])
        plt.xlabel('2$\Theta$ ($\lambda$=' + str(self.wavelength) + ' $\AA$)')
        plt.ylabel('Intensity')
        plt.title('PXRD of ' + self.name)

        if filename is None:
            plt.show()
        else:
            plt.savefig(filename)
            plt.close()


    def plotly_pxrd(self, profile='gaussian', minimum_I=0.01, res=0.02, FWHM=0.1, height=450, html=None):
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
            if I > minimum_I:
                label = '<br>2&#952;: {:6.2f}<br>d: {:6.4f}<br>'.format(theta2, d)
                label += 'I: {:6.4f}</br>hkl: ({:d}{:d}{:d})'.format(I, h, k, l)
                x.append(theta2)
                y.append(-0.1)
                labels.append(label)

        trace1 = go.Bar(x=x, y=y, text=labels,
                        hovertemplate = "%{text}",
                        width=0.5, name='hkl indices')
        if profile is None:
            fig = go.Figure(data=[trace1])
        else:
            spectra = self.get_profile(method=profile, res=res, user_kwargs={"FWHM": FWHM})
            trace2 = go.Scatter(x=spectra[0], y=spectra[1], name='Profile: ' + profile)
            fig = go.Figure(data=[trace2, trace1])

        fig.update_layout(height=height,
                          xaxis_title = '2&#952; ({:.4f} &#8491;)'.format(self.wavelength),
                          yaxis_title = 'Intensity',
                          title = 'PXRD of '+self.name)

        if os.environ['_'].find('jupyter') == -1:
            if html is None:
                return fig.to_html()
            else:
                fig.write_html(html)
        else:
            print("This is running on Jupyter Notebook")
            return fig


    def get_profile(self, method='gaussian', res=0.01, user_kwargs=None):
        """
        return the profile detail
        """

        return Profile(method, res, user_kwargs).get_profile(self.theta2, \
                self.xrd_intensity, np.degrees(self.min2theta), np.degrees(self.max2theta))


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

    def __init__(self, method='mod_pseudo-voigt', res=0.02, user_kwargs=None):

        self.method = method
        self.user_kwargs = user_kwargs
        self.res = res
        kwargs = {}

        if method == 'mod_pseudo-voigt':
            _kwargs = {
                      'U': 5.776410E-03,
                      'V': -1.673830E-03,
                      'W': 5.668770E-03,
                      'A': 1.03944,
                      'eta_h': 0.504656,
                      'eta_l': 0.611844,
                     }
        elif method in ['gaussian', 'lorentzian', 'pseudo-voigt']:
            _kwargs = {'FWHM': 0.1}

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

        N = int((max2theta-min2theta)/self.res)
        px = np.linspace(min2theta, max2theta, N)
        py = np.zeros((N))

        for two_theta, intensity in zip(two_thetas, intensities):
            if self.method == 'gaussian':
                fwhm = self.kwargs['FWHM']
                tmp = gaussian(two_theta, px, fwhm)

            elif self.method == 'lorentzian':
                fwhm = self.kwargs['FWHM']
                tmp = lorentzian(two_theta, px, fwhm)

            elif self.method == 'pseudo-voigt':
                try:
                    fwhm_g = self.kwargs['FWHM-G']
                    fwhm_l = self.kwargs['FWHM-L']
                except:
                    fwhm_g = self.kwargs['FWHM']
                    fwhm_l = self.kwargs['FWHM']

                fwhm = (fwhm_g**5 + 2.69269*fwhm_g**4*fwhm_l + 2.42843*fwhm_g**3*fwhm_l**2 +
                       4.47163*fwhm_g**2*fwhm_l**3 + 0.07842*fwhm_g*fwhm_l**4 + fwhm_l**5)**(1/5)
                eta = 1.36603*fwhm_l/fwhm - 0.47719*(fwhm_l/fwhm)**2 + 0.11116*(fwhm_l/fwhm)**3
                tmp = pseudo_voigt(two_theta, px, fwhm, eta)

            elif self.method == 'mod_pseudo-voigt':
                U = self.kwargs['U']
                V = self.kwargs['V']
                W = self.kwargs['W']
                A = self.kwargs['A']
                eta_h = self.kwargs['eta_h']
                eta_l = self.kwargs['eta_l']

                fwhm = np.sqrt(U*np.tan(np.pi*two_theta/2/180)**2 + V*np.tan(np.pi*two_theta/2/180) + W)
                x = px - two_theta
                tmp = mod_pseudo_voigt(x, fwhm, A, eta_h, eta_l, N)

            py += intensity * tmp

        py /= np.max(py)

        self.spectra = np.vstack((px,py))
        return self.spectra


# ------------------------------ Similarity between two XRDs ---------------------------------
class Similarity():

    def __init__(self, f, g, N = None, x_range = None, l = 2.0, weight = 'cosine'):

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
        res1 = (fx[-1] - fx[0])/len(fx)
        res2 = (gx[-1] - gx[0])/len(gx)
        self.resolution = min([res1, res2])/3 # improve the resolution

        if N is None:
            self.N = int(2*self.l/self.resolution)
        else:
            self.N = N
        self.r = np.linspace(-self.l, self.l, self.N)

        if x_range is None:
            x_min = max(np.min(fx), np.min(gx))
            x_max = min(np.max(fx), np.max(gx))
        else:
            x_min, x_max = x_range[0], x_range[1]

        self.x_range = [x_min,x_max]

        f_inter = interp1d(fx, fy, 'cubic', fill_value = 'extrapolate')
        g_inter = interp1d(gx, gy, 'cubic', fill_value = 'extrapolate')

        fgx_new = np.linspace(x_min, x_max, int((x_max-x_min)/self.resolution)+1)
        fy_new = f_inter(fgx_new)
        gy_new = g_inter(fgx_new)

        self.fx, self.gx, self.fy, self.gy = fgx_new, fgx_new, fy_new, gy_new
        self.weight = weight
        if self.weight == 'triangle':
            w = self.triangleFunction()
        elif self.weight == 'cosine':
            w = self.cosineFunction()
        else:
            msg = self.weight + 'is not supported'
            raise NotImplementedError(msg)

        Npts = len(self.fx)
        d = self.fx[1] - self.fx[0]

        self.value = similarity_calculate(self.r, w, d, Npts, self.fy, self.gy)

    def __str__(self):
        s = "The similarity between two PXRDs is {:.4f}".format(self.value)
        return s

    def __repr__(self):
        return str(self)

    def triangleFunction(self):

        """
        Triangle function to weight correlations
        """
        w = 1 - np.abs(self.r/self.l)
        ids = (np.abs(self.r) > self.l)
        w[ids] = 0

        return w

    def cosineFunction(self):

        """
        cosine function to weight correlations
        """

        w = 0.5 * (np.cos(np.pi * self.r/self.l) + 1.)
        ids = (np.abs(self.r) > self.l)
        w[ids] = 0

        return w

    def show(self, filename=None, fontsize=None, labels=["profile 1", "profile 2"]):

        """
        show the comparison plot

        Args:
            filename (None): name of the xrd plot. If None, show the plot
            labels [A, B]: labels of each plot
        """
        import matplotlib.pyplot as plt
        import matplotlib
        if fontsize is not None:
            matplotlib.rcParams.update({'font.size': fontsize})

        fig1 = plt.figure(1, figsize=(15, 6))
        fig1.add_axes((.1,.3,.8,.6))

        plt.plot(self.fx, self.fy, label=labels[0])
        plt.plot(self.fx, -self.gy, label=labels[1])
        plt.legend()

        # Residual plot
        residuals = self.gy - self.fy
        fig1.add_axes((.1,.1,.8,.2))
        plt.plot(self.gx, residuals, '.r', markersize = 0.5)
        plt.title("{:6f}".format(self.value))

        if filename is None:
            plt.show()
        else:
            plt.savefig(filename)
            plt.close()

@nb.njit(nb.f8[:](nb.f8[:], nb.f8, nb.f8, nb.f8, nb.f8, nb.i8), cache = True)
def mod_pseudo_voigt(x, fwhm, A, eta_h, eta_l, N):

    """
    A modified split-type pseudo-Voigt function for profiling peaks
    - Izumi, F., & Ikeda, T. (2000).
    """

    tmp = np.zeros((N))
    for xi, dx in enumerate(x):
        if dx < 0:
            A = A
            eta_l = eta_l
            eta_h = eta_h
        else:
            A = 1/A
            eta_l = eta_h
            eta_h = eta_l

        tmp[xi] = ((1+A)*(eta_h + np.sqrt(np.pi*np.log(2))*(1-eta_h))) /\
            (eta_l + np.sqrt(np.pi*np.log(2)) * (1-eta_l) + A*(eta_h +\
            np.sqrt(np.pi*np.log(2))*(1-eta_h))) * (eta_l*2/(np.pi*fwhm) *\
            (1+((1+A)/A)**2 * (dx/fwhm)**2)**(-1) + (1-eta_l)*np.sqrt(np.log(2)/np.pi) *\
            2/fwhm *np.exp(-np.log(2) * ((1+A)/A)**2 * (dx/fwhm)**2))
    return tmp

@nb.njit(nb.f8[:](nb.f8, nb.f8[:], nb.f8), cache = True)
def gaussian(theta2, alpha, fwhm):

    """
    Gaussian function for profiling peaks
    """

    tmp = ((alpha - theta2)/fwhm)**2
    return np.exp(-4*np.log(2)*tmp)

@nb.njit(nb.f8[:](nb.f8, nb.f8[:], nb.f8), cache = True)
def lorentzian(theta2, alpha, fwhm):

    """
    Lorentzian function for profiling peaks
    """

    tmp = 1 + 4*((alpha - theta2)/fwhm)**2
    return 1/tmp

def pseudo_voigt(theta2, alpha, fwhm, eta):

    """
    Original Pseudo-Voigt function for profiling peaks
    - Thompson, D. E. Cox & J. B. Hastings (1986).
    """

    L = lorentzian(theta2, alpha, fwhm)
    G = gaussian(theta2, alpha, fwhm)
    return eta * L + (1 - eta) * G

@nb.njit(nb.f8(nb.f8[:], nb.f8[:], nb.f8, nb.i8, nb.f8[:], nb.f8[:]))
def similarity_calculate(r, w, d, Npts, fy, gy):

    """
    Compute the similarity between the pair of spectra f, g
    """

    xCorrfg_w, aCorrff_w, aCorrgg_w = 0, 0, 0
    for r0, w0 in zip(r, w):
        Corrfg, Corrff, Corrgg = 0, 0, 0
        shift = int(r0/d)
        for i in range(Npts):
            if 0 <= i + shift <= Npts-1:
                Corrfg += fy[i]*gy[i+shift]*d
                Corrff += fy[i]*fy[i+shift]*d
                Corrgg += gy[i]*gy[i+shift]*d

        xCorrfg_w += w0*Corrfg*d
        aCorrff_w += w0*Corrff*d
        aCorrgg_w += w0*Corrgg*d

    return np.abs(xCorrfg_w / np.sqrt(aCorrff_w * aCorrgg_w))


def create_index():
    """
    shortcut to get the index
    """
    hkl_index = []
    for i in [-1,0,1]:
        for j in [-1,0,1]:
            for k in [-1,0,1]:
                hkl = np.array([i,j,k])
                if sum(hkl*hkl)>0:
                    hkl_index.append(hkl)
    hkl_index = np.array(hkl_index).reshape([len(hkl_index), 3])
    return hkl_index
