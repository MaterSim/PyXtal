#Standard Libraries
from math import acos, pi, ceil, sin,cos,sqrt
import numpy as np
import re
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
plt.style.use("ggplot")


#PyXtal imports
from pyxtal.database.element import Element

def angle(a,b):
    """ calculate the angle between vector a and b """
    return acos(np.dot(a,b)/np.linalg.norm(a)/np.linalg.norm(b))


class crystal(object):
    """a class of crystal structure. 
    Attributes:
        cell_para: a,b,c, alpha, beta, gamma
        cell_matrix: 3*3 matrix
        rec_matrix: reciprocal of cell matrix
        atom_type:  elemental type (e.g. Na Cl)
        composition: chemical composition (e.g., [1,1])
        coordinate: atomic positions (e.g., [[0,0,0],[0.5,0.5,0.5]])
    """

    def __init__(self, fileformat='POSCAR', filename=None, \
                 lattice=None, atom_type=None, composition=None, coordinate=None):
        """Return a structure object with the proper structures info"""
        if fileformat == 'POSCAR':
           self.from_POSCAR(filename)
        elif fileformat == 'cif':
           self.from_cif(filename)
        else:
           self.from_dict(lattice, atom_type, composition, coordinate)
    
    def from_cif(self, filename):
        cif_struc = cif(filename)
        lattice = self.para2matrix(cif_struc.cell_para)
        composition = cif_struc.composition
        coordinate = cif_struc.coordinate
        atom_type = cif_struc.atom_type
        self.from_dict(lattice, atom_type, composition, coordinate)

    def from_POSCAR(self, filename):

        f = open(filename)

        tag = f.readline()
        lattice_constant = float(f.readline().split()[0])

        # Now the lattice vectors
        a = []
        for ii in range(3):
            s = f.readline().split()
            floatvect = float(s[0]), float(s[1]), float(s[2])
            a.append(floatvect)
        lattice = np.array(a) * lattice_constant

        # Number of atoms. 
        atom_type = f.readline().split()
        comp = f.readline().split()
        composition = []
        if len(atom_type)==len(comp):
           for num in comp:
               composition.append(int(num))
        else:
           print('Value Error POSCAR symbol and composition is inconsistent')
        ac_type = f.readline().split()
        # Check if atom coordinates are cartesian or direct
        cartesian = ac_type[0].lower() == "c" or ac_type[0].lower() == "k"
        tot_natoms = sum(composition)
        coordinate = np.empty((tot_natoms, 3))
        for atom in range(tot_natoms):
            ac = f.readline().split()
            coordinate[atom] = (float(ac[0]), float(ac[1]), float(ac[2]))
        # Done with all reading
        f.close()
        if cartesian:
            coordinate *= lattice_constant
        cell_para = []
        self.from_dict(lattice, atom_type, composition, coordinate)

    def from_dict(self, lattice, atom_type, composition, coordinate):
        self.cell_matrix = np.array(lattice) 
        self.atom_type = atom_type
        self.composition = np.array(composition)
        self.coordinate = np.array(coordinate)
        self.cell_para = self.matrix2para(self.cell_matrix)
        self.rec_matrix = self.rec_lat(self.cell_matrix)
        self.name = ''
        for ele, num in zip(self.atom_type, self.composition):
            self.name += ele
            if num > 1:
               self.name += str(num)
   
    @staticmethod
    def rec_lat(matrix):
        """ calculate the reciprocal lattice """
        rec_lat = np.zeros([3,3])
        V = np.linalg.det(matrix)
        rec_lat[0] = np.cross(matrix[1], matrix[2])/V
        rec_lat[1] = np.cross(matrix[2], matrix[0])/V
        rec_lat[2] = np.cross(matrix[0], matrix[1])/V
        return  rec_lat #* 2 * pi

    @staticmethod
    def matrix2para(matrix):
        """ 3x3 representation -> 1x6 (a, b, c, alpha, beta, gamma)"""
        cell_para = np.zeros(6)
        cell_para[0] = np.linalg.norm(matrix[0])
        cell_para[1] = np.linalg.norm(matrix[1])
        cell_para[2] = np.linalg.norm(matrix[2])
    
        cell_para[5] = angle(matrix[0], matrix[1])
        cell_para[4] = angle(matrix[0], matrix[2])
        cell_para[3] = angle(matrix[1], matrix[2])

        return cell_para

    @staticmethod
    def para2matrix(cell_para):
        """ 1x6 (a, b, c, alpha, beta, gamma) -> 3x3 representation -> """
        matrix = np.zeros([3,3])
        matrix[0][0] = cell_para[0]
        matrix[1][0] = cell_para[1]*cos(cell_para[5])
        matrix[1][1] = cell_para[1]*sin(cell_para[5])
        matrix[2][0] = cell_para[2]*cos(cell_para[4])
        matrix[2][1] = cell_para[2]*cos(cell_para[3])*sin(cell_para[4])
        matrix[2][2] = sqrt(cell_para[2]**2 - matrix[2][0]**2 - matrix[2][1]**2)
        
        return matrix
    
class cif(object):
    """a class of cif reader
    Attributes:
        wavelength: default: 1.54181a, namely Cu-Ka
        max2theta: the range of 2theta angle
        intensity: intensities for all hkl planes
        pxrd: powder diffraction data
    """

    def __init__(self, filename):
        """Return a XRD object with the proper info"""
        self.from_file(filename)
        self.parse_cell()
        self.parse_atom()
        self.apply_symops()

    def from_file(self, filename):
        cif = np.genfromtxt(filename, dtype=str, delimiter='\n')
        
        # 3 modes in each flag:  
        # 0: not started; 
        # 1: reading; 
        # 2: done
        flags = {'cell':0, 'symops':0, 'atom':0}

        atom = {}
        cell = {}
        symops = {'string':[], 'matrix':[]}

        for lines in cif:

            if 'loop_' in lines:  
                #if a _loop lines starts, the current reading flag switch to 0
                for item in flags.keys():
                    if flags[item] == 1:
                        flags[item] = 2

            elif '_cell_length_' in lines or '_cell_angle_' in lines:
                #_cell_length_a          4.77985

                flags['cell'] = 1
                cell_str = lines.split()
                item = cell_str[0].replace(' ','')
                value = float(cell_str[1].split("(")[0])
                cell[item] = value

            elif '_symmetry_equiv_pos_as_xyz' in lines:
                #_symmetry_equiv_pos_as_xyz
                flags['symops'] = 1
      
            elif '_space_group_symop_operation_xyz' in lines:
                #_space_group_symop_operation_xyz
                flags['symops'] = 1
                
            elif flags['symops'] == 1:
                #1, 'x, y, z'
                #    x, -y, z
                raw_line = lines.strip().strip("'").split(' ', 1)
                if raw_line[0].isdigit():     
                    sym_str = raw_line[1].strip("'")
                else:
                    sym_str = lines.strip().strip("'").replace(' ', '')
                sym_str = sym_str.replace("'","")
                symops['string'].append(sym_str)
                symops['matrix'].append(self.xyz2sym_ops(sym_str))

            elif '_atom_site' in lines: 
                flags['atom'] = 1
                atom_str = lines.replace(' ','')
                item = atom_str
                atom[item] = []

            elif flags['atom'] == 1:
                raw_line = lines.split()
                for i, item in enumerate(atom.keys()):
                    raw_text = raw_line[i]
                    
                    if item.find('fract')>0:
                       value = float(raw_text.split("(")[0])
                    elif item.find('symbol')>0:
                       m_symbol = re.compile("([A-Z]+[a-z]*)")
                       value = str(m_symbol.findall(raw_text)).strip("[]").strip("''")
                       #print(raw_text, value)
                    else:
                       value = raw_text
                       
                    atom[item].append(value)

            elif flags['cell'] + flags['symops'] + flags['atom'] == 6:
                break

        self.cell = cell
        self.atom = atom
        self.symops = symops
   
    def parse_cell(self):
        cell_para = np.zeros(6)
        cell = self.cell
        for item in cell.keys():
            if item.find('_length_a') > 0:
                cell_para[0] = cell[item]
            elif item.find('_length_b') > 0:
                cell_para[1] = cell[item]
            elif item.find('_length_c') > 0:
                cell_para[2] = cell[item]
            elif item.find('_angle_alpha') > 0:
                cell_para[3] = np.radians(cell[item])
            elif item.find('_angle_beta') > 0:
                cell_para[4] = np.radians(cell[item])
            elif item.find('_angle_gamma') > 0:
                cell_para[5] = np.radians(cell[item])
        self.cell_para = cell_para

    def parse_atom(self):
        atom = self.atom
        N_atom = len(atom['_atom_site_fract_x'])
        cif_xyz = np.zeros([N_atom, 3])

        for item in atom.keys():
            if item.find('_fract_x') > 0:
                cif_xyz[:,0] = np.array(atom[item])
            elif item.find('_fract_y') > 0:
                cif_xyz[:,1] = np.array(atom[item])
            elif item.find('_fract_z') > 0:
                cif_xyz[:,2] = np.array(atom[item])

        self.cif_xyz = cif_xyz

    #generates all coordinates from rotation matrices and translation vectors
    def apply_symops(self):
        fract_xyz = self.cif_xyz
        symops_matrix = self.symops['matrix']
        atom_type = self.atom['_atom_site_type_symbol']
        sym_coordinates = {}
        
        for item in atom_type:
            sym_coordinates[item] = []


        for ii,item in enumerate(atom_type):
            for mat_vec in symops_matrix:
                sym_temp = np.dot(mat_vec[0], fract_xyz[ii].transpose()) + mat_vec[1]
                sym_coordinates[item].append(sym_temp)
        self.coordinate, self.composition, self.atom_type = \
                      self.remove_duplicate(sym_coordinates)

    #remove equivalent points and keep the unique ones
    #get the numbers of atoms per species
    @staticmethod
    def remove_duplicate(sym_coordinates):
        coordinate = []
        composition = []
        atom_type = []
        for item in sym_coordinates.keys():
            atom_type.append(item)
            raw_equiv = np.array(sym_coordinates[item])
            raw_equiv = raw_equiv - np.floor(raw_equiv)
            raw_equiv = np.around(raw_equiv, 4)
            raw_equiv = np.unique(raw_equiv, axis=0)
            composition.append(len(raw_equiv))
            if coordinate == []:
                coordinate = raw_equiv
            else:
                coordinate = np.concatenate((coordinate,raw_equiv),axis=0)

        return coordinate, composition, atom_type


    #function generates rotation matrices and translation vectors from equivalent points
    @staticmethod
    def xyz2sym_ops(string):
        #rotational matrix dictionary
        rot_dic = {}
        rot_dic['x'] = np.array([1.0,0,0])
        rot_dic['y'] = np.array([0,1.0,0])
        rot_dic['z'] = np.array([0,0,1.0])
        parts = string.strip().replace(' ','').lower().split(',')
        rot_mat = []
        rot_temp = np.array([0.,0.,0.])
        trans_vec = np.array([0.,0.,0.])
        #use re module to read xyz strings
        m_rot = re.compile(r"([+-]?)([\d\.]*)/?([\d\.]*)([x-z])")
        m_trans = re.compile(r"([+-]?)([\d\.]+)/?([\d\.]*)(?![x-z])")
        for jj,item in enumerate(parts):
            #rotation matrix
            for ii,m in enumerate(m_rot.finditer(item)):
                coef = -1 if m.group(1) == '-' else 1
                if m.group(2) != '':
                    if m.group(3) != '':
                        coef *= float(m.group(2))/float(m.group(3))
                    else:
                        coef *= float(m.group(2))
                if ii == 0:                  
                    rot_temp = rot_dic[m.group(4)]*coef
                else:
                    rot_temp += rot_dic[m.group(4)]*coef
            rot_mat.append(rot_temp)
            #translation vector
            for m in m_trans.finditer(item):
                coef = -1 if m.group(1) == '-' else 1
                if m.group(3) != '':
                    coef = float(m.group(2))/float(m.group(3))
                else:
                    coef = float(m.group(2))
                trans_vec[jj] = 1.0*coef
        return (rot_mat, trans_vec)
         

class XRD(object):
    """a class of crystal structure. 
    Attributes:
        cell_para: a,b,c, alpha, beta, gamma
        cell_matrix: 3*3 matrix
        rec_matrix: reciprocal of cell matrix
        atom_type:  elemental type (e.g. Na Cl)
        composition: chemical composition (e.g., [1,1])
        coordinate: atomic positions (e.g., [[0,0,0],[0.5,0.5,0.5]])
    """

    def __init__(self, crystal, wavelength=1.54184, max2theta=180):
        """Return a XRD object with the proper info"""
        self.wavelength = wavelength
        self.max2theta = np.radians(max2theta)
        self.name = crystal.name
        self.all_dhkl(crystal)
        self.atom_scatter(crystal)
        self.structure_factor(crystal)
        self.rec_matrix = crystal.rec_matrix
        self.intensity()
        self.pxrd()

    def by_hkl(self, hkl):
        """ d for any give abitray [h,k,l] index """
        id1 = np.where(np.all(self.hkl_list == np.array(hkl), axis=1 ))
        if id1 is None:
           print('This hkl is not in the given 2theta range')
        else:
           print('  2theta     d_hkl     hkl       Intensity')
           for i in id1[0]:
              #print(len(i), self.xrd_intensity[i])
              print('%8.3f  %8.3f   [%2d %2d %2d] %8.2f' % \
                    (np.degrees(self.theta2[i]), self.d_hkl[i], \
                     self.hkl_list[i,0], self.hkl_list[i,1], self.hkl_list[i,2], \
                     self.xrd_intensity[i] ))        
           #return np.degrees(self.theta2[id1]), self.d_hkl[id1], self.xrd_intensity[id1] 

    def all_dhkl(self, crystal):
        """ 3x3 representation -> 1x6 (a, b, c, alpha, beta, gamma)"""
        #d_min = self.wavelength/self.max2theta*pi/2
        d_min = self.wavelength/sin(self.max2theta/2)/2
  
        # This block is to find the shortest d_hkl, 
        # for all basic directions (1,0,0), (0,1,0), (1,1,0), (1,-1,0) and so on, 26 in total 
        hkl_max = np.array([1,1,1])
        for h1 in [-1, 0, 1]:
            for k1 in [-1, 0, 1]:
                for l1 in [-1, 0, 1]:
                    hkl_index = np.array([[h1,k1,l1]])
                    d = float(np.linalg.norm( np.dot(hkl_index, crystal.rec_matrix), axis=1))
                    if d>0:
                       multiple = 1/d/d_min
                       hkl_index *= round(multiple)
                       for i in range(len(hkl_max)):
                           if hkl_max[i] < hkl_index[0,i]:
                              hkl_max[i] = hkl_index[0,i]
        #h1 = 2*ceil(np.linalg.norm(crystal.cell_para[0])/d_min)
        #k1 = 2*ceil(np.linalg.norm(crystal.cell_para[1])/d_min)
        #l1 = 2*ceil(np.linalg.norm(crystal.cell_para[2])/d_min)
        h1, k1, l1 = hkl_max
        h = np.arange(-h1,h1+1)
        k = np.arange(-k1,k1+1)
        l = np.arange(-l1,l1+1)
        
        hkl = np.array((np.meshgrid(h,k,l))).transpose()
        hkl_list = np.reshape(hkl, [len(h)*len(k)*len(l),3])
        hkl_list = hkl_list[np.where(hkl_list.any(axis=1))[0]]
        d_hkl = 1/np.linalg.norm( np.dot(hkl_list, crystal.rec_matrix), axis=1)
        #for ix, a in enumerate(hkl_list):
        #    if np.array_equal(a, np.array([1,-1,3])) is True:
        #       print(a)
        #       break
        #
        #print(ix, hkl_list[ix], d_hkl[ix], d_min)

        shortlist = d_hkl > (d_min)
        d_hkl = d_hkl[shortlist]
        hkl_list = hkl_list[shortlist]
        sintheta = self.wavelength/2/d_hkl

        self.theta = np.arcsin(sintheta)
        self.hkl_list = hkl_list
        self.d_hkl = d_hkl
        
        #return hkl_list, d_hkl, sintheta

    def atom_scatter(self, crystal):
        """ N*M array; N: atoms, M: N_hkl"""
        f = np.zeros([sum(crystal.composition), len(self.d_hkl)])
        d0 = 1/2/self.d_hkl
        count = 0
        for i, ele in enumerate(crystal.atom_type):
            c = Element(ele).scatter
            f_tmp = c[0]*np.exp(-c[4]*d0) + \
                    c[1]*np.exp(-c[5]*d0) + \
                    c[2]*np.exp(-c[6]*d0) + \
                    c[3]*np.exp(-c[7]*d0) + c[8]
            for j in range(count,count+crystal.composition[i]):
                f[j] = f_tmp
            count += crystal.composition[i]

        self.f = f
   
    def structure_factor(self, crystal):
        """ N*1 array"""
        F = []
        for fj, hkl in zip(self.f.transpose(), self.hkl_list):
            F_tmp = np.exp(2*pi*1j*np.dot(crystal.coordinate, hkl.transpose()))
            F.append(np.dot(fj, F_tmp))

        self.F = np.array(F)

    def intensity(self):
        """" Calculate intensity, return N*1 array """
        LP = 1/np.sin(self.theta)**2/np.cos(self.theta)
        P = 1 + np.cos(2*self.theta)**2
        I = (np.abs(self.F))**2*LP*P
        self.xrd_intensity = I
        self.theta2 = 2*self.theta
        rank = np.argsort(self.theta2)
        self.theta2 = self.theta2[rank]
        self.hkl_list = self.hkl_list[rank]
        self.d_hkl = self.d_hkl[rank]
        self.xrd_intensity = self.xrd_intensity[rank]

    def pxrd(self):
        """ Group the equivalent hkl planes together by 2\theta angle
            N*6 arrays, Angle, d_hkl, h, k, l, intensity
        """
        rank = range(len(self.theta2)) #np.argsort(self.theta2)
        PL = []
        last = []
        for i in rank:
            if self.xrd_intensity[i] > 0.01:
               angle = np.degrees(self.theta2[i])
               if PL is None:
                  PL.append([angle, self.d_hkl[i], \
                             self.hkl_list[i,0], self.hkl_list[i,1], self.hkl_list[i,2], \
                             self.xrd_intensity[i]])
               elif abs(angle-last) < 1e-4:
                  PL[-1][-1] += self.xrd_intensity[i]
               else:
                  PL.append([angle, self.d_hkl[i], \
                             self.hkl_list[i,0], self.hkl_list[i,1], self.hkl_list[i,2], \
                             self.xrd_intensity[i]])
               last = angle
        PL = (np.array(PL))
        PL[:,-1] = PL[:,-1]/max(PL[:,-1])
        self.pxrd = PL

    def plot_pxrd(self, filename=None, minimum_I = 0.01, show_hkl=True):
        """ plot PXRD """

        #print('  2theta     d_hkl     hkl       Intensity')
        dx = np.degrees(self.max2theta)
        for i in self.pxrd:
            plt.bar(i[0],i[-1], color='b', width=dx/180)
            if i[-1] > minimum_I:
               if show_hkl:
                  label = self.draw_hkl(i[2:5])
                  plt.text(i[0]-dx/40, i[-1], label[0]+label[1]+label[2])
   
        ax=plt.gca()
        plt.grid()
        plt.xlim(0,dx)
        plt.xlabel('2θ')
        plt.ylabel('Intensity')
        plt.title('PXRD of '+self.name+ ', $\lambda$='+str(self.wavelength)+'$\AA$')
        if filename is None:
           plt.show()
        else:
           plt.savefig(filename)
           plt.close()

    @staticmethod
    def draw_hkl(hkl):
        """turn negative numbers in hkl to overbar"""
        hkl_str= []
        for i in hkl:
            if i<0:
               label = str(int(-i))
               label = r"$\bar{" + label + '}$'
               hkl_str.append(str(label))
            else:
               hkl_str.append(str(int(i)))

        return hkl_str

from optparse import OptionParser
import pandas as pd
from tabulate import tabulate

if __name__ == "__main__":
    #-------------------------------- Options -------------------------
    parser = OptionParser()
    parser.add_option("-m", "--hkl", dest="hkl", metavar='hkl index',
                      help="show hkl_index info, e.g., [1,0,0]")
    parser.add_option("-a", "--angle", dest="max2theta", default=180, type='float',
                      help="2theta angle range, default=180", metavar="angle")
    parser.add_option("-t", "--transform", dest="trans", metavar="files",
                      help="export file in different format")
    parser.add_option("-p", "--plot", dest="plot", 
                      help="plot pxrd", metavar="plot")
    parser.add_option("-w", "--wavelength", dest="wavelength", default=1.54184, type='float',
                      help="wavelength: 1.54184", metavar="wave")
    parser.add_option("-c", "--crystal", dest="structure",default='',
                      help="crystal from file, cif or poscar, REQUIRED", metavar="crystal")
    parser.add_option("-f", "--full", dest="full", action='store_true',
                      help="show full hkl reflections", metavar="full")
    parser.add_option("-i", "--intensity", dest="minimum_I", default=0.01, type='float',
                      help="the minimum intensity to show, default 0.01", metavar="intensity")
    
    
    
    (options, args) = parser.parse_args()    
    if options.structure.find('cif') > 0:
        fileformat = 'cif'
    else:
        fileformat = 'POSCAR'
    
    test = crystal(fileformat, filename=options.structure)
    xrd = XRD(test, wavelength=options.wavelength, \
                   max2theta=options.max2theta)   
    if options.full: 
        col_name = {'2theta': xrd.pxrd[:,0], \
                  'd_hkl':  xrd.pxrd[:,1], \
                  'h': xrd.pxrd[:,2], \
                  'k': xrd.pxrd[:,3], \
                  'l': xrd.pxrd[:,4], \
                  'Intensity':xrd.pxrd[:,5]}
    else:
        rank1 = xrd.xrd_intensity > options.minimum_I
        col_name = {'2theta':    np.degrees(xrd.theta2[rank1]), \
                  'd_hkl':     xrd.d_hkl[rank1],\
                  'h':         xrd.hkl_list[rank1,0], \
                  'k':         xrd.hkl_list[rank1,1], \
                  'l':         xrd.hkl_list[rank1,2], \
                  'Intensity': xrd.xrd_intensity[rank1] }
    
    df = pd.DataFrame(col_name)
    print(tabulate(df, headers='keys')) #, tablefmt='psql'))
    if options.plot is not None:
        xrd.plot_pxrd(filename=options.plot, minimum_I = options.minimum_I)
