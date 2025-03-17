import contextlib
import os
import shutil
import subprocess
import numpy as np

class CHARMM:
    """
    A calculator to perform oragnic crystal structure optimization in CHARMM.

    Args:
        - struc: pyxtal.molecular_crystal.molecular_crystal object
        - label (str): label for this calculation
        - algo (str): "abnr"
        - lat_mut (bool): mutate lattice or not
        - rotate (bool): rotate lattice or not
        - prefix (str): prefix of this calculation
        - atom_info (dict): atom_labels
        - folder (str): folder name
        - opt (str): 'conv', 'conp', 'single'
        - steps (int): optimization steps
        - exe (str): charmm executable
        - input (str): charmm input file
        - output (str): charmm output file
        - dump (str): charmm dump structure
    """

    def __init__(
        self,
        struc,
        label="_",
        algo="abnr",
        lat_mut=False,
        rotate=False,
        prefix="pyxtal",
        atom_info=None,
        folder=".",
        opt="conp",
        steps=None,
        exe="charmm",
        input="charmm.in",
        output="charmm.log",
        dump="result.pdb",
        debug=False,
        timeout=20,
    ):
        self.errorE = 1e+5
        self.error = False
        if steps is None:
            steps = [2000, 1000]
        self.debug = debug
        self.timeout = timeout
        # check charmm Executable
        if shutil.which(exe) is None:
            raise BaseException(f"{exe} is not installed")
        else:
            self.exe = exe

        self.atom_info = atom_info

        # Optimization setting
        self.opt = opt
        self.optlat = False
        self.algo = algo
        if type(steps) is list:
            self.steps = steps
        elif type(steps) is int:
            self.steps = [2000, steps]
        self.opt = opt

        # Files IO
        self.folder = folder
        self.prefix = prefix
        self.label = label
        self.rtf = self.prefix + ".rtf"
        self.prm = self.prefix + ".prm"
        self.psf = self.label + self.prefix + ".psf"
        self.crd = self.label + self.prefix + ".crd"
        self.input = self.label + input
        self.output = self.label + output
        self.dump = self.label.lower() + dump  # charmm only output lower

        # For parsing the output
        self.positions = None
        self.forces = None
        self.optimized = False
        self.cell = None
        self.cputime = 0.0

        # Structure Manipulation
        struc.resort()
        self.structure = struc
        try:
            self.structure.optimize_lattice()
        except:
            print("bug in Lattice")
            print(self.structure)
            print(self.structure.lattice)
            print(self.structure.lattice.matrix)
            #raise ValueError("Problem in Lattice")
            self.error = True
        # print("\nbeginining lattice: ", struc.lattice)

        self.lat_mut = lat_mut
        self.rotate = rotate

    def run(self, clean=True):
        """
        Only run calc if it makes sense
        """
        if not self.error:
            os.makedirs(self.folder, exist_ok=True)
            cwd = os.getcwd()
            os.chdir(self.folder)

            self.write()  # ; print("write", time()-t0)
            res = self.execute()  # ; print("exe", time()-t0)
            if res is not None:
                self.read()  # ; print("read", self.structure.energy)
            else:
                self.structure.energy = self.errorE
                self.error = True
            if clean:
                self.clean()

            os.chdir(cwd)

    def execute(self):
        cmd = self.exe + " < " + self.input + " > " + self.output
        # os.system(cmd)
        with open(os.devnull, 'w') as devnull:
            try:
                # Run the external command with a timeout
                result = subprocess.run(
                    cmd, shell=True, timeout=self.timeout, check=True, stderr=devnull)
                return result.returncode  # Or handle the result as needed
            except subprocess.CalledProcessError as e:
                print(f"Command '{cmd}' failed with return code {e.returncode}.")
                return None
            except subprocess.TimeoutExpired:
                print(f"External command {cmd} timed out.")
                return None

    def clean(self):
        os.remove(self.input)
        os.remove(self.output)
        os.remove(self.crd) if os.path.exists(self.crd) else None
        os.remove(self.psf) if os.path.exists(self.psf) else None
        os.remove(self.dump) if os.path.exists(self.dump) else None

    def write(self):
        """
        setup the necessary files for charmm calculation
        """
        lat = self.structure.lattice
        if self.lat_mut:
            lat = lat.mutate()

        a, b, c, alpha, beta, gamma = lat.get_para(degree=True)
        ltype = lat.ltype
        if ltype in ['trigonal', 'Trigonal']:
            ltype = 'hexagonal'

        fft = self.FFTGrid(np.array([a, b, c]))

        with open(self.input, "w") as f:
            # General input
            f.write("! Automated Charmm calculation\n\n")
            f.write("bomlev -1\n")
            f.write("! top and par\n")
            f.write(f"read rtf card name {self.rtf:s}\n")
            f.write(f"read para card name {self.prm:s}\n")
            f.write("Read sequence card\n")
            f.write(f"{len(self.structure.mol_sites):5d}\n")
            atom_count = []
            for site in self.structure.mol_sites:
                atom_count.append(len(site.molecule.mol))
                if self.atom_info is None:
                    f.write(f"U0{site.type:d} ")
                else:
                    f.write("{:s} ".format(
                        self.atom_info["resName"][site.type]))

            f.write("\ngenerate main first none last none setup warn\n")
            f.write("Read coor card free\n")
            f.write("* Residues coordinate\n*\n")
            f.write(f"{sum(atom_count):5d}\n")
            for i, site in enumerate(self.structure.mol_sites):
                res_name = f"U0{site.type:d}" if self.atom_info is None else self.atom_info[
                    "resName"][site.type]

                # reset lattice if needed (to move out later)
                site.lattice = lat

                # Rotate if needed
                if self.rotate:
                    site.perturbate(lat.matrix, trans=0.5, rot="random")
                coords, species = site._get_coords_and_species(first=True)
                if i == 0:
                    count = 0
                else:
                    count += atom_count[i - 1]

                for j, coord in enumerate(coords):
                    if self.atom_info is None:
                        label = f"{species[j]:s}{j + 1:d}"
                    else:
                        label = self.atom_info["label"][site.type][j]
                    f.write(
                        "{:5d}{:5d}{:>4s}  {:<4s}{:10.5f}{:10.5f}{:10.5f}\n".format(
                            j + 1 + count, i + 1, res_name, label, *coord
                        )
                    )
                    # quickly check if
                    if abs(coord).max() > 500.0:
                        print("Unexpectedly large input coordinates, stop and debug")
                        print(self.structure)
                        self.structure.to_file('bug.cif')
                        import sys
                        sys.exit()

            f.write(f"write psf card name {self.psf:s}\n")
            f.write(f"write coor crd card name {self.crd:s}\n")
            f.write(f"read psf card name {self.psf:s}\n")
            f.write(f"read coor card name {self.crd:s}\n")

            # Structure info
            f.write("\n! crystal parameters\n")
            f.write(f"set shape {ltype:s}\n")
            f.write(f"set a     {a:12.6f}\n")
            f.write(f"set b     {b:12.6f}\n")
            f.write(f"set c     {c:12.6f}\n")
            f.write(f"set alpha {alpha:12.6f}\n")
            f.write(f"set beta  {beta:12.6f}\n")
            f.write(f"set gamma {gamma:12.6f}\n")
            f.write("coor conv FRAC SYMM @a @b @c @alpha @beta @gamma\n")
            f.write("coor stat select all end\n")
            f.write("Crystal Define @shape @a @b @c @alpha @beta @gamma\n")
            site0 = self.structure.mol_sites[0]
            f.write(
                f"Crystal Build cutoff 14.0 noperations {len(site0.wp.ops) - 1:d}\n")
            for i, op in enumerate(site0.wp.ops):
                if i > 0:
                    f.write(f"({op.as_xyz_str():s})\n")

            f.write(
                "image byres xcen ?xave ycen ?yave zcen ?zave sele resn LIG end\n")
            f.write("set 7 fswitch\n")
            f.write("set 8 atom\n")
            f.write("set 9 vatom\n")
            f.write("Update inbfrq 10 imgfrq 10 ihbfrq 10 -\n")
            f.write(
                "ewald pmewald lrc fftx {:d} ffty {:d} fftz {:d} -\n".format(*fft))
            f.write("kappa 0.34 order 6 CTOFNB 12.0 CUTNB 14.0 QCOR 1.0 -\n")
            f.write("@7 @8 @9 vfswitch !\n")
            f.write(f"mini {self.algo:s} nstep {self.steps[0]:d}\n")
            if len(self.steps) > 1:
                f.write(
                    f"mini {self.algo:s} lattice nstep {self.steps[1]:d} \n")
            if len(self.steps) > 2:
                f.write(f"mini {self.algo:s} nstep {self.steps[2]:d}\n")

            f.write(
                "coor conv SYMM FRAC ?xtla ?xtlb ?xtlc ?xtlalpha ?xtlbeta ?xtlgamma\n")  #
            f.write(f"\nwrite coor pdb name {self.dump:s}\n")  #
            f.write("*CELL :  ?xtla  ?xtlb  ?xtlc ?xtlalpha ?xtlbeta ?xtlgamma\n")  #
            f.write(f"*Z = {len(site0.wp):d}\n")
            f.write("*Energy(kcal): ?ener\n")
            f.write("stop\n")
        # print("STOP")
        # import sys
        # sys.exit()

    def read(self):
        with open(self.output) as f:
            lines = f.readlines()
            self.version = lines[2]
            if lines[-1].find("CPU TIME") != -1:
                self.optimized = True
                self.cputime = float(lines[-2].split()[-2])
                for line in lines:
                    if line.find("MINI> ") != -1:
                        with contextlib.suppress(Exception):
                            self.structure.iter = int(line.split()[1])
                            # raise RuntimeError("Something is wrong in the charmm output")
                    elif line.find("ABNORMAL TERMINATION") != -1:
                        self.optimized = False
                        break

        if self.optimized:
            with open(self.dump) as f:
                lines = f.readlines()
                positions = []
                for line in lines:
                    if line.find("REMARK ENERGY") != -1:
                        tmp = line.split(":")
                        self.structure.energy = float(tmp[-1])
                    elif line.find("REMARK Z") != -1:
                        tmp = line.split(":")[-1].split()
                        Z = int(tmp[-1])
                    elif line.find("REMARK CELL") != -1:
                        tmp = line.split(":")[-1].split()
                        tmp = [float(x) for x in tmp]
                        self.structure.lattice.set_para(tmp)
                    elif line.find("ATOM") != -1:
                        try:
                            xyz = line.split()[5:8]
                            XYZ = [float(x) for x in xyz]
                            positions.append(XYZ)
                        except:
                            # print("Warning: BAD charmm output: " + line)
                            pass
                positions = np.array(positions)
                self.structure.energy *= Z

            count = 0
            # if True:
            try:
                for _i, site in enumerate(self.structure.mol_sites):
                    coords = positions[count: count + len(site.molecule.mol)]
                    site.update(coords, self.structure.lattice)
                    count += len(site.molecule.mol)
                # print("after relaxation  : ", self.structure.lattice, "iter: ", self.structure.iter)
                self.structure.optimize_lattice()
                self.structure.update_wyckoffs()
                # print("after latticeopt  : ", self.structure.lattice, self.structure.check_distance()); import sys; sys.exit()
            except:
                # molecular connectivity or lattice optimization
                self.structure.energy = self.errorE
                self.error = True
                if self.debug:
                    print("Cannot retrieve Structure after optimization")
                    print("lattice", self.structure.lattice)
                    self.structure.to_file("1.cif")
                    #print("Check 1.cif in ", os.getcwd())
                    pairs = self.structure.check_short_distances()
                    #if len(pairs) > 0:
                    #    print(self.structure.to_file())
                    #    print("short distance pair", pairs)

        else:
            self.structure.energy = self.errorE
            self.error = True
            if self.debug:
                print(self.structure)
                import sys; sys.exit()

    def FFTGrid(self, ABC):
        """
        determine the grid used for
        """

        grid = np.array(
            [
                1,
                2,
                3,
                4,
                5,
                6,
                8,
                9,
                10,
                12,
                15,
                16,
                18,
                20,
                24,
                25,
                27,
                30,
                32,
                36,
                40,
                45,
                48,
                50,
                54,
                60,
                72,
                75,
                80,
                81,
                90,
                96,
                100,
                108,
                120,
                125,
                135,
                144,
                150,
                160,
                162,
                180,
                200,
                216,
                225,
                240,
                243,
                250,
                270,
                288,
                300,
                324,
                360,
                375,
                400,
            ]
        )
        fftxyz = [0, 0, 0]
        for i, l in enumerate(ABC):
            tmp = grid[grid > l]
            fftxyz[i] = max([32, tmp[0]])  # minimum is 32
        return fftxyz


def check_prm(path):
    """
    Todo: move it to charmm interface
    Sometimes there are something wrong with prm file (RUBGIH),
    """
    with open(path) as f:
        lines = f.readlines()
        pairs = []
        triplets = []
        imphis = []
        ids = []
        do_angle = False
        do_bond = False
        do_imphi = False
        for i, l in enumerate(lines):
            tmp = l.split()
            if l.find("ANGLE") == 0:
                do_angle = True
            elif l.find("BOND") == 0:
                do_bond = True
            elif l.find("DIHEDRAL") == 0:
                pass
            elif l.find("IMPHI") == 0:
                do_imphi = True
            elif do_bond:
                if len(tmp) > 0:
                    pair1 = tmp[0] + " " + tmp[1]
                    pair2 = tmp[1] + " " + tmp[0]
                    if (pair1 in pairs) or (pair2 in pairs):
                        ids.append(i)
                        print("Duplicate bonds: ", l[:-2])
                    else:
                        pairs.append(pair1)
                else:
                    do_bond = False
            elif do_angle:
                if len(tmp) > 0:
                    pair1 = tmp[0] + " " + tmp[1] + " " + tmp[2]
                    pair2 = tmp[2] + " " + tmp[1] + " " + tmp[0]
                    if (pair1 in triplets) or (pair2 in triplets):
                        ids.append(i)
                        print("Duplicate angles: ", l[:-2])
                    else:
                        triplets.append(pair1)
                else:
                    do_angle = False
            # elif do_dihedral:
            #    if len(tmp) > 0:
            #        pair1 = tmp[0]+ ' ' + tmp[1] + ' ' + tmp[2] + ' ' +tmp[3] + ' ' + tmp[5]
            #        pair2 = tmp[3]+ ' ' + tmp[1] + ' ' + tmp[2] + ' ' +tmp[0] + ' ' + tmp[5]
            #        if (pair1 in quads) or (pair2 in quads):
            #            ids.append(i)
            #            #print(pair1)
            #            #print(pair2)
            #            #print(quads)
            #            print('Duplicate dihedral angles: ', l[:-2])
            #            #import sys
            #            #sys.exit()
            #        else:
            #            quads.append(pair1)
            #    else:
            #        do_dihedral = False
            elif do_imphi:
                if len(tmp) > 0:
                    pair1 = tmp[0] + " " + tmp[1] + " " + tmp[2] + " " + tmp[3]
                    # pair2 = tmp[0] + ' ' + tmp[1] + ' ' + tmp[2] + ' ' +tmp[0]
                    if pair1 in imphis:  # or (pair2 in imphis):
                        ids.append(i)
                        print("Duplicate imphi angles: ", l[:-2])
                    else:
                        imphis.append(pair1)
                else:
                    do_imphi = False
                    break

    lines = [lines[i] for i in range(len(lines)) if i not in ids]
    with open(path, "w") as f:
        f.writelines(lines)


class RTF:
    """
    Parse, convert, merge the RTF files for CHARMM
    """

    def __init__(self, input):
        self.keywords = [
            "MASS",
            "RESI",
            "CHARGE",
            "ATOM",
            "BOND",
            "ANGL",
            "DIHE",
            "IMPH",
        ]
        if type(input) == str:
            print("Converting RTF from the file", input)
            f = open(input)
            rtf = f.read().split("\n")
            f.close()
        else:
            rtf = input

        self.parse(rtf)

    def parse(self, rtf):
        """
        parse rtf list to get the dictionary
        """
        mass = []
        self.residues = []

        do_resi = False
        residue = None
        for l in rtf:
            tmp = l.split()
            if l.find("MASS") == 0:
                mass.append(l)
                do_resi = False
            elif l.find("RESI") == 0:
                if residue is not None:
                    self.residues.append(residue)
                do_resi = True
                residue = {
                    "NAME": tmp[1],
                    "CHARGE": float(tmp[2]),
                    "ATOM": [],
                    "BOND": [],
                    "ANGL": [],
                    "DIHE": [],
                    "IMPH": [],
                }

            if do_resi:
                if l.find("ATOM") == 0:  # len(tmp) > 0:
                    residue["ATOM"].append(l)
                elif l.find("BOND") == 0:  # len(tmp) > 0:
                    residue["BOND"].append(l)
                elif l.find("ANGL") == 0:  # len(tmp) > 0:
                    residue["ANGL"].append(l)
                elif l.find("DIHE") == 0:  # len(tmp) > 0:
                    residue["DIHE"].append(l)
                elif l.find("IMPH") == 0:  # len(tmp) > 0:
                    residue["IMPH"].append(l)
        self.parse_mass(mass)
        self.residues.append(residue)

    def to_string(self):
        """
        convert the RTF object to string
        """

        strs = "* Topology File.\n"
        strs += "*\n   99   1\n"
        for l, m in zip(self.labels, self.mass):
            strs += f"MASS   -1 {l:4s} {m:12.6f}\n"

        for res in self.residues:
            # print(res)
            strs += "\nRESI {:3s} {:5.3f}\n".format(res["NAME"], res["CHARGE"])
            strs += "GROUP\n"
            for a in res["ATOM"]:
                strs += f"{a:s}\n"
            strs += "\n"
            for b in res["BOND"]:
                strs += f"{b:s}\n"
            strs += "\n"
            for a in res["ANGL"]:
                strs += f"{a:s}\n"
            strs += "\n"
            for d in res["DIHE"]:
                strs += f"{d:s}\n"
            strs += "\n"
            for i in res["IMPH"]:
                strs += f"{i:s}\n"

        return strs

    def to_file(self, filename="test.rtf"):
        """
        Write the PRM file
        """

        f = open(filename, "w")
        f.writelines(self.to_string())
        f.close()

    def parse_mass(self, mass):
        self.labels = []
        self.mass = []
        for a in mass:
            tmp = a.split()
            self.labels.append(tmp[2])
            self.mass.append(float(tmp[3]))

    def merge(self, rtf1=None, single=None):
        """
        merge two RTFs into one
        """

        if rtf1 is not None:
            for l, m in zip(rtf1.labels, rtf1.mass):
                if l not in self.labels:
                    self.labels.append(l)
                    self.mass.append(m)

            self.residues += rtf1.residues
            for i, res in enumerate(self.residues):
                residue = {
                    "NAME": str(i) + res["NAME"][:2],
                    "CHARGE": res["CHARGE"],
                    "ATOM": [],
                    "BOND": [],
                    "ANGL": [],
                    "DIHE": [],
                    "IMPH": [],
                }
                for a in res["ATOM"]:
                    tmp = a.split()
                    a1 = str(i) + tmp[1]
                    a = f"ATOM {a1:6s} {tmp[2]:2s}{float(tmp[3]):12.6f}"
                    residue["ATOM"].append(a)
                for a in res["BOND"]:
                    tmp = a.split("!")
                    tmp1 = tmp[0].split()
                    a1, a2 = str(i) + tmp1[1], str(i) + tmp1[2]
                    a = f"BOND {a1:6s} {a2:6s}"
                    if len(tmp) > 1:
                        a += f"!{tmp[-1]:12s}"
                    residue["BOND"].append(a)
                for a in res["ANGL"]:
                    tmp = a.split("!")
                    tmp1 = tmp[0].split()
                    a1, a2, a3 = str(
                        i) + tmp1[1], str(i) + tmp1[2], str(i) + tmp1[3]
                    a = f"ANGL {a1:6s} {a2:6s} {a3:6s} "
                    if len(tmp) > 1:
                        a += f"!{tmp[-1]:12s}"
                    residue["ANGL"].append(a)
                for a in res["DIHE"]:
                    tmp = a.split("!")
                    tmp1 = tmp[0].split()
                    a1, a2, a3, a4 = (
                        str(i) + tmp1[1],
                        str(i) + tmp1[2],
                        str(i) + tmp1[3],
                        str(i) + tmp1[4],
                    )
                    a = f"DIHE {a1:6s} {a2:6s} {a3:6s} {a4:6s} "
                    if len(tmp) > 1:
                        a += f"!{tmp[-1]:12s}"
                    residue["DIHE"].append(a)
                for a in res["IMPH"]:
                    tmp = a.split("!")
                    tmp1 = tmp[0].split()
                    a1, a2, a3, a4 = (
                        str(i) + tmp1[1],
                        str(i) + tmp1[2],
                        str(i) + tmp1[3],
                        str(i) + tmp1[4],
                    )
                    a = f"IMPH {a1:6s} {a2:6s} {a3:6s} {a4:6s}"
                    if len(tmp) > 1:
                        a += f"!{tmp[-1]:12s}"
                    residue["IMPH"].append(a)
                self.residues[i] = residue

        elif single is not None:
            count = str(len(self.residues))
            l = single["label"]
            self.labels.append(l)
            self.mass.append(single["mass"])
            self.residues += [
                {
                    "NAME": single["name"],
                    "CHARGE": single["charge"],
                    "ATOM": ["ATOM {:6s} {:2s}{:12.6f}".format(count + l, l, float(single["charge"]))],
                    "BOND": [],
                    "ANGL": [],
                    "DIHE": [],
                    "IMPH": [],
                }
            ]


class PRM:
    """
    Parse, convert and merge PRM files for CHARMM
    """

    def __init__(self, input):
        self.keywords = ["BOND", "ANGLE", "DIHEDRAL", "IMPHI", "NONBOND"]
        if type(input) == str:
            print("Converting PRM from the file", input)
            f = open(input)
            prm = f.read().split("\n")
            f.close()
        else:
            prm = input

        self.dict = self.parse(prm)

    def parse(self, prm):
        """
        parse prm list to get the dictionary
        """
        dicts = {}
        for key in self.keywords:
            dicts[key] = []

        do_bond = False
        do_angle = False
        do_dihedral = False
        do_imphi = False
        do_nonbond = False

        for l in prm:
            tmp = l.split()
            if l.find("ANGLE") == 0:
                do_angle = True
            elif l.find("BOND") == 0:
                do_bond = True
            elif l.find("DIHEDRAL") == 0:
                do_dihedral = True
            elif l.find("IMPHI") == 0:
                do_imphi = True
            elif l.find("NONBOND") == 0:
                do_nonbond = True

            elif do_bond:
                if len(tmp) > 0:
                    dicts["BOND"].append(l)
                else:
                    do_bond = False

            elif do_angle:
                if len(tmp) > 0:
                    dicts["ANGLE"].append(l)
                else:
                    do_angle = False

            elif do_dihedral:
                if len(tmp) > 0:
                    dicts["DIHEDRAL"].append(l)
                else:
                    do_dihedral = False

            elif do_imphi:
                if len(tmp) > 0:
                    dicts["IMPHI"].append(l)
                else:
                    do_imphi = False

            elif do_nonbond:
                if len(tmp) > 0:
                    if len(tmp) == 7 and tmp[0] != "!":
                        dicts["NONBOND"].append(l)
                else:
                    do_nonbond = False

        return dicts

    def to_string(self):
        """
        convert the PRM object to string
        """

        strs = "* Force Field Parameter File.\n"
        strs += "*"
        for key in self.keywords:
            if key == "NONBOND":
                strs += "\n\nNONBONDED  NBXMOD 5  GROUP SWITCH CDIEL - \n"
                strs += "CUTNB 14.0  CTOFNB 12.0  CTONNB 10.0  EPS 1.0  E14FAC 0.83333333  WMIN 1.4\n"
                strs += "!                Emin     Rmin/2              Emin/2     Rmin  (for 1-4's)\n"
                strs += "!             (kcal/mol)    (A)\n"
            else:
                strs += f"\n\n{key:s}\n"
            for l in self.dict[key]:
                strs += f"{l:s}\n"

        return strs

    def to_file(self, filename="test.prm"):
        """
        Write the PRM file
        """

        f = open(filename, "w")
        f.writelines(self.to_string())
        f.close()

    def merge(self, prm1=None, single=None):
        """
        merge PRM into one
        """
        if prm1 is not None:
            for key in self.keywords:
                for data in prm1.dict[key]:
                    if data not in self.dict[key]:
                        self.dict[key].append(data)
        elif single is not None:
            self.dict["NONBOND"] += [single["nonbond"]]
            # add the nonbonded parameters


if __name__ == "__main__":
    from pyxtal.db import database

    w_dir = "tmp"
    if not os.path.exists(w_dir):
        os.makedirs(w_dir)

    db = database("benchmarks/test.db")
    row = db.get_row("ACSALA")
    struc = db.get_pyxtal("ACSALA")
    c_info = row.data["charmm_info"]
    prm = open(w_dir + "/pyxtal.prm", "w")
    prm.write(c_info["prm"])
    prm.close()
    rtf = open(w_dir + "/pyxtal.rtf", "w")
    rtf.write(c_info["rtf"])
    rtf.close()

    calc = CHARMM(struc, atom_info=c_info, folder=w_dir)
    print(calc.structure.lattice)
    calc.run(clean=False)
    print(calc.structure.energy)
    print(calc.structure.lattice)
