# This is a collections of examples and tutorials to demonstrate how to use PyXtal
```
├── README.md
├── example_01_3D_VASP.py   
├── example_02_LJ_cluster.py
├── example_03_LJ_optalg.py
├── example_04_LJ_38.py
├── example_05_LJ_4D.py
├── example_06_C_2D_lammps.py
├── example_07_3D_ICE_lammps.py
└── tutorials_notebook
    ├── 01_atomic_crystals.ipynb
    └── 02_molecular_crystals.ipynb
```

# Installation
Examples 06 and 07 require the installation of LAMMPS and GULP as follows

### LAMMPS
download the most recent version of lammps (e.g., `lammps-11Aug17`) from http://lammps.sandia.gov/download.html

```
$ tar -xf lammps-stable.tar.gz
$ cd lammps-11Aug17/src
$ make yes-class2
$ make yes-kspace
$ make yes-manybody
$ make yes-molecule
$ make yes-python
$ make yes-rigid
$ make yes-user-misc
$ make mpi
$ make mpi mode=shlib
```
It works up to `LAMMPS (3 Mar 2020)`
This should create the file `liblammps_mpi.so` in the src directory, as well as a soft link liblammps.so, which is what the Python wrapper will load by default.
Then one just need to add the path of src to the `.bashrc` file as follows,

```
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/scratch/qzhu/soft/lammps/src #please modify the path if necessary
```
then in lammps/src

```
make install-python
```

### GULP
Download GULP from http://gulp.curtin.edu.au/gulp/
After the installation, export the following environment variables
```
$ export PATH=/scratch/qzhu/pkgs/gulp-5.2/Src:$PATH
$ export GULP_LIB=/scratch/qzhu/pkgs/gulp-5.2/Libraries
```

### Test
To test if it works correctly, run the following,
```
$ python test_installation.py 
PyXtal:  0.0.9
ase:  3.19.1
Using PyXtal to generate structure
Convert PyXtal structure to ASE
launch the LAMMPS calculator
-5.348217915131606
launch the GULP calculator
-5.34821791
```
