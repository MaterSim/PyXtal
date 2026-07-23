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
├── example_08_QRS_known_cell.py
├── example_09_QRS_conf.py
├── example_10_mlp_relax.py
└── tutorials_notebook
    ├── 01_atomic_crystals.ipynb
    └── 02_molecular_crystals.ipynb
    └── 03_pxrd.ipynb
    └── 04_box.ipynb
    └── 05-crystal-packing.ipynb
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

# Sample script to run example_09 and example_10
```
#!/bin/sh -l
#SBATCH --partition=Apus
#SBATCH -J 0822-Z2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --time=72:00:00
export OMP_NUM_THREADS=1
nproc=$SLURM_CPUS_PER_TASK

# Same settings as Tests-0807/CEQGEL:
#   no ASU / no fix-translation, dl=1.2, soft-clash -0.5, default max-grid-product 1e9
#python example_09_QRS_conf.py --code ACRDIN05 BOQQUT01 CEQGEL FLUANT JOHSOP UJIRIO02 XAFPAY02 XATMOV ZZZDKE01 \
#python example_09_QRS_conf.py --code ACSALA ACSALA13 ADAMAN01 AFIGIH AVIBEN AXIDER AXOSOW01 BENZEN BETFOV BIPHEN BOQQUT BZPHAN01 \
#python example_09_QRS_conf.py --code HAJJIN HAMTIZ01 JAYDUI JOHSOP JUFRIO KEKJEQ KONTIQ KONTIQ09 KUPWOJ KUPWUP LEVJON MERQIM MERQOS MERRAF MIVDEC MUVMIA \
#python example_09_QRS_conf.py --code NACJAF NAPHTA15 NICLAN OBEQET OBEQIX OBEQOD OBEQUJ OKUPUG OXALAC02 OXALAC11 PAHYON01 PHENAZ01 PYRZIN01 QAXMEH53 QQQCIG04 QUATER10 QUPHEN \
#python example_09_QRS_conf.py --code RESORA03 SIRMIQ01 SITJUC SUCXIZ TBZPYR TIDFES UJIRIO01 UJIRIO05 UREAXX02 UVAGUU UVAHAB UVAHEF VOBYAN WEXBOS WICZUF WIDBAO \
#python example_09_QRS_conf.py --code XAFPAY XAFPAY01 XAFPAY03 XAFQAZ XAFQIH XAFQON XATJOT XELYUJ XULDUD XULDUD01 YIHVUI YOKBIK \
#python example_09_QRS_conf.py --code CAYKUJ COUMAR01 CRYSEN CYANAM01 DBZCOR DOCXEA DURNAH DUTGIK ECENAD EVIDEV FIQSEG FOJTUU FOJVAC FORMAM FUNZOE GUFJOG \
#python example_09_QRS_conf.py --code ACEMID02 MERQUY PYRENE07 \
#  --out-dir Tests-0822 \
#  --nproc $nproc --ngen 200 --npop 96 #--no-check-stable #--no-fix-translation --no-asu-clamp --no-check-stable #--delta-length 1.0 #--no-soft-clash #--no-order-identical-sites
python example_10_mlp_relax.py Tests-0822 --calculator MACEOFF --max-unique 100
```
