import os
from pyxtal.db import database
from pyxtal.interface.charmm import CHARMM
from ost.parameters import ForceFieldParameters

# Get pyxtal
db = database('pyxtal/database/test.db')
xtal = db.get_pyxtal('ACSALA') # missing data??
smiles = [mol.smile for mol in xtal.molecules]

# setup simulation
for i, method in enumerate(['am1bcc']):
    print('\n', method)
    dir_name = 'test1'
    if not os.path.exists(dir_name): os.mkdir(dir_name)
    params = ForceFieldParameters(smiles, 'gaff')
    parameters = params.params_init.copy()

    os.chdir(dir_name)
    xtal.to_file('init.cif')
    atoms = xtal.get_forcefield('gaff', code='charmm', chargemethod=method, parameters=parameters)
    atom_info = atoms.get_atom_info()
    calc = CHARMM(xtal, prefix='charmm', atom_info=atom_info)
    #print(calc.structure.lattice)
    calc.run() #clean=False)
    print(calc.structure.energy)
    print(calc.structure.lattice)
    calc.structure.to_file('opt.cif')
