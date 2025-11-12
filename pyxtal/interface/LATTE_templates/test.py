from pyxtal import pyxtal
from ase.io import read
from latte_calc import latte_single_optimize

#s = pyxtal()
#s.from_random(3, 19, ["C"], [4])
# Load CIF file into ASE Atoms
s = read("alpha_match.cif")

s_opt, epa, t, err = latte_single_optimize(
    s,
    exe="./LATTE_DOUBLE",
    workdir=".",
    log_name="latte.log",
    control_relax=1,        # RELAX=1
    control_maxiter=30,    # MAXITER=300
    control_rlxftol=0.1, # RLXFTOL=1.0E-03
    control_restart=0,      # RESTART=1
)

print("Error:", err)
print("Energy per atom:", epa)
print("Optimized structure (pyxtal):", s_opt is not None)


'''
from pyxtal import pyxtal
from latte_calc import LATTECalc

s = pyxtal()
s.from_random(3, 19, ["C"], [4])

calc = LATTECalc(
    s,
    exe="./LATTE_DOUBLE",
    workdir=".",
    log_name="latte.log",
    control_relax=1,
    control_maxiter=300,
    control_rlxftol=1.0e-3,
    control_restart=1,
)
calc.run()

print("Error:", calc.error)
print("Energy per atom:", calc.energy_per_atom)
print("Optimized structure:", calc.optimized)
'''