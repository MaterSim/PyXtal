import psutil
import os

def print_memory_usage():
    process = psutil.Process(os.getpid())
    mem_info = process.memory_info()
    print(f"Memory usage: {mem_info.rss / 1024 ** 2:.2f} MB")  # RSS is the Resident Set Size

if __name__ == "__main__":
    print_memory_usage()
    from pyxtal.symmetry import Group
    print_memory_usage()
    g = Group(4, quick=True)
    print_memory_usage()
    g = Group(227)
    print_memory_usage()

    from pyxtal import pyxtal
    xtal = pyxtal()
    xtal.from_spg_wps_rep(194, ['2c', '2b'], [2.46, 6.70])
    print_memory_usage()

    xtal.subgroup_once()
    print_memory_usage()
