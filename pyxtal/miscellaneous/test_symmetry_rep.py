from pyxtal.symmetry import Group
from pyxtal.operations import OperationAnalyzer
for g in range(1, 231):
    spg = Group(g)
    print('\n\n\n', spg.number, spg.symbol)
    ss = spg.get_spg_symmetry_object()
    ss.to_one_hot(verbose=True)
    for wp in spg:
        ss = wp.get_site_symmetry_object()
        print('\n', g, wp.get_label(), ss.name)
        ss.to_one_hot(verbose=True)
