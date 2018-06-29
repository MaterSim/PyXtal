from structure import get_wyckoff_generators
from structure import get_wyckoffs
from structure import are_equal

for sg in range(1, 231):
    wyckoffs = get_wyckoffs(sg)
    generators = get_wyckoff_generators()
    for wp_index, wp in enumerate(wyckoffs):
        generators = g
