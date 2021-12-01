import warnings
warnings.filterwarnings('ignore')

from pyxtal import pyxtal
from pyxtal.util import search_csd_entries_by_code

c = pyxtal(molecular=True)
for csd in ['ACSALA', 'TROXAN', 'TRIZIN', 'HXMTAM', 'PYRZIN', 
            'PYRZOL', 'CYHEXO', 'CYTSIN', 'ADAMAN01', 'IMAZOL01', 
            'URACIL', 'ACETAC', 'CYANAM01', 'FORMAM', 'SUCACB02',
            'ECARBM01', 'XAFQAZ', 'KONTIQ', 'XATJOT', 'XAFQON']:
    codes = search_csd_entries_by_code(csd)
    for code in codes:
        c.from_CSD(code)
        print(c)
