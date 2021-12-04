import warnings
warnings.filterwarnings('ignore')

from ase.db import connect
from pyxtal import pyxtal
from pyxtal.util import search_csd_entries_by_code

strucs = []
for csd in [#'GLYCIN01',  'HUSVUM', 'HIHDUZ', 'DOPSIK', 'TABQOG03', 'DUPXOD', 
            #'GLYMSN10', 'GOBFAF', 'GOBZOO', 'ARAZEZ','HMBENZ21', 'BEMQUF', 'KUZJIB', 
            #NUMRET long
            'SUZGUS', 'MABCOL', 'UMEQUB', 'XACXOS', 'AMBZPH',
            'AMISAS', 'APAPOY', 'HIWZIZ', 'OHUXEV', 'QOJPAI',
            'HOWTOF', 'NUGCOI', 'NUJMAH', 'NUNJIP', 'ZUZCEF',
            'VITRUL', 'VOLCOP', 'QOLXEX02', 'UKULAQ',
            'ACEMID02', 'AFIPIP', 'BUTHEE', 'CELKEK', 'HAHDID',
            'HEXWIQ01', 'HIFWOJ', 'HIFZIF', 'HIRBAM', 'EHINOB',
            'HUKYIW', 'HYQUIN06', 'ISIVIR', 'JATBIQ', 'JUDBUI',
            'AFUVAZ', 'BZCBNL01', 'NUDREK',
            'HEVRUV', 'JAPCIM', 'TCYETY02', 'ZZZWOU01', 'ELIFOX',
            'GADLOQ', 'GAJLAI', 'MAMFIT', 'OKIXIS', 'NOLRUC',
            'PACNOA', 'SABQEV', 'SAJWEJ', 'HAMQOC', 'NUMRET',
            'HCHXDO', 'HETVUY', 'HIRYOY', 'HIYLOQ01', 'HONWIQ'
            'TROXAN', 'PYRZIN', 'ACETAC', 'ADAMAN01', 
            'TRIZIN', 'HXMTAM', 'PYRZOL', 'CYHEXO', 'CYTSIN', 
            'IMAZOL01', 'URACIL', 'CYANAM01', 'FORMAM', 'SUCACB02',
            'ECARBM01', 'XAFQAZ', 'KONTIQ', 'XATJOT', 'XAFQON'
            ]:
    codes = search_csd_entries_by_code(csd)
    for code in codes:
        c = pyxtal(molecular=True)
        c.from_CSD(code)
        strucs.append(c)
        print(c)
        #if code == 'TRIZIN04': print(c.to_file()); print(c); import sys; sys.exit()

with connect('test.db') as db:
    for xtal in strucs:
        if xtal.tag['ccdc_number'] is None: xtal.tag['ccdc_number']=1240839
        kvp = {
                "csd_code": xtal.tag['csd_code'],
                "mol_smi": xtal.tag['smiles'],
                "ccdc_number": xtal.tag['ccdc_number'],
                #"publication": xtal.tag['publication'],
                }
        print(kvp)
        db.write(xtal.to_ase(), key_value_pairs=kvp)
