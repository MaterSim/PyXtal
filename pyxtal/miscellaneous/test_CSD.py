from pyxtal import pyxtal
s = pyxtal(molecular=True)
codes = ['CABCOB02','CEFMEN',  'DACSOP',  'EHIVIC',  'EWETOQ',  
         'JEFWAS',  'KIVXEU',  'KUCPON',  'MUMMOY',  'OCEPUK', 
         'TABREX01','TBTYAC01','XEZGUF','ZIWTUX']
for code in codes:
    s.from_CSD(code)
    sph = s.get_spherical_images()
    sph.plot_sph_images(figname='t.png', molecule=True)
