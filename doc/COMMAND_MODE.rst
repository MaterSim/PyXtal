Run PyXtal executables
==============================

Currently, we provide several utilities to the users so that they can run the code from command line with Python scripting. 
They include:

- ``Pyxtal``: a tool to generate atomic crystals
- ``Pyxtal_m``: a tool to generate molecular crystals
- ``Pyxtal_test``: a tool to test all modules

All of them can be accessed by invoking the ``-h`` command:
::

    $ pyxtal -h

                 ______       _    _          _   
                (_____ \     \ \  / /        | |   
                 _____) )   _ \ \/ / |_  ____| |  
                |  ____/ | | | )  (|  _)/ _  | | 
                | |    | |_| |/ /\ \ |_( ( | | |___
                |_|     \__  /_/  \_\___)_||_|_(___
                       (____/      
    
    
    ----------------------(version 0.0.1 )----------------------
    
    A Python package for random crystal generation
    The source code is available at https://github.com/qzhu2017/pyxtal
    Developed by Zhu's group at University of Nevada Las Vegas
    
    
    Usage: pyxtal [options]
    
    Options:
      -h, --help            show this help message and exit
      -s sg, --symmetry=sg  desired symmetry, number or string, e.g., 36, Pbca, Ih
      -e element, --element=element
                            desired elements: e.g., Li
      -n numIons, --numIons=numIons
                            desired numbers of atoms: 16
      -f factor, --factor=factor
                            volume factor: default 1.0
      -v verbosity, --verbosity=verbosity
                            verbosity: default 0; higher values print more
                            information
      -a attempts, --attempts=attempts
                            number of crystals to generate: default 1
      -o outdir, --outdir=outdir
                            Directory for storing output cif files: default 'out'
      -d dimension, --dimension=dimension
                            desired dimension: (3, 2, 1, 0): default 3
      -t thickness, --thickness=thickness
                            Thickness, in Angstroms, of a 2D crystal, or area of a
                            1D crystal, None generates a value automatically:
                            default None
    
