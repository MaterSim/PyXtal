Installation and Setup
======================

Dependencies
------------

Versions indicated are those used during development. Other versions may be compatible, but have not yet been tested.

  * `SciPy 1.0.1 <https://www.scipy.org/install.html>`_  
  * `NumPy 1.14.3 <https://www.scipy.org/scipylib/download.html>`_  
  * `Pandas 0.20.3 <https://pandas.pydata.org/getpandas.html>`_  
  * `Pymatgen 2017.9.3 <http://pymatgen.org/#getting-pymatgen>`_  
  * `SpgLib for Python 1.9.9.44 <https://atztogo.github.io/spglib/python-spglib.html#installation>`_  

Optional:

  * `Openbabel 2.4.1 <http://openbabel.org/wiki/Category:Installation>`_  

Openbabel is not necessary, and only adds additional file format support for importing molecules. You must install the C++ pacakge before installing the Python bindings. For Debian based systems, your distribution may already have installable packages:

``sudo apt-get install openbabel``  

``pip install openbabel``

Also note that the openbabel Python bindings require swig to install:

``sudo apt-get install swig``  

Older version of swig (before 2.0) will not work. For other systems, you must compile the openbabel Python bindings yourself. There are tutorials for this on the `openbabel website
<https://openbabel.org/docs/dev/UseTheLibrary/PythonInstall.html>`_, as well as in the `pymatgen documentation
<http://pymatgen.org/installation.html#openbabel-mac-os-x-tested-on-v2-3-2>`_.

Installation
------------

To install PyXtal, first install all dependencies, then make a copy of the source code:

``git clone https://github.com/qzhu2017/pyxtal``

Then, inside of the downloaded directory, run

``python setup.py install``

This will install the module. The code can be used within Python via

``pyxtal_test``

You expect to see the following output.

::

     ______       _    _          _   
    (_____ \     \ \  / /        | |   
     _____) )   _ \ \/ / |_  ____| |  
    |  ____/ | | | )  (|  _)/ _  | | 
    | |    | |_| |/ /\ \ |_( ( | | |___
    |_|     \__  /_/  \_\___)_||_|_(___
           (____/      


------------------------(version 0.0.1 )----------------------

A Python package for random crystal generation
The source code is available at https://github.com/qzhu2017/pyxtal
Developed by Zhu's group at University of Nevada Las Vegas


====== Testing functionality for pyXtal version 0.1dev ======
Importing sys...
Success!
Importing numpy...
Success!
Importing pymatgen...
Success!
Importing pandas...
Success!

  
More extensive test can be invoked by running

``pyxtal_test -m all``

Ideally, one should see the completion of all modules in the end.
