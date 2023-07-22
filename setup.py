from distutils.core import setup
import setuptools  # noqa
from os import path

this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, "README.md"), encoding="utf-8") as f:
    long_description = f.read()

exec(open("pyxtal/version.py").read())

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="pyxtal",
    version=__version__,
    author="Scott Fredericks, Qiang Zhu",
    author_email="qiang.zhu@unlv.edu",
    description="Python code for generation of crystal structures based on symmetry constraints.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/qzhu2017/PyXtal",
    packages=[
        "pyxtal",
        "pyxtal.database",
        "pyxtal.interface",
        "pyxtal.optimize",
        "pyxtal.potentials",
        "pyxtal.database.cifs",
    ],
    package_data={
        "pyxtal.database": ["*.csv", "*.json"],
        "pyxtal.database.cifs": ["*.cif", "*.vasp"],
        'pyxtal.potentials': ['*'],
    },
    scripts=[
        "scripts/pyxtal_main.py",
        "scripts/pyxtal_symmetry.py",
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        "spglib>=1.10.4",
        "pymatgen>=2022.0.17",
        "pandas>=0.24.2",
        "networkx>=2.3",
        "py3Dmol>=0.8.0",
        'ase>=3.18.0',  #covered by pymatgen
        'numba>=0.55.2', #now supports numpy 1.22
        'scipy>=1.7.3',
        'importlib_metadata>=1.4',
        'pyshtools>=4.10.3',
        #"openbabel>=3.0.0",
    ],
    python_requires=">=3.7, <=3.12", #add the restriction for now issue #189
    license="MIT",
)
