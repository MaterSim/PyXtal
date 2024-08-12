from distutils.core import setup
from os import path

import setuptools  # noqa

this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, "README.md"), encoding="utf-8") as f:
    long_description = f.read()

with open("README.md") as fh:
    long_description = fh.read()

setup(
    name="pyxtal",
    version="1.0.0",
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
        "pyxtal.database": ["*.csv", "*.json", "*.db"],
        "pyxtal.database.cifs": ["*.cif", "*.vasp"],
        "pyxtal.potentials": ["*"],
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
        "pymatgen>=2024.3.1",
        "pandas>=2.0.2",
        "networkx>=2.3",
        "ase>=3.23.0",
        "scipy>=1.7.3",
        "numpy>=1.26,<2",  # prevent the use of numpy2
        "importlib_metadata>=1.4",
        "typing-extensions>=4.12",
    ],
    extras_require={
        "visualization": ["py3Dmol>=0.8.0"],
        "descriptor": ["pyshtools>=4.10.3"],
        "molecules": ["openbabel", "pybel"],
        "test": ["wheel", "pytest", "coverage", "pytest-cov", "monty>=2024.2.26"],
    },
    python_requires=">=3.9",
    license="MIT",
)
