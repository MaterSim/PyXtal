from setuptools import setup
from setuptools.command.install import install
import os

class CustomInstallCommand(install):
    def run(self):
        # Check if the custom environment variable is set
        if os.getenv("INSTALL_JULIA") == "1":
            print("CustomInstallCommand: Installing Julia packages...")
            try:
                from juliacall import Main as jl
                jl.seval("using Pkg")
                jl.seval('Pkg.add(Pkg.PackageSpec(name="CrystalNets", uuid="7952bbbe-a946-4118-bea0-081a0932faa9"))')
                print("CrystalNets has been installed successfully.")
            except ImportError:
                print("juliacall is not installed. Skipping Julia package installation.")
            except Exception as e:
                print(f"An error occurred during Julia package installation: {e}")
                sys.exit(1)  # Optionally, stop the installation if this step fails
        else:
            print("CustomInstallCommand: Skipping Julia package installation.")

        # Run the standard install process
        install.run(self)

with open("README.md") as fh:
    long_description = fh.read()

setup(
    name="pyxtal",
    version="1.0.7",
    author="Scott Fredericks, Kevin Parrish, Qiang Zhu",
    author_email="alecfans@gmail.com",
    description="Python code for generation of crystal structures based on symmetry constraints.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/MaterSim/PyXtal",
    packages=[
        "pyxtal",
        "pyxtal.database",
        "pyxtal.interface",
        "pyxtal.lego",
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
    setup_requires=["juliacall>=0.9.0"],
    install_requires=[
        "spglib>=2.5.0",
        "pymatgen>=2024.3.1",
        "pandas>=2.0.2",
        "networkx>=2.3",
        "ase>=3.23.0",
        "scipy>=1.7.3",
        #"numpy>=1.26,<2",  # prevent the use of numpy2
        "vasprun-xml>=1.0.4",  # prevent the use of numpy2
        "importlib_metadata>=1.4",
        "typing-extensions>=4.12",
        "pyocse>=0.1.1",
        "psutil",
    ],
    extras_require={
        "visualization": ["py3Dmol>=0.8.0"],
        "descriptor": ["pyshtools>=4.10.3"],
        "molecules": ["openbabel", "pybel"],
        "test": ["wheel", "pytest", "coverage", "pytest-cov", "monty>=2024.2.26"],
    },
    cmdclass={
        'install': CustomInstallCommand,
    },
    python_requires=">=3.9",
    license="MIT",
)

