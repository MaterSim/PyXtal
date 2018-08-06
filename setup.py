from distutils.core import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    include_package_data=True,
    name="pyxtal",
    version="0.1dev",
    author="Scott Fredericks, Qiang Zhu",
    author_email="fredes3@unlv.nevada.edu",
    description="Python code for ab initio generation of crystal structures based on symmetry constraints.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/qzhu2017/pyXtal",
    packages=['pyxtal', 'pyxtal.database'],
    package_data={'pyxtal': ['pyxtal.database/*.csv']},
    classifiers=(
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ),
    requires=['spglib', 'pymatgen', 'numpy', 'scipy', 'ase'],
)
