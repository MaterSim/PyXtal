import setuptools
#from distutils.core import setup

from glob import glob

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
#setup(
    include_package_data=True,
    name="crystallography",
    version="0.1dev",
    author="Scott Fredericks, Qiang Zhu",
    author_email="fredes3@unlv.nevada.edu",
    description="Python code for ab initio generation of crystal structures based on symmetry constraints.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/qzhu2017/crystallography",
    packages=setuptools.find_packages(),

    data_data = [('database', 'database/wyckoff_list.csv'), ('database', 'database/wyckoff_generators.csv'), ('database', 'database/wyckoff_symmetry.csv') ],

    package_data={
    'crystallography': [
        'database/*.csv',
        ],
    },

    classifiers=(
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ),
)
