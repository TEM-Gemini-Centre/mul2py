from setuptools import setup, find_packages

exec(open('mul2py/release_info.py').read()) #collect version info

setup(
    name=name,
    version=str(version),
    license=license,
    author=author,
    author_email=email,
    description="Conversion of MULTEM Multislice image simulation results from MATLAB to Python.",
    long_description=open('README.md').read(),
    classifiers=[
        "Programming Language :: Python :: 3.7",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Natural Language :: English",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Physics",
    ],
    packages=find_packages(),
    install_requires=[
        "numpy",
        "matplotlib",
        "pathlib",
        "hyperspy",
        "h5py",
        "argparse"

    ],
    package_data={
        "": ["LICENSE", "README.md"],
        "": ["*.py"],
        "": ["*.ipynb"],
    },
)
