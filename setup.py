#!/usr/bin/env python
from setuptools import setup, find_packages


setup(
    name='jwst_reffiles',
    version='0.0.0',
    description='Create JWST reference files',
    long_description=('A tool to create CRDS-formatted reference files'
                      'for JWST from a set of input dark current files'
                      'and a set of flat field files.'),
    author='STScI (Rest, Hilbert, Canipe, et al.)',
    author_email='arest@stsci.edu',
    url='https://github.com/spacetelescope/jwst_reffiles',
    license="BSD-3-Clause",
    classifiers=[
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering :: Astronomy',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ],
    packages=find_packages(),
    include_package_data=True,
    scripts=["jwst_reffiles/mkrefs.py"],
    install_requires=[
        'astropy>=4.0',
        'jwst',
        'numpy>=1.16',
        'matplotlib>=1.4.3',
        'scipy>=1.1',
    ],
    extras_require=dict(
        docs=[
            "nbsphinx",
            "ipykernel",
            "sphinx",
            "sphinx_rtd_theme",
        ],
        test=[
            "pytest",
        ]
    ),
    )
