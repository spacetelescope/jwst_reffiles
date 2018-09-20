#!/usr/bin/env python
import os
import subprocess
import sys
from setuptools import setup, find_packages, Extension, Command
from setuptools.command.test import test as TestCommand


# allows you to build sphinx docs from the package
# main directory with "python setup.py build_sphinx"

try:
    from sphinx.cmd.build import build_main
    from sphinx.setup_command import BuildDoc

    class BuildSphinx(BuildDoc):
        """Build Sphinx documentation after compiling C source files"""

        description = 'Build Sphinx documentation'

        def initialize_options(self):
            BuildDoc.initialize_options(self)

        def finalize_options(self):
            BuildDoc.finalize_options(self)

        def run(self):
            build_cmd = self.reinitialize_command('build_ext')
            build_cmd.inplace = 1
            self.run_command('build_ext')
            build_main(['-b', 'html', './docs', './docs/_build/html'])

except ImportError:
    class BuildSphinx(Command):
        user_options = []

        def initialize_options(self):
            pass

        def finalize_options(self):
            pass

        def run(self):
            print('!\n! Sphinx is not installed!\n!', file=sys.stderr)
            exit(1)


class PyTest(TestCommand):
    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = ['tests/test_calib_prep.py']
        self.test_suite = True

    def run_tests(self):
        # import here, cause outside the eggs aren't loaded
        import pytest
        errno = pytest.main(self.test_args)
        sys.exit(errno)


# make sure jwst is available
try:
    import jwst
except ImportError:
    try:
        subprocess.check_call(['git', 'clone',
                               'https://github.com/spacetelescope/jwst.git'])
        sys.path.insert(1, 'jwst')
        # import jwst
    except subprocess.CalledProcessError as e:
        print(e)
        exit(1)


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
    keywords=['astronomy'],
    classifiers=[
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering :: Astronomy',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ],
    packages=find_packages(exclude=["examples"]),
    install_requires=[
        'astropy>=1.2',
        'numpy>=1.9',
        'matplotlib>=1.4.3',
        'asdf>=2.1.0',
        'scipy>=0.17',
        'gwcs>=0.9'
    ],
    include_package_data=True,
    cmdclass={
        'test': PyTest,
        'build_sphinx': BuildSphinx
    },)
