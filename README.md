# JWST_Reffiles

[![Build Status](https://travis-ci.org/spacetelescope/jwst_reffiles.svg?branch=master)](https://travis-ci.org/spacetelescope/jwst_reffiles)

This repository contains code that can be used to generate JWST reference files
used in the data calibration pipeline. The output reference files are in CRDS format
and should be immediately useable in the pipeline.

## Installation

To install `jwst_reffiles`, first clone the repository:
`git clone https://github.com/spacetelescope/jwst_reffiles.git`

Then, install the package:
```
cd jwst_reffiles
pip install .
```

## Dependencies

JWST calibration pipeline:

* [JWST calibration pipeline][d5].

[d5]: https://github.com/spacetelescope/jwst


### Contributing

Prior to contributing to the `jwst_reffiles` development, please review our [style guide](https://github.com/spacetelescope/jwst_reffiles/blob/master/style_guide/style_guide.md).

The following is a bare bones example of a best work flow for contributing to the project:

1. Create a fork off of the `spacetelescope` `jwst_reffiles` repository.
2. Make a local clone of your fork.
3. Ensure your personal fork is pointing `upstream` properly.
4. Create a branch on that personal fork.
5. Make your software changes.
6. Push that branch to your personal GitHub repository (i.e. `origin`).
7. On the `spacetelescope` `jwst_reffiles` repository, create a pull request that merges the branch into `spacetelescope:master`.
8. Assign a reviewer from the team for the pull request.
9. Iterate with the reviewer over any needed changes until the reviewer accepts and merges your branch.
10. Delete your local copy of your branch.


## Code of Conduct

Users and contributors to the `jwst_reffiles` repository should adhere to the [Code of Conduct](https://github.com/spacetelescope/jwst_reffiles/blob/master/CODE_OF_CONDUCT.md).  Any issues or violations pertaining to the Code of Conduct should be brought to the attention of a `jwst_reffiles` team member or to `conduct@stsci.edu`.


## Questions

Any questions about the `jwst_reffiles` project or its software can be directed to `arest@stsci.edu` or `hilbert@stsci.edu`.


## Current Development Team
- Armin Rest [@arest](https://github.com/arest)
- Bryan Hilbert [@bhilbert4](https://github.com/bhilbert4)
- Alicia Canipe [@aliciacanipe](https://github.com/aliciacanipe)
