PDELab
======

This is the development version of the dune-pdelab tutorials.

[PDELab][0] is a PDE solver toolbox built on top of DUNE, the [Distributed and Unified
Numerics Environment][1].

This package contains a number of individual tutorials that present different aspects of
the PDELab framework.

If you need help, please ask on our [mailinglist][2]. Bugs can also be submitted
to the [bugtracker][3] instead.

Dependencies
------------

The current development version of dune-pdelab-tutorials depends on the following software:

* PDELab and DUNE core libraries (dune-common, dune-geometry, dune-grid, dune-istl,
  dune-localfunctions, dune-typetree, dune-functions, dune-pdelab) version 2.5.0, and their respective
  dependencies.

* The [dune-uggrid][4] module (specifically the releases/2.5 branch).
  This supersedes the old external library UG.

* The [dune-alugrid][5] library (specifically the releases/2.5 branch), and its respective
  dependencies.

* A compiler with support for C++14, at least GCC 4.9 or clang 3.5.

* CMake 3.1.

License
-------

All of the contained tutorials and example programs are free open-source documentation
and software.

See the file [LICENSE.md][6] for full copying permissions.

Installation
------------

For installation instructions please see the [DUNE website][1].

Links
-----

 [0]: http://www.dune-project.org/pdelab/
 [1]: http://www.dune-project.org
 [2]: http://lists.dune-project.org/mailman/listinfo/dune-pdelab
 [3]: https://gitlab.dune-project.org/pdelab/dune-pdelab-tutorials/issues
 [4]: https://gitlab.dune-project.org/staging/dune-uggrid
 [5]: https://gitlab.dune-project.org/extensions/dune-alugrid
 [6]: LICENSE.md
