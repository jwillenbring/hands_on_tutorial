# Trilinos Hands-on Tutorial #

## Reference and tutorial material ##

  * [Trilinos website](http://trilinos.org) (click on logo above)
  * Each Trilinos package has [Doxygen documentation](http://trilinos.org/packages/).
  * [Trilinos "Getting Started"](http://trilinos.org/?page_id=95) page.
  * See links below to slides and video recordings of tutorials from previous Trilinos User Group meetings.

## Options for running examples ##

<a href='Hidden comment: 
=== Student shell accounts (VECPAR 2014 only) ===

* Tutorial material lives on [https://github.com/jwillenbring/Trilinos_tutorial Github]
* Steps to get started
# Log in to wopr.nic.uoregon.edu
# ssh nucN where N is one of 01, 02, ..., 09, 12, ..., 15 (omit 10 and 11)
# cd Trilinos_tutorial
# source ./setup.sh  (loads MPI and Trilinos modules)
# cd cmake_build
# ./live-cmake  (builds all the examples)
# Change into build subdirectories to run examples by hand
* We may do the last one or two steps for you; stay tuned!
* You may also use the Trilinos tutorial on your computer:
# git clone https://github.com/jwillenbring/Trilinos_tutorial.git
# cd Trilinos_tutorial
# Modify ./live-cmake as necessary for your build environment
# Run ./live-cmake and continue from there

===  HPCLinux virtual machine (VM) ===

* HPCLinux VM with same build environment as the student shell accounts
* We won"t cover this today, but you can download it and try it at home
* *NOTE*: The VM file is *QUITE LARGE* (11116475904 bytes)
* Download locations:
* [http://tau.uoregon.edu/livedvd/latest.ova University of Oregon (HTTP)]
* [ftp://paratools08.rrp.net/hpclinux/latest.ova ParaTools, Inc. (FTP)]
* wopr.nic.uoregon.edu:~livetau/HPCLinux_June14.2.ova
* VM [http://www.paratools.com/HPCLinux setup instructions] from ParaTools
'></a>
### **WebTrilinos**: Web interface to edit, build, and run C++ code ###
  * [Click here](http://trilinos.csbsju.edu/WebTrilinosMPI-shared-12.0/c++/index.html) to access WebTrilinos (Trilinos version 12.0 C++ API).
  * [Click here](http://trilinos.csbsju.edu/WebTrilinosMPI-shared-12.0/index.html) to access WebTrilinos (Trilinos version 12.0, General access).
  * [Click here](http://trilinos.csbsju.edu/WebTrilinosMPI-OpenMP-shared-11.14/c++/index.html) to access WebTrilinos (Trilinos version 11.14).
  * Access to this site is password protected.  Login information will be given during live tutorials as needed.
  * This page gives you a text box in which you can paste, type, or edit C++ code.
    * That code will compile and link against a recent release of Trilinos, and run.
    * The web page will show you the resulting output.
    * Use Ctrl+A to highlight all the example code, Ctrl+C to copy it, and Ctrl+V to paste it in the WebTrilinos window.
    * You can't read or write files, but you can embed input data in your program as a string.
  * WebTrilinos is installed on a server at St. John's University, MN.

### MueLu virtual machine (VM) ###

  * VM and accompanying PDF delivering the [MueLu tutorial](http://trilinos.org/packages/muelu/muelu-tutorial/)

### Build Trilinos yourself on your computer ###

Please read the BuildScript page to learn how to set Trilinos' build configuration options.

  1. Prerequisites
    * C++ and C compiler
    * [CMake](http://www.cmake.org/cmake/resources/software.html) version >= 2.8
    * BLAS and LAPACK libraries
    * MPI (Message Passing Interface) (optional)
  1. Download [Trilinos](http://trilinos.org/download)
  1. Find a configuration script suitable for your system.
  1. Copy the script into your build directory and modify it if necessary.
    * You may modify the script by hand, use the CMake graphical user interface, or use the console-based interface `ccmake`.
  1. Use the script to run CMake.
  1. Run `make` and `make install`.  You may specify the `-j <N>` flag for a parallel build.

Once you have built and installed Trilinos, you may build your program with Trilinos.  You have two options: use CMake, or use a Makefile.
  * CMake: See [our example](https://code.google.com/p/trilinos/wiki/CMakeFindPackageTrilinosExample)
  * Makefile
    1. Get an example Makefile from [here](http://trilinos.org/oldsite/Export_Makefile_example.txt).
      * Learn about the Makefile and the Makefile.export system [here](http://trilinos.org/oldsite/Export_Makefile.txt).
      * You can find the Makefile.export.package\_name files once you have built and installed Trilinos.
      * They will be in the include directory of your installation directory.
      * Example: `TrilinosInstall/include/Makefile.export.Epetra`.
    1. Customize the Makefile to your situation.
      * The example Makefile uses only the Epetra package, but `Makefile.export.$PACKAGE_NAME` files are available for all packages.
      * You may also include `Makefile.export.Trilinos` to get all Trilinos packages.
    1. Add your code into a file called main.cpp and build away!

## Examples illustrating Kokkos thread-scalable expressions to generate Tpetra objects ##

> If using cmake to build your own application on an installed trilinos, use the top level [CMakeLists.txt](CMakeLists_Kokkos01.md) and create two subdirectories for the examples below with the files [CMakeLists.txt file](CMakeLists_Kokkos02.md) for KokkosBasic directory and [CMakeLists.txt file](CMakeLists_Kokkos03.md) for KokkosTpetra directory.

  1. **Learn the basics of [Kokkos](http://trilinos.org/packages/kokkos)**
    * [Lesson 1a](KokkosTutorial01.md): Learn how to Initialize Kokkos using functors.
    * [Lesson 1b](KokkosTutorial01b.md): Learn how to Initialize Kokkos using lambdas.
    * [Lesson 2a](KokkosTutorial02.md): Learn how to perform a reduction using functors.
    * [Lesson 2b](KokkosTutorial02b.md): Learn how to perform a reduction using lambdas.
    * [Lesson 3](KokkosTutorial03.md): Learn how to createa simple view.
    * [Lesson 4](KokkosTutorial04.md): Learn the basic of memory spaces.
    * [Lesson 5](KokkosTutorial05.md): Learn how to use atomics.
  1. **Learn how to create and use [Kokkos](http://trilinos.org/packages/kokkos) with [Tpetra](http://trilinos.org/packages/tpetra).**
    * [Lesson 1](KokkosExample01.md): Learn how to initialize Kokkos.
    * [Lesson 2](KokkosExample02.md): Learn how to write a simple parallel for loop.
    * [Lesson 3](KokkosExample03.md): Learn how to construct a local tridiagonal matrix using a thread-parallel approach.
    * [Lesson 4](KokkosExample04.md): Learn how to construct Tpetra matrices and vectors from existing Kokkos arrays.
    * [Lesson 5](KokkosExample05.md): Learn how to write a simple conjugate gradients solver with Tpetra/Kokkos.

## Examples illustrating the Tpetra-based solver stack ##

  1. **Learn how to create and use [Tpetra](http://trilinos.org/packages/tpetra) dense and sparse linear algebra objects.**
    * [Lesson 1](http://trilinos.org/docs/dev/packages/tpetra/doc/html/Tpetra_Lesson01.html): "Hello world!"  Learn different ways to initialize MPI (or not) and pass an MPI communicator to Tpetra.
    * [Lesson 2](http://trilinos.org/docs/dev/packages/tpetra/doc/html/Tpetra_Lesson02.html): Learn how to make a Tpetra vector, given a communicator from Lesson 1.
    * [Lesson 3](http://trilinos.org/docs/dev/packages/tpetra/doc/html/Tpetra_Lesson03.html): Learn how to implement a simple numerical algorithm (the power method) using Tpetra sparse matrices and vectors.
    * [Lesson 4](http://trilinos.org/docs/dev/packages/tpetra/doc/html/Tpetra_Lesson04.html): Learn different ways to construct a Tpetra sparse matrix.
    * [Lesson 5](http://trilinos.org/docs/dev/packages/tpetra/doc/html/Tpetra_Lesson05.html): Learn how to migrate data in a Tpetra object between two different parallel data distributions.
    * [Advanced exercise](Tpetra_Exercises_Advanced_CrsMatrix_ExplicitTranspose.md): Learn how to compute the explicit transpose of a sparse matrix.
  1. **Learn how to solve linear systems using the [Belos](http://trilinos.org/packages/belos) package of iterative linear solvers, and the [Ifpack2](http://trilinos.org/packages/ifpack2) package of preconditioners**
    * [Create an Ifpack2 preconditioner](Ifpack2CreatePreconditioner.md)
    * [Solve a linear system using Belos and Ifpack2](Tpetra_Belos_CreateSolver.md)

## Examples illustrating the Epetra-based solver stack ##

### Create and use Epetra dense and sparse linear algebra objects ###
  * [Lesson 1](EpetraLesson01.md): "Hello world!"  Learn different ways to initialize MPI (or not) and pass an MPI communicator to Epetra.
  * [Lesson 2](EpetraLesson02.md): Learn how to make an Epetra vector, given a communicator from Lesson 1.
  * [Lesson 3](EpetraLesson03.md): Learn how to implement a simple numerical algorithm (the power method) using Epetra sparse matrices and vectors.
  * [Lesson 4](EpetraLesson04.md): Learn different ways to construct an Epetra sparse matrix.
  * [Lesson 5](EpetraLesson05.md): Learn how to migrate data in an Epetra object between two different parallel data distributions.

### Generate test linear systems using the [Galeri](http://trilinos.sandia.gov/packages/galeri) package ###

  * [Generate a matrix, discretized 2D Laplacian on a Cartesian grid.](GaleriLinearSystem.md)
  * Try generating matrices for some different operators.  The list of supported operators is [here](http://trilinos.org/docs/dev/packages/galeri/doc/html/gl_GalleryCrsMatrix.html).

### Create an algebraic preconditioner using the [Ifpack](http://trilinos.org/packages/ifpack) package ###

  * [Create a preconditioner using the Ifpack preconditioner factory.](IfpackFactory.md)
  * Try generating different preconditioners:
    * The preconditioners supported in the factory interface are: "IC", "ICT", "ILU", "ILUT", and "Amesos".
    * The list of supported parameters for the factory is [here](http://trilinos.org/docs/dev/packages/ifpack/doc/html/index.html#ifp_params).

### Solve a linear system using Amesos ###

The [Amesos](http://trilinos.org/packages/amesos) package
implements direct linear solvers.  It wraps provides one native sparse direct solver called KLU and provides access to other community sparse direct solvers such as [SuperLU](http://crd-legacy.lbl.gov/~xiaoye/SuperLU/).  Amesos accepts Epetra data objects like any Trilinos solver.

**[Solve a linear system using Amesos with SuperLU](AmesosSuperLU.md)**

### Solve a linear system using AztecOO ###

The [AztecOO](http://trilinos.org/packages/aztecoo) package
implements iterative linear solvers.  It wraps an earlier library
which provided its own linear algebra implementation.  !AztecOO can
also work with Epetra matrices and vectors, and any preconditioner
that implements `Epetra_Operator`.  The latter includes Ifpack
preconditioners.

  * [Solve a linear system using Ifpack and AztecOO](IfpackAztecOO.md)

### Solve a linear system using ML and AztecOO ###

The [ML](http://trilinos.org/packages/ml) package implements
multilevel solvers, including algebraic multigrid.  You may use ML's
solvers as preconditioners if you like; this is the way most of our
users use ML.  The following examples show how to do this with the
iterative linear solvers provided by the
[AztecOO](http://trilinos.org/packages/aztecoo/) package.

  * [Use ML as a black-box preconditioner with AztecOO](MLAztecOO.md)
  * [Use ML as a preconditioner with AztecOO, but set some multigrid options](MLAztecOO2.md)

### Solve a linear system using Belos ###

The [Belos](http://trilinos.org/packages/belos) package
implements iterative linear solvers.  Unlike AztecOO, Belos can work
with just about any linear algebra implementation.  Belos also
provides block solvers and other algorithmic optimizations, like
subspace recycling.  Once you learn how to use Belos with Epetra
objects, it is easy to learn how to use Belos with Tpetra or other
linear algebra implementations.

  * [Solve a linear system using Ifpack and Belos](IfpackBelos2.md)

### Solve a nonlinear system using NOX ###

The [NOX](http://trilinos.org/packages/nox) package implements iterative nonlinear solvers.

  * [Example 1: Newton's method](NOXNewton.md)
  * [Example 2: Newton's method](NOXNewton2.md)

### Solve eigenvalue problems using Anasazi ###

The [Anasazi](http://trilinos.org/packages/anasazi) package
implements several different iterative solvers for both standard (A x
= \lambda x) and generalized (K x = \lambda M x) eigenvalue problems.
The first two examples show the simple use case of finding a few of
the eigenvalues of largest magnitude:

  * [Compute the largest eigenpairs of an eigenvalue problem using block Davidson.](AnasaziBlockDavidson.md)
  * [Compute the largest eigenpairs of an eigenvalue problem using LOBPCG.](AnasaziLOBPCG.md)

The next two examples show how to use inverse iteration with Block
Krylov-Schur to find a few of the eigenvalues of smallest magnitude.
You may use just about any linear solver for inverse iteration.  The
following examples illustrate this for two different Trilinos linear
solvers.

  * [Inverse iteration using the KLU sparse direct solver through Amesos.](AnasaziBlockKrylovSchurEpetraExGenAmesos.md)
  * [Inverse iteration using an AztecOO iterative linear solver with an Ifpack preconditioner.](AnasaziBlockKrylovSchurEpetraExGenAztecOO.md)

### Zoltan tutorial ###

The Zoltan developers generously contributed a [hands-on tutorial](ZoltanHandsOnTutorial.md) of their own.

## Examples illustrating other packages ##

  1. **[PyTrilinos](http://trilinos.org/packages/pytrilinos) tutorial materials.**
    * [PyTrilinos tutorial slides](http://trilinos.org/oldsite/packages/pytrilinos/PyTrilinosTutorial.pdf)
    * [WebTrilinos for Python](https://www.users.csbsju.edu/trilinos/WebTrilinosMPI/python/index.html)
    * [Solve a simple 2D Laplace problem using PyTrilinos](PyTrilinosExample.md)
  1. **Learn how to use utilities in the [Teuchos](http://trilinos.sandia.gov/packages/teuchos) package.**
    * **Basic Support Tools**
      * [Build a parameter list (used to pass parameters to all Trilinos packages).](TeuchosPL.md)
      * [Build a reference-counted pointer (used to eliminate memory leak issues in most Trilinos packages).](TeuchosRCP.md)
      * [Build a command-line parser (tool for changing runtime behavior of program, providing documentation for options).](TeuchosCLP.md) [[unhighlighted](TeuchosCLP_unhighlighted.md)]
      * [Build a time monitor (tool for timing individual methods or operations in a program).](TeuchosTime.md)
    * **Template Support Tools**
      * [Create a templated BLAS wrapper.](TeuchosBLAS.md)
      * [Create a templated LAPACK wrapper.](TeuchosLAPACK.md)
      * [Create a templated serial, dense matrix.](TeuchosSDM.md)

## Other examples ##

Some examples don't work with the web tutorial, since they read from files.  You can try them out by downloading Trilinos and looking in the examples in the source tree.  For example, the [Intrepid](http://trilinos.org/packages/intrepid) discretizations package has examples in the `packages/trilinoscouplings/examples/scaling/` directory.  The following might be of interest:
  * example\_CurlLSFEM.cpp: driver for solving div-curl first order system in 3D with tangential boundary condition using curl-conforming elements
  * example\_DivLSFEM.cpp: driver for the same system with normal boundary condition and div-conforming elements
  * example\_Poisson.cpp: solving the Poisson equation using a Galerkin finite element method

## Learning more: ##

  * Web sites of Trilinos User Group meetings
    * [2013](http://trilinos.org/community/events/trilinos-user-group-2013/)
    * [2012](http://trilinos.org/community/events/trilinos-user-group-meeting-2012/)
    * [2011](http://trilinos.org/community/events/trilinos-user-group-meeting-2011/)
    * [2010](http://trilinos.org/community/events/trilinos-user-group-meeting-2010/)
    * [2009](http://trilinos.org/community/events/trilinos-user-group-meeting-2009/)
    * [2008](http://trilinos.org/community/events/trilinos-user-group-meeting-2008/)
    * [2007](http://trilinos.org/community/events/trilinos-user-group-meeting-2007/)
  * [Mail lists.](http://trilinos.org/community/mail-lists/)