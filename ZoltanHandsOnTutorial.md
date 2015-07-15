# Zoltan Hands-on Tutorial #

## First things first: ##

  * The Zoltan library is a collection of data management services for parallel, unstructured, adaptive, and dynamic applications. It simplifies the load-balancing, data movement, unstructured communication, and memory usage difficulties that arise in dynamic applications such as adaptive finite-element methods, particle methods, and crash simulations.
  * Zoltan can be downloaded as a standalone library from the http://www.cs.sandia.gov/~web1400/1400_download.html. Zoltan is also part of the Trilinos project. For information on the Trilinos project visit  [website](http://trilinos.sandia.gov) (click on logo above).
  * For future reference, a lot of information can be found on the [Zoltan](http://www.cs.sandia.gov/Zoltan/Zoltan.html) page.

## Let's dive in, you have two options: ##

  * **WebTrilinos**
    * [WebTrilinos](http://trilinos.sandia.gov/packages/webtrilinos) is a scientific portal, a web-based environment to use several Trilinos packages through the web.
    * **Note:** Only packages that Trilinos was configured with can be used through the WebTrilinos interfaces.
    * Installed at St. John's University, MN.  The following versions of the C++ interface are available (the examples below use the C++ interface):
      * [Trilinos development branch (10.7), MPI-enabled](https://www.users.csbsju.edu/trilinos/WebTrilinosMPI/c++/index.html)
      * [Trilinos 10.6, MPI-enabled](https://www.users.csbsju.edu/trilinos/WebTrilinosMPI-shared-10.6/c++/index.html)
      * [Trilinos 10.4, serial](https://www.users.csbsju.edu/trilinos/WebTrilinosSERIAL-10.4/c++/index.html)

<a href='Hidden comment: 
* [https://www.users.csbsju.edu/trilinos/WebTrilinosSERIAL-CUDA-10.6.1/c++/index.html Trilinos 10.7 (dev), CUDA-enabled]
'></a>
  * In addition to the C++ interfaces, a [Matrix Portal interface](https://www.users.csbsju.edu/trilinos/WebTrilinosMPI-shared-10.4/MatrixPortal/index.html) and a [Python interface](https://www.users.csbsju.edu/trilinos/WebTrilinosMPI-shared-10.4/python/index.html) are also available.
  * Access to this site is password protected.  Login information will be given during live tutorials as needed.
  * **Reminder:** Use Ctrl+A to highlight all the example code, Ctrl+C to copy it, and Ctrl+V to paste it in the WebTrilinos window.
  * **Download and build the examples on your machine**
    1. Download [Trilinos](http://trilinos.sandia.gov/download)
    1. Download and install CMake (if not already done).
      * The latest CMake release can be downloaded from [here](http://www.cmake.org/cmake/resources/software.html).
    1. Download and install Clapack (if you don't have LAPACK and BLAS already).
      * Get Clapack from [here](http://www.netlib.org/clapack/clapack-3.2.1-CMAKE.tgz).
      * Use CMake to build Clapack.
    1. Build and install Trilinos.
      * An example script is [here](BuildScript.md).
      * You can find other scripts in the Trilinos distribution Trilinos/sample\_scripts directory.
      * You can also use the CMake gui, or the text-gui ccmake.
    1. Get Makefile from [here](http://trilinos.sandia.gov/Export_Makefile_example.txt).
      * Learn about the Makefile and the Makefile.export system [here](http://trilinos.sandia.gov/Export_Makefile.txt).
      * You can find the Makefile.export.package\_name files once you have built and installed Trilinos.
      * They will be in the include directory of your installation directory.
      * Example: `TrilinosInstall/include/Makefile.export.Epetra.`
    1. Customize the Makefile to your situation (the example Makefile uses only the Epetra package, but Makefile.export.package\_name files are available for all packages).
    1. Add Example code into a file called main.cpp and build away!

## Now for some examples: ##

  1. **Load balancing using Zoltan**
    * [Example that uses Recursive Coordinate Bisection in Zoltan to partition a mesh input](ZoltanSimpleRCB.md)
    * [Example that uses Graph Partitioning in Zoltan to partition a graph input](ZoltanSimpleGraph.md)
    * [Example that uses Hypergraph Partitooning in Zoltan to partition a hypergraph input](ZoltanSimplePHG.md)
    * [Example that uses Graph Repartioning in Zoltan to repartition a graph input](ZoltanGraphRepartition.md)
    * [Example that uses Graph Partitioning in Zoltan to partition a graph input and migrate the data](ZoltanGraphMigrate.md)

  1. **Matrix Partitioning using Isorropia**
    * [Example that partitions an Epetra matrix in Trilinos and redistributes the matrix](IsorropiaMatrixPartition.md)


## Learning more: ##

  * [Zoltan website](http://www.cs.sandia.gov/Zoltan/Zoltan.html).
  * [Zoltan User Guide](http://www.cs.sandia.gov/Zoltan/ug_html/ug.html).
  * [Zoltan Related Papers and Presentations](http://www.cs.sandia.gov/Zoltan/Zoltan_pubs.html)
  * [Isorropia website](http://trilinos.sandia.gov/packages/docs/dev/packages/isorropia/doc/html/index.html)