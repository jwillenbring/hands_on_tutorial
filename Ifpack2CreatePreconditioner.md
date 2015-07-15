# Lesson topics #

This lesson shows how to create an Ifpack2 preconditioner from a Tpetra sparse matrix.  Ifpack2 is a package that can compute various kinds of preconditioners or smoothers, including the following:
  * Relaxations (Jacobi, Gauss-Seidel, and Symmetric Gauss-Seidel)
  * Chebyshev (a polynomial preconditioner)
  * ILUT (incomplete LU with threshold and fill level parameters)
  * RILUK (incomplete LU with threshold, fill level, and overlap parameters)
  * Diagonal (useful mainly only for comparison)

Ifpack2's preconditioners are defined on the calling process.  RILUK also allows overlap between process' domains.

Ifpack2 accepts a [Tpetra](http://trilinos.sandia.gov/packages/tpetra/) sparse matrix object (CrsMatrix) as input.  The preconditioners that it computes are Tpetra objects that fully exploit "MPI+X" parallelism (X = OpenMP, Pthreads, Intel's Threading Building Blocks, ...).  They may be used either as preconditioners in the iterative linear solvers package [Belos](http://trilinos.sandia.gov/packages/belos), as smoothers, or in your own iterative solvers or algebraic multigrid implementations.

# Related lessons #

This example will make more sense if you have done the Tpetra lessons.  In particular, you should understand the [Tpetra power method](TpetraPowerMethod.md) example, that creates a Tpetra sparse matrix and vectors and uses them for a simple numerical algorithm.

The [Belos iterative solver example](Tpetra_Belos_CreateSolver.md) shows how to use an Ifpack2 conditioner in an iterative solve.  It uses this lesson's example code to create an Ifpack2 preconditioner, so you can apply what you've learned about Ifpack2 to that lesson.

# Requirements #

Building this example requires using the 10.8 release or later of Trilinos, because it uses a new feature of [Teuchos::TimeMonitor](http://trilinos.sandia.gov/packages/docs/dev/packages/teuchos/doc/html/classTeuchos_1_1TimeMonitor.html).  If you would like to run this example with Trilinos 10.6, you will need to change the `getNewCounter` calls to `getNewTimer`.

# Exercises: Things to change #

## The matrix to factor ##

We include a very simple test linear system: a 1-D finite difference discretization of the Poisson equation with Dirichlet boundary conditions.  This means that the convergence rate should be relatively insensitive to preconditioner parameters.  You might like to change the test problem to make it harder to solve.  If you are running Trilinos directly rather than through the web interface, you might like to try reading in a file.  Otherwise, you might like to make the sparsity pattern of the matrix more complicated, so that an incomplete factorization with zero fill is not an exact LU factorization, or make the problem nonsymmetric, or make it more ill conditioned.  What happens if you make the matrix numerically or even structurally singular?

## The type of preconditioner ##

The example below computes an [ILUT preconditioner](http://www-users.cs.umn.edu/~saad/PDF/umsi-92-38.pdf).  ILUT is an incomplete LU factorization with an upper bound on the fill level, and two drop thresholds (relative and absolute).  It might be interesting to adjust parameters like the fill level.  Zero fill means that the structure of the L and U factors is exactly the same as the structure of their corresponding parts of the input matrix A.

You might also like to try a different preconditioner.  For example, the following implementation of `getPrecondTypeAndParameters()` (function in the full code example below) creates a Jacobi relaxation preconditioner.
```
void 
getPrecondTypeAndParameters (std::string& precondType, Teuchos::ParameterList& pl)
{
  using Teuchos::ParameterList;

  // The name of the type of preconditioner to use.  "RELAXATION"
  // includes Jacobi and variants of Gauss-Seidel.  You can pick which
  // one of these to use in the ParameterList.
  precondType = "RELAXATION";

  // The type of relaxation.  Options include "Jacobi",
  // "Gauss-Seidel", and "Symmetric Gauss-Seidel".
  const std::string relaxationType ("Jacobi");

  // The number of relaxation sweeps to do.
  const int numSweeps = 1;

  // Damping factor for Jacobi.  Default is 1.0.  
  //
  // NOTE: In general, 'double' should be replaced with a ScalarType
  // template parameter (the type of entries in the sparse matrix).
  // Complex-valued damping factors are legal if ScalarType is complex.
  const double dampingFactor = 1.0;

  // Diagonal elements with absolute value smaller than this are set to this.
  // Default is 0.0, meaning that zero diagonal elements are left alone.
  // Remember that relaxation methods divide by the diagonal in each row.
  //
  // NOTE: In general, 'double' should be replaced with
  // Teuchos::ScalarTraits<ScalarType>::magnitudeType.  Absolute values
  // (a.k.a. magnitudes) have "magnitude type," which is real if
  // ScalarType is complex.
  const double minDiagVal = 1.0e-8;

  pl.set ("relaxation: type", relaxationType);
  pl.set ("relaxation: damping factor", dampingFactor);
  pl.set ("relaxation: min diagonal value", minDiagVal);
}
```
Note that the "compute condition number estimate" feature of Ifpack2 doesn't work so well with a relaxation preconditioner.

## Note on performance ##

If you are running this example on the demo web server, please be aware that performance results are not necessarily indicative of performance "in the field."  This is because all of you are sharing the server, so the available memory bandwidth and CPU cycles per person is less than it would normally be if you were running on a dedicated compute node or on your own lightly loaded workstation.

# Code example #

```
//
// Example: Create an Ifpack2 preconditioner from a Tpetra::CrsMatrix.
//
#include <Ifpack2_Factory.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_Version.hpp>
#include <iostream>

// Get the Ifpack2 preconditioner type and its parameter list.
// You may modify this function to change the preconditioner type
// and its parameters.
//
// The first output argument is the preconditioner name.  In this
// case, it's "ILUT", for Saad's ILUT incomplete factorization
// preconditioner.
//
// The second output argument is the parameter list for the
// preconditioner.  Give it to the preconditioner's setParameters()
// method.  The parameter list this function returns tells ILUT to
// use fill level 2, drop tolerance 0, and absolute threshold 0.1.
//
// Note that with Ifpack2, the type of preconditioner is separate
// from the ParameterList for that preconditioner.
void 
getPrecondTypeAndParameters (std::string& precondType, Teuchos::ParameterList& pl)
{
  using Teuchos::ParameterList;

  // The name of the type of preconditioner to use.
  precondType = "ILUT";

  // Ifpack2 expects arguments of type 'double' here, regardless of
  // the scalar or magnitude types of the entries of the sparse
  // matrix.
  const double fillLevel = 2.0;
  const double dropTol = 0.0;
  const double absThreshold = 0.1;

  pl.set ("fact: ilut level-of-fill", fillLevel);
  pl.set ("fact: drop tolerance", dropTol);
  pl.set ("fact: absolute threshold", absThreshold);
}

// This function encapsulates creation of an Ifpack2 preconditioner
// from a Tpetra::CrsMatrix.  It returns the preconditioner as a
// Tpetra::Operator, which is the parent class of
// Ifpack2::Preconditioner.
//
// The template parameter TpetraMatrixType must be a specialization of
// Tpetra::CrsMatrix.  We template this function on the matrix type,
// rather than on the five template arguments of Tpetra::CrsMatrix,
// because CrsMatrix has some nice typedefs that let us retrieve those
// template arguments.  It's easier to template on one thing than on
// five things!  Recall that if T is a template parameter, and if T is
// a class with a typedef inside T::U, then you have to use "typename"
// to get at the typedef U.
//
// You don't have to use this function to make an Ifpack2
// preconditioner.  We just find it easier to read the code example if
// we wrap up the preconditioner creation in its own little function.
template<class TpetraMatrixType>
Teuchos::RCP<Tpetra::Operator<typename TpetraMatrixType::scalar_type,
			      typename TpetraMatrixType::local_ordinal_type,
			      typename TpetraMatrixType::global_ordinal_type,
			      typename TpetraMatrixType::node_type> >
createPreconditioner (const Teuchos::RCP<const TpetraMatrixType>& A,
		      const std::string& precondType,
		      const Teuchos::ParameterList& plist,
		      std::ostream& out,
		      std::ostream& err)
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::Time;
  using Teuchos::TimeMonitor;
  using std::endl;

  // Fetch the typedefs defined by Tpetra::CrsMatrix.
  typedef typename TpetraMatrixType::scalar_type scalar_type;
  typedef typename TpetraMatrixType::local_ordinal_type local_ordinal_type;
  typedef typename TpetraMatrixType::global_ordinal_type global_ordinal_type;
  typedef typename TpetraMatrixType::node_type node_type;

  // Ifpack2's generic Preconditioner interface implements
  // Tpetra::Operator.  A Tpetra::Operator is an abstraction of a
  // function mapping a (Multi)Vector to a (Multi)Vector, with the
  // option of applying the transpose or conjugate transpose of the
  // operator.  Tpetra::CrsMatrix implements Operator as well.
  typedef Tpetra::Operator<scalar_type, local_ordinal_type, 
                           global_ordinal_type, node_type> op_type;

  // These are just some convenience typedefs.
  typedef Teuchos::ScalarTraits<scalar_type> STS;
  typedef typename STS::magnitudeType magnitude_type;
  typedef Teuchos::ScalarTraits<magnitude_type> STM;

  // An Ifpack2::Preconditioner is-a Tpetra::Operator.  Ifpack2
  // creates a Preconditioner object, but users of iterative methods
  // want a Tpetra::Operator.  That's why create() returns an Operator
  // instead of a Preconditioner.
  typedef Ifpack2::Preconditioner<scalar_type, local_ordinal_type, 
                                  global_ordinal_type, node_type> prec_type;

  // Create timers to show how long it takes for Ifpack2 to do various operations.
  RCP<Time> initTimer = TimeMonitor::getNewCounter ("Ifpack2::Preconditioner::initialize");
  RCP<Time> computeTimer = TimeMonitor::getNewCounter ("Ifpack2::Preconditioner::compute");
  RCP<Time> condestTimer = TimeMonitor::getNewCounter ("Ifpack2::Preconditioner::condest");

  err << "Creating ILUT preconditioner" << endl 
      << "-- Configuring" << endl;
  //
  // Create the preconditioner and set parameters.
  //
  // This doesn't actually _compute_ the preconditioner.
  // It just sets up the specific type of preconditioner and
  // its associated parameters (which depend on the type).
  //
  RCP<prec_type> prec;
  Ifpack2::Factory factory;
  // Set up the preconditioner of the given type.
  prec = factory.create (precondType, A);
  prec->setParameters (plist);

  err << "-- Initializing" << endl;
  {
    TimeMonitor mon (*initTimer);
    prec->initialize();
  }

  // THIS ACTUALLY COMPUTES THE PRECONDITIONER
  // (e.g., does the incomplete factorization).
  err << "-- Computing" << endl;
  {
    TimeMonitor mon (*computeTimer);
    prec->compute();
  }

  if (precondType != "RELAXATION") {
    err << "-- Estimating condition number" << endl;
    magnitude_type condest = STM::one();
    {
      TimeMonitor mon (*condestTimer);
      condest = prec->computeCondEst (Ifpack2::Cheap);
    }
    out << endl << "Ifpack2 preconditioner's estimated condition number: " << condest << endl;
  }
  return prec;
}

// Create and return a simple example CrsMatrix.
template<class TpetraMatrixType>
Teuchos::RCP<const TpetraMatrixType>
createMatrix (const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
              const Teuchos::RCP<typename TpetraMatrixType::node_type>& node)
{
  using Teuchos::arcp;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::Time;
  using Teuchos::TimeMonitor;
  using Teuchos::tuple;

  typedef TpetraMatrixType matrix_type;

  // Fetch the timer for sparse matrix creation.
  //
  // If you are using Trilinos 10.6 instead of the development branch
  // (10.7), just create a new timer here, and remove the bit in
  // main() that creates a timer.
  //
  RCP<Time> timer = TimeMonitor::lookupCounter ("Sparse matrix creation");
  if (timer.is_null())
    timer = TimeMonitor::getNewCounter ("Sparse matrix creation");

  // Time the whole scope of this routine, not counting timer lookup.
  TimeMonitor monitor (*timer);

  // Fetch typedefs from the Tpetra::CrsMatrix.
  typedef typename TpetraMatrixType::scalar_type scalar_type;
  typedef typename TpetraMatrixType::local_ordinal_type local_ordinal_type;
  typedef typename TpetraMatrixType::global_ordinal_type global_ordinal_type;
  typedef typename TpetraMatrixType::node_type node_type;

  // The type of the Tpetra::Map that describes how the matrix is distributed.
  typedef Tpetra::Map<local_ordinal_type, global_ordinal_type, node_type> map_type;

  // The global number of rows in the matrix A to create.  We scale
  // this relative to the number of (MPI) processes, so that no matter
  // how many MPI processes you run, every process will have 10 rows.
  const Tpetra::global_size_t numGlobalElements = 10 * comm->getSize();
  
  // Construct a Map that puts approximately the same number of
  // equations on each processor.
  const global_ordinal_type indexBase = 0;
  RCP<const map_type > map = 
    rcp (new map_type (numGlobalElements, indexBase, comm, 
                       Tpetra::GloballyDistributed, node));

  // Get update list and the number of equations that this MPI process
  // owns.
  const size_t numMyElements = map->getNodeNumElements();
  ArrayView<const global_ordinal_type> myGlobalElements = map->getNodeElementList();

  // NumNz[i] will be the number of nonzero entries for the i-th
  // global equation on this MPI process.
  ArrayRCP<size_t> NumNz = arcp<size_t> (numMyElements);

  // We are building a tridiagonal matrix where each row is (-1 2 -1),
  // so we need 2 off-diagonal terms (except for the first and last
  // equation).
  for (size_t i = 0; i < numMyElements; ++i) {
    if (myGlobalElements[i] == 0 || static_cast<Tpetra::global_size_t>(myGlobalElements[i]) == numGlobalElements-1) {
      NumNz[i] = 2; // First or last equation
    } else {
      NumNz[i] = 3;
    }
  }

  // Create a Tpetra::Matrix using the Map, with a static allocation
  // dictated by NumNz.  (We know exactly how many elements there will
  // be in each row, so we use static profile for efficiency.)
  RCP<matrix_type> A = rcp (new matrix_type (map, NumNz, Tpetra::StaticProfile));
  
  // We are done with NumNZ; free it.
  NumNz = Teuchos::null;

  // Add rows one at a time.  Off diagonal values will always be -1.
  const scalar_type two    = static_cast<scalar_type> ( 2.0);
  const scalar_type negOne = static_cast<scalar_type> (-1.0);

  for (size_t i = 0; i < numMyElements; i++) {
    if (myGlobalElements[i] == 0) {
      A->insertGlobalValues (myGlobalElements[i], 
                             tuple (myGlobalElements[i], myGlobalElements[i]+1),
                             tuple (two, negOne));
    }
    else if (static_cast<Tpetra::global_size_t> (myGlobalElements[i]) == numGlobalElements-1) {
      A->insertGlobalValues (myGlobalElements[i],
                             tuple (myGlobalElements[i]-1, myGlobalElements[i]),
                             tuple (negOne, two));
    }
    else {
      A->insertGlobalValues (myGlobalElements[i],
                             tuple (myGlobalElements[i]-1, myGlobalElements[i], myGlobalElements[i]+1),
                             tuple (negOne, two, negOne));
    }
  }

  // Finish up the matrix.
  A->fillComplete ();
  return A;
}

int 
main (int argc, char *argv[]) 
{
  using std::endl;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ParameterList;
  using Teuchos::Time;
  using Teuchos::TimeMonitor;

  Teuchos::oblackholestream blackHole;
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &blackHole);
  RCP<const Teuchos::Comm<int> > comm = 
    Tpetra::DefaultPlatform::getDefaultPlatform().getComm();

  // Set up Tpetra typedefs.
  typedef double scalar_type;
  typedef int local_ordinal_type;
  typedef long global_ordinal_type;
  typedef Kokkos::DefaultNode::DefaultNodeType node_type;

  // Create the Kokkos Node instance.  
  // Tpetra objects use this object for intranode parallel operations.
  RCP<node_type> node = Kokkos::DefaultNode::getDefaultNode ();

  const int myRank = comm->getRank();
  const int numProcs = comm->getSize();
  std::ostream& out = (myRank == 0) ? std::cout : blackHole;
  std::ostream& err = (myRank == 0) ? std::cerr : blackHole;

  // Make a timer for sparse matrix creation.
  //
  // If you are using Trilinos 10.6 instead of the development branch
  // (10.7), just delete this line of code, and make the other change
  // mentioned above.
  RCP<Time> sparseMatrixCreationTimer = 
    TimeMonitor::getNewCounter ("Sparse matrix creation");

  // Run the whole example: create the sparse matrix, and compute the preconditioner.
  // Print out the Tpetra software version information.
  out << Tpetra::version() << endl << endl;

  typedef Tpetra::CrsMatrix<scalar_type, local_ordinal_type, global_ordinal_type, node_type> matrix_type;
  typedef Tpetra::Operator<scalar_type, local_ordinal_type, global_ordinal_type, node_type> op_type;

  // Create an example sparse matrix.
  RCP<const matrix_type> A = createMatrix<matrix_type> (comm, node);

  // Get the preconditioner type and its parameter list.  Modify the
  // definition of that function if you want to use a different
  // preconditioner or parameters.
  std::string precondType;
  ParameterList plist;
  getPrecondTypeAndParameters (precondType, plist);

  // Compute the preconditioner using the matrix A.
  // The matrix A itself is not modified.
  RCP<op_type> M = createPreconditioner<matrix_type> (A, precondType, plist, out, err);

  // Summarize global performance timing results, for all timers
  // created using TimeMonitor::getNewCounter().
  TimeMonitor::summarize (out);

  return 0;
}
```