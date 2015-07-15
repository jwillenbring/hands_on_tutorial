# Lesson topics #

The Epetra package provides distributed sparse linear algebra.  It
includes sparse matrices, vectors, and other linear algebra objects,
along with computational kernels. This lesson shows the MPI (or
non-MPI) initialization you need to do in order to start using Epetra.
The initialization procedure differs slightly, depending on whether
you are writing a code from scratch, or introducing Epetra into an
existing code base.  We will give example codes and discussion for the
following three use cases:

  1. A code which only uses MPI through Trilinos
> 2. A code which uses MPI on its own as well as through Trilinos
> 3. A code which does not use MPI

# Initialization for a code that only uses MPI through Trilinos #

This section explains how to set up the distributed-memory parallel
environment for using Epetra, in a code which only uses MPI through
Trilinos.  If you want to introduce Epetra into an existing MPI
application, please see the next section.  This example works whether
or not Trilinos was built with MPI support.

Epetra was written for distributed-memory parallel programming.  It
uses [MPI](http://en.wikipedia.org/wiki/Message_Passing_Interface)
(the Message Passing Interface) for this.  However, Epetra will work
correctly whether or not you have built Trilinos with MPI support.  It
does so by interacting with MPI through an interface called
Epetra\_Comm.  If MPI is enabled, then this wraps an MPI\_Comm.
Otherwise, this is a "serial communicator" with one process, analogous
to MPI\_COMM\_SELF.

Epetra expects that the user call MPI\_Init before using MPI, and call
MPI\_Finalize after using MPI (usually at the end of the program).  You
may either do this manually, or use Teuchos::GlobalMPISession.  The
latter calls MPI\_Init and MPI\_Finalize for you in an MPI build, and
does not call them if you did not build Trilinos with MPI support.
However, you may only use Teuchos::GlobalMPISession if Trilinos was
built with the Teuchos package enabled.  Epetra does not require the
Teuchos package, so the the following example illustrates the standard
idiom for initializing MPI (if available) and getting an Epetra\_Comm
corresponding to MPI\_COMM\_WORLD.  The example works whether or not
Trilinos was build with MPI support.

```
//
// This example includes conditional MPI initialization, getting an
// Epetra communicator wrapper, and printing out Epetra version
// information.
//

// This defines useful macros like HAVE_MPI, which is defined if and
// only if Epetra was built with MPI enabled.
#include <Epetra_config.h>

#ifdef HAVE_MPI
// Your code is an existing MPI code, so it presumably includes mpi.h directly.
#  include <mpi.h>
// Epetra's wrapper for MPI_Comm.  This header file only exists if
// Epetra was built with MPI enabled.
#  include <Epetra_MpiComm.h>
#else
#  include <Epetra_SerialComm.h>
#endif // HAVE_MPI

#include <Epetra_Version.h>


//
// ... Your other include files go here ...
//

// Do something with the given communicator.  In this case, we just
// print Epetra's version to the given output stream, on Process 0.
void
exampleRoutine (const Epetra_Comm& comm,
                std::ostream& out)
{
  if (comm.MyPID () == 0) {
    // On (MPI) Process 0, print out the Epetra software version.
    out << Epetra_Version () << std::endl << std::endl;
  }
}

int
main (int argc, char *argv[])
{
  // These "using" declarations make the code more concise, in that
  // you don't have to write the namespace along with the class or
  // object name.  This is especially helpful with commonly used
  // things like std::endl.
  using std::cout;
  using std::endl;

#ifdef HAVE_MPI
  // Start up MPI, if using MPI.  Trilinos doesn't have to be built
  // with MPI; it's called a "serial" build if you build without MPI.
  //
  // It's bad form to ignore the error codes returned by MPI
  // functions, but we do so here for brevity.
  (void) MPI_Init (&argc, &argv);

  // Wrap MPI_COMM_WORLD in an Epetra communicator wrapper.
  // Epetra_MpiComm is a subclass of Epetra_Comm, so you may use it
  // wherever an Epetra_Comm is required.
  Epetra_MpiComm comm (MPI_COMM_WORLD);
#else
  // Make a "serial" (non-MPI) communicator.  It doesn't actually
  // "communicate," because it only has one process, whose rank is
  // always 0.  Epetra_SerialComm is a subclass of Epetra_Comm, so you
  // may use it wherever an Epetra_Comm is required.
  Epetra_SerialComm comm;
#endif

  // Epetra_Comm has methods that wrap basic MPI functionality.
  // MyPID() is equivalent to MPI_Comm_rank, and NumProc() to
  // MPI_Comm_size.
  //
  // With a "serial" communicator, the rank is always 0, and the
  // number of processes is always 1.
  const int myRank = comm.MyPID ();
  const int numProcs = comm.NumProc ();

  if (myRank == 0) {
    cout << "Total number of processes: " << numProcs << endl;
  }

  // Do something with the new Epetra communicator.
  exampleRoutine (comm, cout);

  // This tells the Trilinos test framework that the test passed.
  if (comm.MyPID () == 0) {
    cout << "End Result: TEST PASSED" << endl;
  }

#ifdef HAVE_MPI
  // Since you called MPI_Init, you are responsible for calling
  // MPI_Finalize after you are done using MPI.
  (void) MPI_Finalize ();
#endif // HAVE_MPI

  return 0;
}
```

# Initialization for an existing MPI code #

Epetra also works fine in an existing MPI code.  For this example, we
assume that your code initializes MPI on its own by calling MPI\_Init,
and calls MPI\_Finalize at the end.  It also must get an MPI\_Comm (an
MPI communicator) somewhere, either by using a predefined communicator
such as MPI\_COMM\_WORLD, or by creating a new one.

```
//
// This example shows how to wrap the MPI_Comm (MPI communicator) that
// you are using, so that Epetra can use it as well.  it includes MPI
// initialization, wrapping your MPI_Comm in an Epetra communicator
// wrapper, and printing out Epetra version information.
//

// This defines useful macros like HAVE_MPI, which is defined if and
// only if Epetra was built with MPI enabled.
#include <Epetra_config.h>

#ifdef HAVE_MPI
// Your code is an existing MPI code, so it presumably includes mpi.h directly.
#  include <mpi.h>
// Epetra's wrapper for MPI_Comm.  This header file only exists if
// Epetra was built with MPI enabled.
#  include <Epetra_MpiComm.h>
#else
#  error "This example requires MPI in order to build."
#endif // HAVE_MPI

#include <Epetra_Version.h>

//
// ... Your other include files go here ...
//


// Do something with the given communicator.  In this case, we just
// print Epetra's version to the given output stream, on Process 0.
void
exampleRoutine (const Epetra_Comm& comm,
                std::ostream& out)
{
  if (comm.MyPID () == 0) {
    // On (MPI) Process 0, print out the Epetra software version.
    out << Epetra_Version () << std::endl << std::endl;
  }
}

int
main (int argc, char *argv[])
{
  // These "using" declarations make the code more concise, in that
  // you don't have to write the namespace along with the class or
  // object name.  This is especially helpful with commonly used
  // things like std::endl.
  using std::cout;
  using std::endl;

  // We assume that your code calls MPI_Init.  It's bad form
  // to ignore the error codes returned by MPI functions, but
  // we do so here for brevity.
  (void) MPI_Init (&argc, &argv);

  // This code takes the place of whatever you do to get an MPI_Comm.
  MPI_Comm yourComm = MPI_COMM_WORLD;

  // If your code plans to use MPI on its own, as well as through
  // Trilinos, you should strongly consider giving Trilinos a copy
  // of your MPI_Comm (created via MPI_Comm_dup).  Trilinos may in
  // the future duplicate the MPI_Comm automatically, but it does
  // not currently do this.

  // Wrap the MPI_Comm.  You are responsible for calling MPI_Comm_free
  // on your MPI_Comm after use, if necessary.  (It's not necessary or
  // legal to do this for built-in communicators like MPI_COMM_WORLD
  // or MPI_COMM_SELF.)
  Epetra_MpiComm comm (yourComm);

  // Epetra_Comm has methods that wrap basic MPI functionality.
  // MyPID() is equivalent to MPI_Comm_rank; it returns my process'
  // rank.  NumProc() is equivalent to MPI_Comm_size; it returns the
  // total number of processes in the communicator.
  const int myRank = comm.MyPID ();
  const int numProcs = comm.NumProc ();

  if (myRank == 0) {
    cout << "Total number of processes: " << numProcs << endl;
  }

  // Do something with the new Epetra communicator.
  exampleRoutine (comm, cout);

  // This tells the Trilinos test framework that the test passed.
  if (myRank == 0) {
    cout << "End Result: TEST PASSED" << endl;
  }

  // If you need to call MPI_Comm_free on your MPI_Comm, now would be
  // the time to do so, before calling MPI_Finalize.

  // Since you called MPI_Init, you are responsible for calling
  // MPI_Finalize after you are done using MPI.
  (void) MPI_Finalize ();
  return 0;
}
```

# Initialization for an existing non-MPI code #

If are using a build of Trilinos that has MPI enabled, but you don't
want to use MPI in your application, you may either imitate the first
example above, or create an Epetra\_SerialComm directly as the
"communicator."  The following example shows how to create an
Epetra\_SerialComm.

```
#include <Epetra_config.h>
// Wrapper for a "communicator" containing only one process.  This
// header file always exists, whether or not Epetra was built with MPI
// enabled.
#include <Epetra_SerialComm.h>
#include <Epetra_Version.h>

#include <cstdlib>
#include <sstream>
#include <stdexcept>

//
// ... Your other include files go here ...
//


// Do something with the given communicator.  In this case, we just
// print Epetra's version to the given output stream, on Process 0.
void
exampleRoutine (const Epetra_Comm& comm,
                std::ostream& out)
{
  if (comm.MyPID () == 0) {
    // On (MPI) Process 0, print out the Epetra software version.
    out << Epetra_Version () << std::endl << std::endl;
  }
}

int
main (int argc, char *argv[])
{
  // These "using" declarations make the code more concise, in that
  // you don't have to write the namespace along with the class or
  // object name.  This is especially helpful with commonly used
  // things like std::endl.
  using std::cout;
  using std::endl;

  // Make a "serial" (non-MPI) communicator.  It doesn't actually
  // "communicate," because it only has one process, whose rank is
  // always 0.  Epetra_SerialComm is a subclass of Epetra_Comm, so you
  // may use it wherever an Epetra_Comm is required.
  Epetra_SerialComm comm;

  // Epetra_Comm has methods that wrap basic MPI functionality.
  // MyPID() is equivalent to MPI_Comm_rank, and NumProc() to
  // MPI_Comm_size.
  //
  // With a "serial" communicator, the rank is always 0, and the
  // number of processes is always 1.
  const int myRank = comm.MyPID ();
  const int numProcs = comm.NumProc ();

  // Test the two assertions in the previous comment.
  if (numProcs != 1) {
    std::ostringstream err;
    err << "This is a serial (non-MPI) example, but the number of processes "
        << "in the Epetra_Comm is " << numProcs << " != 1.  Please report "
        << "this bug.";
    throw std::logic_error (err.str ());
  }
  if (myRank != 0) {
    std::ostringstream err;
    err << "This is a serial (non-MPI) example, but the rank of the calling "
      "process in the Epetra_Comm is " << myRank << " != 0.  Please report "
      "this bug.";
    throw std::logic_error (err.str ());
  }

  // Do something with the new communicator.
  exampleRoutine (comm, cout);

  // This tells the Trilinos test framework that the test passed.
  if (myRank == 0) {
    cout << "End Result: TEST PASSED" << endl;
  }

  return 0;
}
```

# Things we didn't explain above #

## Epetra\_Comm, Epetra\_MpiComm, and Epetra\_SerialComm ##

Epetra\_Comm is Epetra's interface to distributed-memory parallel
communication.  It is an abstract base class.  The Epetra\_MpiComm and
Epetra\_SerialComm classes implement this interface.  As the name
indicates, Epetra\_MpiComm implements Epetra\_Comm by using MPI calls.
Epetra\_SerialComm implements Epetra\_Comm without MPI, as a
"communicator" with only one process, whose rank is always zero.
(This is more or less equivalent to MPI\_COMM\_SELF, except without
actually using MPI.)

Since Epetra\_Comm is an abstract base class, you cannot create it
directly.  You must handle it by pointer or reference.  However, you
may create an instance of a subclass of Epetra\_Comm.  The above
examples show how to do this.