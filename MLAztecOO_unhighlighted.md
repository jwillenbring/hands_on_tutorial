
```

// Two-level domain decomposition preconditioner with AztecOO and ML

#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_IntVector.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Time.h"
#include "AztecOO.h"

// includes required by ML
#include "ml_include.h"
#include "Epetra_LinearProblem.h"
#include "ml_MultiLevelOperator.h"
#include "ml_epetra_utils.h"

#include "Trilinos_Util_CommandLineParser.h"
#include "Trilinos_Util_CrsMatrixGallery.h"

using namespace Trilinos_Util;

// =========== //
// MAIN DRIVER //
// =========== //

int main(int argc, char *argv[])
{

#ifdef EPETRA_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else

  Epetra_SerialComm Comm;
#endif

  // initialize the command line parser
  CommandLineParser CLP(argc,argv);

  // initialize an Gallery object
  CrsMatrixGallery Gallery("", Comm);

  // add default values
  if( CLP.Has("-problem_type") == false ) CLP.Add("-problem_type", "laplace_2d" ); 
  if( CLP.Has("-problem_size") == false ) CLP.Add("-problem_size", "100" ); 

  // initialize the gallery as specified in the command line

  Gallery.Set(CLP);

  // retrive pointers to matrix and linear problem
  Epetra_CrsMatrix * Matrix = Gallery.GetMatrix();

  const Epetra_Map * Map = Gallery.GetMap();

  Epetra_LinearProblem * Problem = Gallery.GetLinearProblem();

  // Construct a solver object for this problem
  AztecOO solver(*Problem);

  // solve with CG (change is matrix is not symmetric) 
  solver.SetAztecOption(AZ_solver, AZ_cg);

  // ==============================  M L   S E C T I O N  ============================ //

  // Create and set an ML multilevel preconditioner
  ML *ml_handle;

  // Maximum number of levels
  int N_levels = 10;

  // output level
  ML_Set_PrintLevel(3);

  ML_Create(&ml_handle,N_levels);

  // wrap Epetra Matrix into ML matrix (data is NOT copied)
  EpetraMatrix2MLMatrix(ml_handle, 0, Matrix);

  // as we are interested in smoothed aggregation, create a ML_Aggregate object
  // to store the aggregates

  ML_Aggregate *agg_object;
  ML_Aggregate_Create(&agg_object);

  // specify max coarse size (ML will not coarse further is the matrix at a given level is
  // smaller than specified here)
  ML_Aggregate_Set_MaxCoarseSize(agg_object,1);

  // generate the hierady
  N_levels = ML_Gen_MGHierarchy_UsingAggregation(ml_handle, 0,
                                                  ML_INCREASING, agg_object);

  // Set a symmetric Gauss-Seidel smoother for the MG method (change
  // if the matrix is not symmetric)
  ML_Gen_Smoother_SymGaussSeidel(ml_handle, ML_ALL_LEVELS,
                                  ML_BOTH, 1, ML_DEFAULT);

  // generate solver
  ML_Gen_Solver    (ml_handle, ML_MGV, 0, N_levels-1);

  // wrap ML_Operator into Epetra_Operator
  ML_Epetra::MultiLevelOperator  MLop(ml_handle,Comm,*Map,*Map);

  // ===================== E N D    O F  M L   S E C T I O N  ========================= //

  // set this operator as preconditioner for AztecOO
  solver.SetPrecOperator(&MLop);

  // solve
  solver.Iterate(1550, 1e-12);

  // verify that residual is really small  
  double residual, diff;

  Gallery.ComputeResidual(&residual);

  Gallery.ComputeDiffBetweenStartingAndExactSolutions(&diff);

if( Comm.MyPID() == 0 ) {

    cout << "||b-Ax||_2 = " << residual << endl;
    cout << "||x_exact - x||_2 = " << diff << endl;
  }

  if (residual > 1e-5)
    exit(EXIT_FAILURE);

#ifdef EPETRA_MPI
  MPI_Finalize() ;
#endif
  exit(EXIT_SUCCESS);
}

```