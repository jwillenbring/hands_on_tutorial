
```

#include "Amesos_ConfigDefs.h"

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_CrsMatrix.h"
#include "Epetra_MultiVector.h"
#include "Epetra_LinearProblem.h"
#include "Galeri_Maps.h"
#include "Galeri_CrsMatrices.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Amesos_Superlu.h"

//=============================================================================
bool CheckError(const Epetra_RowMatrix& A,
		const Epetra_MultiVector& x,
		const Epetra_MultiVector& b,
		const Epetra_MultiVector& x_exact)
{
  std::vector<double> Norm;
  int NumVectors = x.NumVectors();
  Norm.resize(NumVectors);
  Epetra_MultiVector Ax(x);
  A.Multiply(false,x,Ax);
  Ax.Update(1.0,b,-1.0);
  Ax.Norm2(&Norm[0]);
  bool TestPassed = false;
  double TotalNorm = 0.0;
  for (int i = 0 ; i < NumVectors ; ++i) {
    TotalNorm += Norm[i];
  }
  if (A.Comm().MyPID() == 0)
    std::cout << "||Ax - b||  = " << TotalNorm << std::endl;
  return(TestPassed);
}

int main(int argc, char *argv[])
{

#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif

  Teuchos::ParameterList GaleriList;

  // The problem is defined on a 2D grid, global size is nx * nx.
  int nx = 30; 
  GaleriList.set("n", nx * nx);
  GaleriList.set("nx", nx);
  GaleriList.set("ny", nx);
  Teuchos::RCP<Epetra_Map> Map = Teuchos::rcp( Galeri::CreateMap("Linear", Comm, GaleriList) );
  Teuchos::RCP<Epetra_RowMatrix> A = Teuchos::rcp( Galeri::CreateCrsMatrix("Laplace2D", &*Map, GaleriList) );


  // =================================================== //
  // E N D   O F   I F P A C K   C O N S T R U C T I O N //
  // =================================================== //

  // At this point, we need some additional objects
  // to define and solve the linear system.

  // defines LHS and RHS
  Epetra_Vector LHS(A->OperatorDomainMap());
  Epetra_Vector EXACT(A->OperatorDomainMap());
  Epetra_Vector RHS(A->OperatorDomainMap());

  // solution is constant
  EXACT.PutScalar(1.0);
  // now build corresponding RHS
  A->Apply(EXACT,RHS);

  // now randomize the solution
  RHS.Random();

  // need an Epetra_LinearProblem to define AztecOO solver
  Epetra_LinearProblem Problem(&*A,&LHS,&RHS);

  Amesos_Superlu Solver(Problem);

  AMESOS_CHK_ERR(Solver.SymbolicFactorization());
  AMESOS_CHK_ERR(Solver.NumericFactorization());
  AMESOS_CHK_ERR(Solver.Solve());

CheckError(*A, LHS, RHS, EXACT);


#ifdef HAVE_MPI
  MPI_Finalize() ; 
#endif

  return(EXIT_SUCCESS);
}

```