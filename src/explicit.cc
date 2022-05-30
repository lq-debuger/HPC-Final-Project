static char help[] = "Solves a tridiagonal linear system.\n\n";

#include <stdio.h>
#include <petscksp.h>

int main(int argc,char **args)
{
  Vec            y, z;          
  Mat            A;
  KSP            ksp;
  PC             pc;                  
  PetscReal      norm,norm_old=1.0,tol=1000.*PETSC_MACHINE_EPSILON;  
  PetscErrorCode ierr;
  PetscInt       i,n = 10,col[3],rstart,rend,nlocal,maxit=10000;
  PetscScalar    value[3],value_vec,lambda;

  ierr = PetscInitialize(&argc,&args,(char*)0,help);if (ierr) return ierr;
  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);

  ierr = PetscInitialize(&argc,&args,(char*)0,help);if (ierr) return ierr;
  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);

  ierr = VecCreate(PETSC_COMM_WORLD,&z);CHKERRQ(ierr);
  ierr = VecSetSizes(z,PETSC_DECIDE,n);CHKERRQ(ierr);
  ierr = VecSetFromOptions(z);CHKERRQ(ierr);
  ierr = VecDuplicate(z,&y);CHKERRQ(ierr);

  ierr = VecGetOwnershipRange(z,&rstart,&rend);CHKERRQ(ierr);
  ierr = VecGetLocalSize(z,&nlocal);CHKERRQ(ierr);


  ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
  ierr = MatSetSizes(A,nlocal,nlocal,n,n);CHKERRQ(ierr);
  ierr = MatSetFromOptions(A);CHKERRQ(ierr);
  ierr = MatSetUp(A);CHKERRQ(ierr);

  // ierr = MatCreate(PETSC_COMM_WORLD,&B);CHKERRQ(ierr);
  // ierr = MatSetSizesBA,PETSC_DECIDE,nlocal,1,n);CHKERRQ(ierr);
  // ierr = MatSetFromOptions(B);CHKERRQ(ierr);
  // ierr = MatSetUp(B);CHKERRQ(ierr);

  if (!rstart) 
  {
    rstart = 1;
    i      = 0; col[0] = 0; col[1] = 1; value[0] = 2.0; value[1] = -1.0; value_vec=1.0;
    ierr   = MatSetValues(A,1,&i,2,col,value,INSERT_VALUES);CHKERRQ(ierr);
    ierr   = VecSetValues(z,1,&i,&value_vec,INSERT_VALUES);CHKERRQ(ierr);
  }
  
  if (rend == n) 
  {
    rend = n-1;
    i    = n-1; col[0] = n-2; col[1] = n-1; value[0] = -1.0; value[1] = 2.0;
    ierr = MatSetValues(A,1,&i,2,col,value,INSERT_VALUES);CHKERRQ(ierr);
  }

  /* Set entries corresponding to the mesh interior */
  value[0] = -1.0; value[1] = 2.0; value[2] = -1.0;
  for (i=rstart; i<rend; i++) 
  {
    col[0] = i-1; col[1] = i; col[2] = i+1;value_vec=0.0;
    ierr   = MatSetValues(A,1,&i,3,col,value,INSERT_VALUES);CHKERRQ(ierr);
    ierr   = VecSetValues(z,1,&i,&value_vec,INSERT_VALUES);CHKERRQ(ierr);
  }

  /* Assemble the matrix and vec */
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(z);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(z);CHKERRQ(ierr);

  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);

  ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);

 
  ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
  ierr = PCSetType(pc,PCJACOBI);CHKERRQ(ierr);
  ierr = KSPSetTolerances(ksp,1.e-7,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);

  ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
  ierr = PCSetFromOptions(pc);CHKERRQ(ierr);

  
  for ( i = 0; i < maxit; i++)
  {
    ierr = KSPSolve(ksp,z,y);CHKERRQ(ierr);
    ierr = VecNorm(y,NORM_2,&norm);CHKERRQ(ierr);
    ierr = VecScale(y,(PetscScalar)1.0/norm);CHKERRQ(ierr);
    ierr = VecCopy(y,z);CHKERRQ(ierr);
    if(PetscAbsScalar(norm-norm_old) <= tol){
      ierr = PetscPrintf(PETSC_COMM_WORLD,"yk is so close to yk-1\n");CHKERRQ(ierr);
      break;
    }
    norm_old = norm; 
  }

  ierr = PetscPrintf(PETSC_COMM_WORLD,"The eigenvector is :\n");CHKERRQ(ierr);
  ierr = VecView(z, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

  /*caculate eigenvalue */
   ierr = MatMult(A,z,y);CHKERRQ(ierr);
   ierr = VecDot(z,y,&lambda);CHKERRQ(ierr);
   ierr = PetscPrintf(PETSC_COMM_WORLD,"The eigenvector is : %g\n",(double)lambda);CHKERRQ(ierr);

  ierr = VecDestroy(&z);CHKERRQ(ierr); ierr = VecDestroy(&y);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = KSPDestroy(&ksp);CHKERRQ(ierr);

  ierr = PetscFinalize();
  return ierr;
}