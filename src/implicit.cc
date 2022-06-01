static char help[] = "Solves the transient heat equation with implicit scheme.\n\n";

#include <stdio.h>
#include <petscksp.h>
#include <math.h>

int main(int argc,char **args)
{
  Vec            u_old, u;          
  Mat            A;
  KSP            ksp;
  PC             pc;                  
  PetscReal      kappa,dt,dx,rho,c,alpha,beta;  
  PetscErrorCode ierr;
  PetscInt       i,n = 256,col[3],rstart,rend,rank,nlocal,steps;
  PetscScalar    value[3],value_vec,xi;

  kappa= 1.0;
  dt   = 0.0002;
  steps= 20.0/dt;
  dx   = 1.0/n;
  rho  = 1.0;
  c    = 1.0;
  alpha= -(kappa*dt/rho/c/dx/dx);
  beta = 1-2*alpha;
  i    = 0;

  ierr = PetscInitialize(&argc,&args,(char*)0,help);if (ierr) return ierr;
  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);

  ierr = PetscInitialize(&argc,&args,(char*)0,help);if (ierr) return ierr;
  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);

  ierr = VecCreate(PETSC_COMM_WORLD,&u_old);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);CHKERRQ(ierr);
  ierr = VecSetSizes(u_old,PETSC_DECIDE,n);CHKERRQ(ierr);
  ierr = VecSetFromOptions(u_old);CHKERRQ(ierr);
  ierr = VecDuplicate(u_old,&u);CHKERRQ(ierr);

  ierr = VecGetOwnershipRange(u_old,&rstart,&rend);CHKERRQ(ierr);
  ierr = VecGetLocalSize(u_old,&nlocal);CHKERRQ(ierr);

  ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
  ierr = MatSetSizes(A,nlocal,nlocal,n,n);CHKERRQ(ierr);
  ierr = MatSetFromOptions(A);CHKERRQ(ierr);
  ierr = MatSetUp(A);CHKERRQ(ierr);

  /*set value for mat A*/
  if (!rstart) 
  {
    rstart = 1;
    i      = 0; col[0] = 0; col[1] = 1; value[0] = beta; value[1] = alpha;
    ierr   = MatSetValues(A,1,&i,2,col,value,INSERT_VALUES);CHKERRQ(ierr);
  }
  
  if (rend == n) 
  {
    rend = n-1;
    i    = n-1; col[0] = n-2; col[1] = n-1; value[0] = alpha; value[1] = beta;
    ierr = MatSetValues(A,1,&i,2,col,value,INSERT_VALUES);CHKERRQ(ierr);
  }

  /* Set entries corresponding to the mesh interior */
  value[0] = alpha; value[1] = beta; value[2] = alpha;
  for (i=rstart; i<rend; i++) 
  {
    col[0] = i-1; col[1] = i; col[2] = i+1;
    ierr   = MatSetValues(A,1,&i,3,col,value,INSERT_VALUES);CHKERRQ(ierr);
  }

  /* Assemble the matrix and vec */
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  // ierr = MatView(A, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

  /*set value for vec u_old*/
  if(rank == 0)
  {
    value_vec = 0;
    ierr = VecSetValues(u_old, 1, &i, &value_vec, INSERT_VALUES);CHKERRQ(ierr);
    for(i = 1; i < n; i++)
    {
      xi = i*dx;
      value_vec = exp(xi) + dt*sin(xi*M_PI)*dt/rho/c;
	    ierr = VecSetValues(u_old, 1, &i, &value_vec, INSERT_VALUES);CHKERRQ(ierr);
    }
    value_vec = 0;
    ierr = VecSetValues(u_old, 1, &i, &value_vec, INSERT_VALUES);CHKERRQ(ierr);
  }

  ierr = VecAssemblyBegin(u_old);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(u_old);CHKERRQ(ierr);
  ierr = VecView(u_old, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);

  ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);

 
  ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
  ierr = PCSetType(pc,PCJACOBI);CHKERRQ(ierr);
  ierr = KSPSetTolerances(ksp,1.e-7,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);

  ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
  ierr = PCSetFromOptions(pc);CHKERRQ(ierr);

  
  for ( i = 0; i < steps; i++)
  {
    ierr = KSPSolve(ksp,u_old,u);CHKERRQ(ierr);
    ierr = VecCopy(u_old,u);CHKERRQ(ierr); 
  }
  ierr = PetscPrintf(PETSC_COMM_WORLD,"---------------------\n");CHKERRQ(ierr);
  ierr = VecView(u_old, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  // ierr = PetscPrintf(PETSC_COMM_WORLD,"The eigenvector is :\n");CHKERRQ(ierr);
  // ierr = VecView(z, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD,"The program is finished\n");CHKERRQ(ierr);

  ierr = VecDestroy(&u_old);CHKERRQ(ierr); ierr = VecDestroy(&u);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = KSPDestroy(&ksp);CHKERRQ(ierr);

  ierr = PetscFinalize();
  return ierr;
}