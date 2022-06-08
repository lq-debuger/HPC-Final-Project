static char help[] = "Solves the transient heat equation with explicit and implicit scheme.\n\n";

#include <stdio.h>
#include <petscksp.h>
#include <petscviewerhdf5.h>
#include <math.h>

int main(int argc,char **args)
{
  Vec            u_old, u,f,u_exact;          
  Mat            A;
  #ifdef IMPLICIT  
  KSP            ksp;
  PC             pc; 
  #endif                 
  PetscReal      kappa,dt,dx,rho,c,alpha,beta,pi,norm=1.0,tol=1000.*PETSC_MACHINE_EPSILON;;  
  PetscErrorCode ierr;
  PetscInt       i,n = 128,col[3],rstart,rend,rank,nlocal,steps,step;
  PetscScalar    value[3],value_vec=0.0,value_f=0.0,vec_ue=0.0,xi;
  PetscViewer    viewer;
  PetscBool      restart = PETSC_FALSE;

  kappa= 1.0;
  // dt   = 0.00002;
  // dt   = 0.000015;
  // dt   = 0.000017;
  // dt   = 0.000019;
  dt   = 0.00006;
  steps= 100000;
  // steps= 1000000;
  rho  = 1.0;
  c    = 1.0;
  pi   = 4.0*atan(1);

  ierr = PetscInitialize(&argc,&args,(char*)0,help);if (ierr) return ierr;
  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL,NULL,"-restart",&restart,NULL);CHKERRQ(ierr);

  ierr = VecCreate(PETSC_COMM_WORLD,&u_old);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);CHKERRQ(ierr);
  ierr = VecSetSizes(u_old,PETSC_DECIDE,n);CHKERRQ(ierr);
  ierr = VecSetFromOptions(u_old);CHKERRQ(ierr);
  ierr = VecDuplicate(u_old,&u);CHKERRQ(ierr);
  ierr = VecDuplicate(u_old,&f);CHKERRQ(ierr);
  ierr = VecDuplicate(u_old,&u_exact);CHKERRQ(ierr);

  ierr = VecGetOwnershipRange(u_old,&rstart,&rend);CHKERRQ(ierr);
  ierr = VecGetLocalSize(u_old,&nlocal);CHKERRQ(ierr);

  ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
  ierr = MatSetSizes(A,nlocal,nlocal,n,n);CHKERRQ(ierr);
  ierr = MatSetFromOptions(A);CHKERRQ(ierr);
  ierr = MatSetUp(A);CHKERRQ(ierr);

  dx   = 1.0/n;
  step = 1;
  #ifdef IMPLICIT
  alpha= -(kappa*dt/rho/c/dx/dx);
  #elif EXPLICIT
  alpha= kappa*dt/rho/c/dx/dx;
  #endif
  beta = 1-2*alpha;
  i    = 0;
  ierr = PetscPrintf(PETSC_COMM_WORLD,"n=%d,dt=%f,dx=%f,steps=%d,alpha=%f,beta=%f\n",n,dt,dx,steps,alpha,beta);CHKERRQ(ierr);

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

  value[0] = alpha; value[1] = beta; value[2] = alpha;
  for (i=rstart; i<rend; i++) 
  {
    col[0] = i-1; col[1] = i; col[2] = i+1;
    ierr   = MatSetValues(A,1,&i,3,col,value,INSERT_VALUES);CHKERRQ(ierr);
  }

  /* Assemble the matrix A*/
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  // ierr = MatView(A, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

  /*set value for vec u_old */
  if(restart)
  {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Reading u_old from h5 ...\n");CHKERRQ(ierr);
    #ifdef IMPLICIT
    ierr = PetscViewerHDF5Open(PETSC_COMM_WORLD,"./h5/implicit.h5",FILE_MODE_READ,&viewer);CHKERRQ(ierr);
    #endif
    #ifdef EXPLICIT
    ierr = PetscViewerHDF5Open(PETSC_COMM_WORLD,"./h5/explicit.h5",FILE_MODE_READ,&viewer);CHKERRQ(ierr);
    #endif
    ierr = PetscObjectSetName((PetscObject) u_old, "restart_Vec");CHKERRQ(ierr);
    ierr = VecLoad(u_old,viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  }else
  {
    i    = 0;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Do not read u_old from uold.h5 ...\n");CHKERRQ(ierr);
    if(rank == 0)
    {
      value_vec = 0;
      ierr      = VecSetValues(u_old, 1, &i, &value_vec, INSERT_VALUES);CHKERRQ(ierr);
      for(i = 1; i < n-1; i++)
      {
        xi        = i*dx;
        // value_vec = exp(xi);
        value_vec = sin(xi*pi);
        ierr      = VecSetValues(u_old, 1, &i, &value_vec, INSERT_VALUES);CHKERRQ(ierr);
      }
      value_vec   = 0;
      ierr        = VecSetValues(u_old, 1, &i, &value_vec, INSERT_VALUES);CHKERRQ(ierr);
    }
  }
  
  ierr = VecAssemblyBegin(u_old);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(u_old);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"inital u is : \n");CHKERRQ(ierr);
  ierr = VecView(u_old, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

  /*get exact solution*/
  if (rank == 0)
  {
    for (i = 0; i < n; i++)
    {
      xi     = i*dx;
      vec_ue = sin(pi*xi)/pi/pi;
      ierr   = VecSetValues(u_exact, 1, &i, &vec_ue, INSERT_VALUES);CHKERRQ(ierr);
    }
  }
  ierr = VecAssemblyBegin(u_exact);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(u_exact);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"exact solution is : \n");CHKERRQ(ierr);
  ierr = VecView(u_exact, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

   /*set value for vec f*/
  if(rank == 0)
  {
    i         = 0; 
    ierr      = VecSetValues(f, 1, &i, &value_f, INSERT_VALUES);CHKERRQ(ierr);
    for(i = 1; i < n; i++)
    {
      xi        = i*dx;
      value_f   = sin(xi*pi)*dt/rho/c;
      ierr      = VecSetValues(f, 1, &i, &value_f, INSERT_VALUES);CHKERRQ(ierr);
    }
  }

  ierr = VecAssemblyBegin(f);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(f);CHKERRQ(ierr);

  #ifdef IMPLICIT
  /*create a linear solver */
  ierr = PetscPrintf(PETSC_COMM_WORLD,"It's implicit scheme\n");CHKERRQ(ierr);
  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);

  ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
 
  ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
  ierr = PCSetType(pc,PCJACOBI);CHKERRQ(ierr);
  ierr = KSPSetTolerances(ksp,1.e-7,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);

  ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
  ierr = PCSetFromOptions(pc);CHKERRQ(ierr);
  ierr = PetscViewerHDF5Open(PETSC_COMM_WORLD,"./h5/implicit.h5",FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) u, "restart_Vec");CHKERRQ(ierr);
  
  while(PetscAbsReal(norm)> tol && step < steps)
  {
    ierr = VecAXPY(u_old,1.0,f);CHKERRQ(ierr);
    ierr = KSPSolve(ksp,u_old,u);CHKERRQ(ierr);
    ierr = VecCopy(u,u_old);CHKERRQ(ierr);
    ierr = VecNorm(u,NORM_2,&norm);CHKERRQ(ierr);
    /*write value of u to uold.h5*/ 
    if(0 == step % 10)
    {
      ierr = VecView(u,viewer);CHKERRQ(ierr);
    }
    step++;
  }
  ierr = PetscPrintf(PETSC_COMM_WORLD,"step=%d\n",step);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  #endif

  #ifdef EXPLICIT
  ierr = PetscPrintf(PETSC_COMM_WORLD,"It's explicit scheme\n");CHKERRQ(ierr);
  ierr = PetscViewerHDF5Open(PETSC_COMM_WORLD,"./h5/explicit.h5",FILE_MODE_WRITE,&viewer);;CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) u, "restart_Vec");CHKERRQ(ierr);
  while(PetscAbsReal(norm) > tol && step < steps)
  {
    ierr = MatMult(A,u_old,u);CHKERRQ(ierr);
    ierr = VecAXPY(u,(PetscScalar)1.0,f);CHKERRQ(ierr);
    ierr = VecCopy(u,u_old);CHKERRQ(ierr);
    ierr = VecNorm(u,NORM_2,&norm);CHKERRQ(ierr);
    /*write value of u to h5*/
    if(step % 10 == 0)
    {
      ierr = VecView(u,viewer);CHKERRQ(ierr);
    }
    step++;
  }
  ierr = PetscPrintf(PETSC_COMM_WORLD,"step=%d\n",step);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  #endif

  ierr = PetscPrintf(PETSC_COMM_WORLD,"---------------------\n");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"numerical result u is : \n");CHKERRQ(ierr);
  ierr = VecView(u_old, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"The program is finished\n");CHKERRQ(ierr);

  ierr = VecDestroy(&u_old);CHKERRQ(ierr); 
  ierr = VecDestroy(&u);CHKERRQ(ierr);
  ierr = VecDestroy(&f);CHKERRQ(ierr);
  ierr = VecDestroy(&u_exact);CHKERRQ(ierr);
  // ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  #ifdef IMPLICIT
  ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
  #endif

  ierr = PetscFinalize();
  return ierr;
}