/* Program usage:  mpiexec -n <procs> scat [-help] [all PETSc options] */

static char help[] = "Solves particle-scattering problem in parallel.\n\
References:\n\
[595]  A.G. Ramm,  Wave scattering by many small bodies and creating materials with a desired refraction coefficient, Afrika Matematika, 22, N1, (2011), 33-55.\n\
Input parameters include:\n\
  -a <particle_radius>         : particle radius\n\
  -d <particles_distance>      : distance between neighboring particles, 1 >> d >> a, default value cuberoot[a^(2-Kappa)]\n\
  -n <total_particles>         : total number of particles, n = O(1/(a^(2-Kappa)))\n\
  -p <total_cubes>             : total number of embeded small cubes in the domain containing all particles (for solving the reduced system)\n\
  -c <total_collocation_points>: total number of collocation points in the domain containing all particles (for solving the integral equation)\n\
  -k <kappa>                   : power constant with respect to the radius of particles, kappa is in [0,1), default value 0.99\n\
  -vol <volume>                : volume of the domain that contains all particles, default value 1\n\
  -ori <original_refraction>   : original refraction coefficient, default value 1\n\
  -des <desired_refraction>    : desired refraction coefficient, default value sqrt(0.2)\n\
  -dis <distribution>          : distribution of particles, default value 1 for uniform distribution\n\
  -view_RHS                    : write RHS vector to stdout\n\
  -view_solution               : write solution vector to stdout\n\
  -standard                    : solve the original system in the standard way, not using convolution\n\n";

/*
  Include "petscksp.h" so that we can use KSP solvers.  Note that this file
  automatically includes:
     petscsys.h    - base PETSc routines   petscvec.h - vectors
     petscmat.h    - matrices
     petscis.h     - index sets            petscksp.h - Krylov subspace methods
     petscviewer.h - viewers               petscpc.h  - preconditioners
*/

//#include <petscksp.h>
#include "scattering3DS.h"
#include "scattering3DP.h"
#include "scatteringIE_Riemann.h"

using namespace std;

#define NORMTYPE NORM_1

extern void GetTime(bool);


Scattering3DS<std::vector<PetscReal>,std::vector<PetscScalar>,PetscScalar,PetscReal,PetscInt,PetscInt> 		Scat3DS;
Scattering3DP<std::vector<PetscReal>,std::vector<PetscScalar>,PetscScalar,PetscReal,PetscInt,PetscInt> 		Scat3DP;
ScatteringIE_Riemann<std::vector<PetscReal>,std::vector<PetscScalar>,PetscScalar,PetscReal,PetscInt,PetscInt> 	ScatIE;
PetscLogDouble  mem;

//MPI:
int           rank;
int           size;
time_t        time1, time2;  //Time usage 


#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **args)
{
  Vec            FFTx,y,z,u0,u0P,u0IE,Distance,nParticlePerCube;        //approx solution, RHS, solution-difference 
  //Vec          FFTr,r,rP,rIE;                                         //residual 
  Mat            B,IE,AFFT;           		                        //linear system matrix 
  KSP            kspP,kspIE,kspFFT;      	                  	//linear solver context 
  //PC           pc;          		                                //Preconditioner context 
  PetscReal      norm, Tol=1e-3, aTol;          		           
  PetscInt       its;           		                        //number of iterations reached 
  PetscInt       Istart,Iend;
  PetscErrorCode ierr;
  KSPConvergedReason reason;
  KSPType 	 ksptype;
  PetscScalar    *xa;
  PetscBool      flg = PETSC_FALSE;
  
  //Scattering3D:
  PetscInt       M = 0;                 			  	   //Total particles
  PetscReal      a = 0;                 			  	   //Particle radius
  PetscReal      ParticleDistance = 0;
  PetscReal      VolQ = 0;              			  	   //Volume of the domain Q that contains all particles
  PetscReal      Kappa = -1;             			  	   //Power const with respect to the radius of particles: Kappa in [0,1)
  std::vector<PetscReal> WaveDirection(3,0);WaveDirection[0] = 1; 	   // WaveDirection is a unit vector that indicates the direction of plane wave
  PetscScalar    OriginalRefractionCoef = 0;
  PetscScalar    DesiredRefractionCoef = 0;
  PetscReal      Distribution = 0;
  PetscInt       TotalCubes = 0;
  PetscInt       N = 0;                 			  	   //Number of sub domains used for solving IE
  PetscScalar    tmp;
  PetscInt       s,j,HM;
  PetscReal      max;


#if defined(PETSC_USE_LOG)
  PetscLogStage  stage;
#endif    

  PetscInitialize(&argc,&args,(char *)0,help);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  //printf ("Hello from task %d! Total cores: %d\n",rank,size);

  time(&time1);

  ierr = PetscOptionsGetReal(PETSC_NULL,"-a",&a,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(PETSC_NULL,"-d",&ParticleDistance,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(PETSC_NULL,"-n",&M,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(PETSC_NULL,"-p",&TotalCubes,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(PETSC_NULL,"-c",&N,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(PETSC_NULL,"-k",&Kappa,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(PETSC_NULL,"-vol",&VolQ,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetScalar(PETSC_NULL,"-ori",&OriginalRefractionCoef,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetScalar(PETSC_NULL,"-des",&DesiredRefractionCoef,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(PETSC_NULL,"-dis",&Distribution,PETSC_NULL);CHKERRQ(ierr);

  //Scatering init 
  ValidateInput(Kappa,VolQ,a,ParticleDistance,M,OriginalRefractionCoef,DesiredRefractionCoef,Distribution,TotalCubes,N);
  Scat3DS.Input(a,Kappa,WaveDirection,ParticleDistance,M,OriginalRefractionCoef,DesiredRefractionCoef,Distribution,VolQ);
  Scat3DP.Input(a,Kappa,WaveDirection,ParticleDistance,M,OriginalRefractionCoef,DesiredRefractionCoef,Distribution,VolQ,TotalCubes,N); 
  ScatIE.Input(a,Kappa,WaveDirection,ParticleDistance,M,OriginalRefractionCoef,DesiredRefractionCoef,Distribution,VolQ,N); 
  Scat3DS.Init();
  Scat3DP.Init();
  ScatIE.Init();
  if(rank==0)
  {
      Output<std::vector<PetscReal>,PetscScalar,PetscReal,PetscInt,PetscInt>(Kappa,VolQ,a,ParticleDistance,M,WaveDirection,\
                                                                             OriginalRefractionCoef,DesiredRefractionCoef,Distribution,TotalCubes,Scat3DS.BoundaryImpedance,N);
  }

  PetscPrintf(PETSC_COMM_WORLD,"\nInitializing done:");
  PetscMemoryGetCurrentUsage(&mem);PetscPrintf(PETSC_COMM_WORLD,"\tMem used: %G \t",mem);GetTime(0);

  //----------------------------------------COMPUTE THE MATRIX NORM---------------------------------------------------
  ierr = PetscLogStageRegister("Norm", &stage);CHKERRQ(ierr);
  ierr = PetscLogStagePush(stage);CHKERRQ(ierr); 

  ierr = VecCreate(PETSC_COMM_WORLD,&u0);CHKERRQ(ierr);
  ierr = VecSetSizes(u0,PETSC_DECIDE,M);CHKERRQ(ierr);
  ierr = VecSetFromOptions(u0);CHKERRQ(ierr);

  ierr = VecGetOwnershipRange(u0,&Istart,&Iend);CHKERRQ(ierr);

  max = 0;
  HM = M/2;
  for(j = 0; j < HM; j++)
  {
        ierr = VecGetArray(u0,&xa);CHKERRQ(ierr);
     	 for(s = Istart; s <Iend ; s++)
     	 {           
            xa[s-Istart] = Scat3DS.CoefMat(j,s);      
     	 }
        if((Istart<=j) && (j<Iend))
           xa[j-Istart] = 1;
        ierr = VecRestoreArray(u0,&xa);CHKERRQ(ierr);
        
        ierr = VecNorm(u0,NORMTYPE,&norm);CHKERRQ(ierr);
        if(max < norm)
        {
           max = norm;
           PetscPrintf(PETSC_COMM_WORLD,"Matrix infinity norm:\t%G\t\t\tMax Row:\t%d\n",max,j);
        }
  }

  PetscPrintf(PETSC_COMM_WORLD,"Computing matrix norm:");
  PetscMemoryGetCurrentUsage(&mem);PetscPrintf(PETSC_COMM_WORLD,"\tMem used: %G \t",mem);GetTime(0);
  PetscPrintf(PETSC_COMM_WORLD,"Matrix infinity norm:\t%G\n",max);

  ierr = PetscLogStagePop();CHKERRQ(ierr); 

  GetTime(1);
  
  ierr = VecDestroy(&u0);CHKERRQ(ierr);

  PetscFinalize();
  return 0;
}

////////////////////////////////////////////////////////////////////////////////////

#undef __FUNCT__
#define __FUNCT__ "GetTime"
void GetTime(bool final)
{
    //Check total time used:
    if(rank==0)  
    {
       time(&time2);  
       checkTime(time1,time2,final);
    }
}