#ifndef COMMON_H
#define COMMON_H

#include <petscksp.h>
#include <fftw3-mpi.h>
#include <vector>
#include <math.h>
#include <iostream>
#include <complex.h>
#include <stdio.h>
//#include <string>
//#include <sstream>


#define _USE_MATH_DEFINES
#define ROOT 0
#define ZERO 1.0e-10
#define PI   (3.14159265359)
#define PI4  (M_PI*4)

/*
//In Optics
#define C   (3*1.0e+10)   // Speed of light cm/s
#define F   (1.0e+14)     // Frequency Hz
#define K   (2*PI*F/C)    // Wave number K = 2PI/lambda
*/
//Acoustic waves
#define C   (34400)       // Speed of light cm/s
#define F   (1000)        // Frequency Hz
#define K   (2*PI*F/C)    // Wave number K = 2PI/lambda


using namespace std;


template<class Vector,class Real, class Integer>
Real Norm(Vector v, Integer norm)
{
    Real sum = 0;
    Integer n = v.size();

    for(Integer i=0;i<n;i++)
    {
        sum += pow(v[i],norm);
    }
    return pow(sum,1.0/norm);
}

template<class Real, class Vector, class Integer>
Real Dot(Vector v, Vector u)
{
    Real sum = 0;
    Integer n = v.size();

    for(Integer i=0;i<n;i++)
    {
        sum += v[i]*u[i];
    }
    return sum;
}

template<class Real, class Complex>
Real cnorm(Complex c)
{
    return sqrt(creal(c)*creal(c)+cimag(c)*cimag(c));
}

template<class Integer>
bool is_nth_power(Integer a, Integer n) 
{
  if(n <= 0)
    return false;
  if((a < 0) && (n % 2 == 0))
    return false;
  a = abs(a);

  Integer b = pow(a, 1. / n);
  return (pow((double) b, n) == a || pow((double) (b+1), n) == a);
}

template <class Complex, class Real, class Long, class Integer> 
void ValidateInput(Real& Kappa, Real& VolQ, Real& ParRadius, Real& ParticleDistance, Long& TotalParticles,\
                  Complex& OriginalRefractionCoef, Complex& DesiredRefractionCoef, Real& Distribution, Integer& TotalCubes, Integer& TotalCollocationPoints)
{
  Real DomainSize;
  Real SubcubeSize;

  if( (Kappa<0) || (Kappa>=1) )
  {
     Kappa = 0.99;
     PetscPrintf(PETSC_COMM_WORLD,"\nWARNING: Kappa should be in (0,1). Set the default value!");
  }
  if(VolQ<=0)
  {
     VolQ = 1;
     PetscPrintf(PETSC_COMM_WORLD,"\nWARNING: Volume of the domain containing particles should be positive. Set the default value!");
  }
  DomainSize = cbrt(VolQ);

  if(TotalParticles<=0)
  {
     TotalParticles = pow(80,3);    
     PetscPrintf(PETSC_COMM_WORLD,"\nWARNING: Total number of particles should be positive. Set the default value!");
  }
  if( (ParRadius<=0) || (ParRadius>=DomainSize) )
  {
     ParRadius = pow(1.0/TotalParticles,1.0/(2-Kappa));
     PetscPrintf(PETSC_COMM_WORLD,"\nWARNING: Particle radius should be positive and less than the domain size. Set the default value!");
  }
  if( (ParticleDistance<=0) || (ParticleDistance>=DomainSize) || (ParticleDistance<ParRadius) )
  {
     ParticleDistance = pow(ParRadius,(2-Kappa)/3);   //ParticleDistance = O(ParticleRadius^(1/3))  
     PetscPrintf(PETSC_COMM_WORLD,"\nWARNING: Distance between neighboring particles should be positive, greater than the radius of a particle, and less than the domain size. Set the default value!");
  }
  if(cnorm<PetscReal,PetscScalar>(OriginalRefractionCoef)<=0)
  {
     OriginalRefractionCoef = 1;
     PetscPrintf(PETSC_COMM_WORLD,"\nWARNING: Original refraction should not be zero. Set the default value!");
  }
  if(cnorm<PetscReal,PetscScalar>(DesiredRefractionCoef)<=0)
  {
     DesiredRefractionCoef = -1+I*0.001;
     PetscPrintf(PETSC_COMM_WORLD,"\nWARNING: Desired refraction should not be zero. Set the default value!");
  }
  if(Distribution<=0)
  {
     Distribution = TotalParticles*pow(ParRadius,2-Kappa)*VolQ;
     PetscPrintf(PETSC_COMM_WORLD,"\nWARNING: Distribution of particles should be positive. Set the default value!");
  }

  SubcubeSize = DomainSize/cbrt(TotalCubes);
  if( (TotalCubes<=0) || (TotalCubes>TotalParticles) || !is_nth_power<Integer>(TotalCubes,3) || (SubcubeSize<ParticleDistance))
  {
     TotalCubes = pow(round(pow(TotalParticles,1.0/3)/10),3);
     if(TotalCubes<=0)
        TotalCubes = 1;
     PetscPrintf(PETSC_COMM_WORLD,"\nWARNING: Total number of subcubes should be positive, cubic and less than total particles. The size of a subcube must be greater than the distance between neighboring particles. Set the default value!");
  }
  if((TotalCollocationPoints>TotalParticles) || (TotalCollocationPoints<TotalCubes) || !is_nth_power<Integer>(TotalCollocationPoints,3))
  {
     TotalCollocationPoints = TotalCubes;
     PetscPrintf(PETSC_COMM_WORLD,"\nWARNING: Total number of collocation points should be positive, cubic, and in (total cubes,total particles). Set the default value!");
  }
}

template <class RealVector, class Complex, class Real, class Long, class Integer> 
void Output(Real kappa, Real VolQ, Real ParRadius, Real ParDist, Long TotalParticles, RealVector WaveDirection,\
            Complex OriginalRefractionCoef, Complex DesiredRefractionCoef, Real Distribution, Integer TotalCubes, Complex BoundaryImpedance, Integer& TotalCollocationPoints)
{
    //ostringstream stream;
    string s;

    cout.precision(6);
    cout<< "\n----------------------- SOLVING PARTICLE SCATTERING PROBLEM -----------------------\n";

    cout<< "\nSpeed of wave, v:\t\t\t\t\t";
    cout<< C;
    cout<< "\nFrequency, f:\t\t\t\t\t\t";
    cout<< F;
    cout<< "\nWave number, k:\t\t\t\t\t\t";
    cout<< K;    
    cout<< "\nDirection of plane wave, alpha:\t\t\t\t(";
    cout<< WaveDirection[0]<<", "<<WaveDirection[1]<<", "<<WaveDirection[2]<<")";
    cout<< "\nKappa:\t\t\t\t\t\t\t";
    cout<< kappa;
    cout<< "\nVolume of the domain that contains all particles, |D|:\t";
    cout<< VolQ;
    cout<< "\nOriginal refraction coefficient, n0:\t\t\t";
    s = (cimag(OriginalRefractionCoef)>=0)?"+":"";
    cout<< creal(OriginalRefractionCoef)<<s<<cimag(OriginalRefractionCoef)<<"i";
    cout<< "\nDesired refraction coefficient, n:\t\t\t";
    s = (cimag(DesiredRefractionCoef)>=0)?"+":"";
    cout<< creal(DesiredRefractionCoef)<<s<<cimag(DesiredRefractionCoef)<<"i";    
    cout<< "\nDistribution of particles, N:\t\t\t\t";
    cout<< Distribution;
    cout<< "\nBoundary impedance, zeta:\t\t\t\t";
    s = (cimag(BoundaryImpedance)>=0)?"+":"";
    cout<< creal(BoundaryImpedance)<<s<<cimag(BoundaryImpedance)<<"i"; 
    cout<< "\nRadius of one particle, a:\t\t\t\t";
    cout<< scientific<<ParRadius;
    cout<< "\nDistance between two neighboring particles, d:\t\t";
    cout<< scientific<<ParDist;
    cout<< "\nNumber of particles, M:\t\t\t\t\t";
    cout<< TotalParticles<<" ("<<scientific<<Real(TotalParticles)<<")";
    cout<< "\nNumber of sub cubes after partitioning the domain, P:\t";
    cout<< TotalCubes;
    cout<< "\nTotal collocation points for solving integral equation:\t";
    cout<< TotalCollocationPoints;
    
    cout<<endl;
} 

#undef __FUNCT__
#define __FUNCT__ "VecView3D"
PetscErrorCode VecView3D(Vec x,PetscInt dim)
{
    PetscInt    slap1,slap2,j,s,t;
    PetscInt    Istart,Iend;
    PetscScalar ***xa;
    string str;

    PetscFunctionBegin;
    
    VecGetOwnershipRange(x,&Istart,&Iend);
    slap1 = Istart/(dim*dim);
    slap2 = Iend/(dim*dim);
    VecGetArray3d(x,dim,dim,slap2-slap1,0,0,slap1,&xa);

    for(t=slap1;t<slap2;t++)
    {
       for(j=0;j<dim;j++)
       {
	     for(s=0;s<dim;s++)
	     {				
		 str = (cimag(xa[j][s][t])>=0)?"+":"";
		 cout<<creal(xa[j][s][t])<<str<<cimag(xa[j][s][t])<<"i\t";	  		
	     }
	     cout<<"\n";
       }
	cout<<"\n";
    }
    VecRestoreArray3d(x,dim,dim,slap2-slap1,0,0,slap1,&xa);

    PetscFunctionReturn(0);
}

/* SLOW - USES LESS MEMORY
#undef __FUNCT__
#define __FUNCT__ "VecCopyN"
inline const PetscErrorCode VecCopyN(Vec xx,Vec yy,PetscInt N)
{
    PetscInt    j;
    PetscInt    Istart,Iend;
    PetscScalar tmp;

    PetscFunctionBegin;

    VecGetOwnershipRange(xx,&Istart,&Iend);
    if(Istart>=0 && Istart<=N)
    {
	if(Iend<0)
	   Iend = N;
       //VecZeroEntries(yy);
       Iend = (Iend<N)?(Iend):(N);
       for(j=Istart;j<Iend;j++)
       {
	  VecGetValues(xx,1,&j,&tmp);
         VecSetValues(yy,1,&j,&tmp,INSERT_VALUES);		  	
       }
    }
    VecAssemblyBegin(yy);
    VecAssemblyEnd(yy);

    PetscFunctionReturn(0);
}
*/
/*
#undef __FUNCT__
#define __FUNCT__ "VecCopyN"
inline const PetscErrorCode VecCopyN(Vec xx,Vec yy,PetscInt N)
{
    PetscInt    j,t;
    PetscInt    Istart,Iend,*pos;
    PetscScalar *val;

    PetscFunctionBegin;

    VecGetOwnershipRange(xx,&Istart,&Iend);
    if(Istart>=0 && Istart<=N)
    {
	if(Iend<0)
	   Iend = N;

       Iend = (Iend<N)?(Iend):(N);
       PetscMalloc((Iend-Istart)*sizeof(PetscScalar),&val); 
       PetscMalloc((Iend-Istart)*sizeof(PetscInt),&pos);
	t = 0;
       for(j=Istart;j<Iend;j++)
       {
         pos[t] = t+Istart;
         t++;	  	
       }
       VecGetValues(xx,Iend-Istart,pos,val);
       VecSetValues(yy,Iend-Istart,pos,val,INSERT_VALUES);
       PetscFree(val);
       PetscFree(pos);
    }
    VecAssemblyBegin(yy);
    VecAssemblyEnd(yy);

    PetscFunctionReturn(0);
}
*/

#undef __FUNCT__
#define __FUNCT__ "VecCopyN"
inline const PetscErrorCode VecCopyN(Vec xx,Vec yy,PetscInt N)
{
    PetscInt    j,t;
    PetscInt    Istart,Iend,*pos;
    PetscScalar *xa;

    PetscFunctionBegin;

    VecGetOwnershipRange(xx,&Istart,&Iend);
    if(Istart>=0 && Istart<=N)
    {
	if(Iend<0)
	   Iend = N;

       Iend = (Iend<N)?(Iend):(N);
       PetscMalloc((Iend-Istart)*sizeof(PetscInt),&pos);
	t = 0;
       for(j=Istart;j<Iend;j++)
       {
         pos[t] = t+Istart;
         t++;	  	
       }
	VecGetArray(xx,&xa);
       VecSetValues(yy,Iend-Istart,pos,xa,INSERT_VALUES);
	VecRestoreArray(xx,&xa);
       PetscFree(pos);      
    }
    VecAssemblyBegin(yy);
    VecAssemblyEnd(yy);

    PetscFunctionReturn(0);
}

void checkTime(time_t StartTime, time_t EndTime, bool final)
{
    int hour, min, sec, time;

    if(final)
    {
	cout<<"\nStarted on:\t"<<ctime(&StartTime);
	cout<<"Finished on:\t"<<ctime(&EndTime);   
    } 	
    time = EndTime-StartTime;
    hour=time/3600;
    time=time%3600;
    min=time/60;
    time=time%60;
    sec=time;
    cout<<"Elapsed time:\t"<<hour<<":"<<min<<":"<<sec;
    cout<<endl;
}

/*
#undef __FUNCT__
#define __FUNCT__ "FFTPaddedGreenCube"
PetscErrorCode FFTPaddedGreenCube(PetscInt dim[])
{
    PetscScalar    tmp;
    PetscReal      scale,norm;
    PetscInt       Istart,Iend,j;
    PetscErrorCode ierr;
    PetscBool      flg = PETSC_FALSE;
    bool	     check = false;

    PetscFunctionBegin;

    //Create RealSpaceCube
    VecGetOwnershipRange(RealSpaceCube,&Istart,&Iend);
    for(j=Istart;j<Iend;j++)
    {
	tmp = Scat3DS.GreenVec[j];
       VecSetValues(RealSpaceCube,1,&j,&tmp,INSERT_VALUES);		  
    }
    VecAssemblyBegin(RealSpaceCube);
    VecAssemblyEnd(RealSpaceCube);

    flg  = PETSC_FALSE;
    ierr = PetscOptionsGetBool(PETSC_NULL,"-view_Green",&flg,PETSC_NULL);CHKERRQ(ierr);
    if (flg) {
       ierr = PetscPrintf(PETSC_COMM_WORLD,"\nPadded Green cube:\n");CHKERRQ(ierr);
       ierr = VecView(RealSpaceCube,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);    
	//VecView3D(RealSpaceCube,dim[0]);
    }

    // Apply FFTW_FORWARD for the Green cube:
    ierr = MatMult(FFTMatrix,RealSpaceCube,FreqSpaceCube);CHKERRQ(ierr);

    //FFTW computes an unnormalized DFT, need to scale
    scale = dim[0]*dim[1]*dim[2];
    scale = 1.0/(PetscReal)scale;
    ierr = VecScale(FreqSpaceCube,scale);CHKERRQ(ierr);

    flg  = PETSC_FALSE;
    ierr = PetscOptionsGetBool(PETSC_NULL,"-view_GreenFFT",&flg,PETSC_NULL);CHKERRQ(ierr);
    if (flg) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"\nFFT of padded Green cube:\n");CHKERRQ(ierr);
      ierr = VecView(FreqSpaceCube,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);    
      //VecView3D(FreqSpaceCube,dim[0]);
    }

    if(check) //CHECK INVERSED FFT:
    {
	// Apply FFTW_BACKWARD 
	ierr = MatMultTranspose(FFTMatrix,FreqSpaceCube,ReconstructCube);CHKERRQ(ierr);
	
	// Compare RealSpaceCube and ReconstructCube
	ierr = VecAXPY(ReconstructCube,-1.0,RealSpaceCube);CHKERRQ(ierr);
	ierr = VecNorm(ReconstructCube,NORMTYPE,&norm);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"\nError norm of |RealSpaceCube-ReconstructCube| = %e\n",norm);CHKERRQ(ierr);   
    }

    ierr = VecDuplicate(FreqSpaceCube,&FFTPaddedGreen);CHKERRQ(ierr);
    VecCopy(FreqSpaceCube,FFTPaddedGreen);

    PetscFunctionReturn(0);
}
*/

/*
#undef __FUNCT__
#define __FUNCT__ "FFTPaddedGreenCube"
PetscErrorCode FFTPaddedGreenCube(PetscInt dim[])
{
    PetscInt 	     slap1,slap2,j,s,t;
    PetscReal      scale;
    PetscReal      norm;
    PetscScalar    ***xa;
    PetscInt       Istart,Iend;
    PetscErrorCode ierr;
    PetscBool      flg = PETSC_FALSE;
    bool	     check = false;

    PetscFunctionBegin;

    //Create RealSpaceCube of size dim: 
    VecGetOwnershipRange(RealSpaceCube,&Istart,&Iend);
    slap1 = Istart/(dim[0]*dim[1]);
    slap2 = Iend/(dim[0]*dim[1]);
    ierr = VecGetArray3d(RealSpaceCube,dim[0],dim[1],slap2-slap1,0,0,slap1,&xa);CHKERRQ(ierr);
    for(j=0;j<dim[0];j++)
    {
       for(s=0;s<dim[1];s++)
       {
	     for(t=slap1;t<slap2;t++)
	     {
 		  xa[j][s][t] = Scat3DS.GreenCube[j][s][t];
	     }
       }
    }
    ierr = VecRestoreArray3d(RealSpaceCube,dim[0],dim[1],slap2-slap1,0,0,slap1,&xa);CHKERRQ(ierr);
    
    flg  = PETSC_FALSE;
    ierr = PetscOptionsGetBool(PETSC_NULL,"-view_Green",&flg,PETSC_NULL);CHKERRQ(ierr);
    if (flg) {
       ierr = PetscPrintf(PETSC_COMM_WORLD,"\nPadded Green cube:\n");CHKERRQ(ierr);
       ierr = VecView(RealSpaceCube,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);    
	//VecView3D(RealSpaceCube,dim[0]);
    }

    // Apply FFTW_FORWARD for the Green cube:
    ierr = MatMult(FFTMatrix,RealSpaceCube,FreqSpaceCube);CHKERRQ(ierr);

    //FFTW computes an unnormalized DFT, need to scale
    scale = dim[0]*dim[1]*dim[2];
    scale = 1.0/(PetscReal)scale;
    ierr = VecScale(FreqSpaceCube,scale);CHKERRQ(ierr);

    flg  = PETSC_FALSE;
    ierr = PetscOptionsGetBool(PETSC_NULL,"-view_GreenFFT",&flg,PETSC_NULL);CHKERRQ(ierr);
    if (flg) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"\nFFT of padded Green cube:\n");CHKERRQ(ierr);
      ierr = VecView(FreqSpaceCube,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);    
      //VecView3D(FreqSpaceCube,dim[0]);
    }

    if(check) //CHECK INVERSED FFT:
    {
	// Apply FFTW_BACKWARD 
	ierr = MatMultTranspose(FFTMatrix,FreqSpaceCube,ReconstructCube);CHKERRQ(ierr);
	
	// Compare RealSpaceCube and ReconstructCube
	ierr = VecAXPY(ReconstructCube,-1.0,RealSpaceCube);CHKERRQ(ierr);
	ierr = VecNorm(ReconstructCube,NORMTYPE,&norm);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"\nError norm of |RealSpaceCube-ReconstructCube| = %e\n",norm);CHKERRQ(ierr);   
    }

    ierr = VecDuplicate(FreqSpaceCube,&FFTPaddedGreen);CHKERRQ(ierr);
    VecCopy(FreqSpaceCube,FFTPaddedGreen);

    PetscFunctionReturn(0);
}
*/

/*
#undef __FUNCT__
#define __FUNCT__ "Convert1Dto3D"
inline const PetscErrorCode Convert1Dto3D(Vec xx,Vec yy,PetscInt N)
{
    PetscInt    j,s,t,m,dim[3],slap1,slap2;
    PetscInt    n,Istart,Iend;
    PetscScalar ***xa;
    Vec         vec;
    VecScatter  scat;
    IS          is;

    PetscFunctionBegin;

    VecGetSize(xx,&n);
    VecCreate(PETSC_COMM_SELF,&vec);
    VecSetSizes(vec,PETSC_DECIDE,n);
    VecSetFromOptions(vec);

    ISCreateStride(PETSC_COMM_WORLD,n,0,1,&is);
    VecScatterCreate(xx,PETSC_NULL,vec,is,&scat);
    VecScatterBegin(scat,xx,vec,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(scat,xx,vec,INSERT_VALUES,SCATTER_FORWARD);

    VecZeroEntries(yy);
    VecGetOwnershipRange(yy,&Istart,&Iend);
    dim[0]=2*N-2;
    dim[1]=dim[0];
    dim[2]=dim[0];
    slap1 = Istart/(dim[0]*dim[1]);
    slap2 = Iend/(dim[0]*dim[1]);
    VecGetArray3d(yy,dim[0],dim[1],slap2-slap1,0,0,slap1,&xa);
    for(j=0;j<N;j++)
    {
       for(s=0;s<N;s++)
       {
	     for(t=slap1;t<slap2;t++)
	     {
		if(t<N)
		{
		  m = Scat3DS.Index2Order(j,s,t);
		  VecGetValues(vec,1,&m,&xa[j][s][t]);
		}
	     }
       }
    }
    VecRestoreArray3d(yy,dim[0],dim[1],slap2-slap1,0,0,slap1,&xa);

    VecDestroy(&vec);
    VecScatterDestroy(&scat);
    ISDestroy(&is);

    PetscFunctionReturn(0);
}
*/

/*
#undef __FUNCT__
#define __FUNCT__ "Convert3Dto1D"
inline const PetscErrorCode Convert3Dto1D(Vec xx,Vec yy,PetscInt N)
{
    PetscInt    j,s,t,m,dim[3],slap1,slap2;
    PetscInt    Istart,Iend;
    PetscScalar ***xa;

    PetscFunctionBegin;  

    VecGetOwnershipRange(xx,&Istart,&Iend);
    dim[0]=2*N-2;
    dim[1]=dim[0];
    dim[2]=dim[0];
    slap1 = Istart/(dim[0]*dim[1]);
    slap2 = Iend/(dim[0]*dim[1]);
    VecGetArray3d(xx,dim[0],dim[1],slap2-slap1,0,0,slap1,&xa);
    for(j=0;j<N;j++)
    {
       for(s=0;s<N;s++)
       {
	     for(t=slap1;t<slap2;t++)
	     {
		if(t<N)
		{
		   m = Scat3DS.Index2Order(j,s,t);
	          VecSetValues(yy,1,&m,&xa[j][s][t],INSERT_VALUES);		  
		}
	     }
       }
    }
    VecRestoreArray3d(xx,dim[0],dim[1],slap2-slap1,0,0,slap1,&xa);

    VecAssemblyBegin(yy);
    VecAssemblyEnd(yy);

    PetscFunctionReturn(0);
}
*/

#endif //COMMON_H
