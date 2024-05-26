#include <math.h>
//#include <stdio.h>
#include "mex.h"

/*
 * powdersimu.c - simulation of the powder average for g- and dH anisotropies
 *
 * multiplies an input scalar times an input matrix and outputs a
 * matrix 
 *
 * This is a MEX-file for MATLAB.
 * Copyright 13.03.2010 - 24.03.2014   Anton Potocnik @ IJS F5  
 */

/* $Revision: 3.0 $ */

void powdersimudL(double *Y, double *X, mwSize N, double freq, double gx, double gy, double gz, double dHx, double dHy, double dHz, double phase, mwSize m)
{
  unsigned int i,j,k;
  double fi, costh, sinth, cosfi, sinfi, g, dH, Hc, cosph, sinph, maxi=-100000, temp, temp1;
  
  cosph = cos(phase);
  sinph = sin(phase);
	
  for (i=0; i<N; i++) {   /* fi */
	for (j=0; j<N; j++) {   /* th */
		fi = 6.283185*i/(N-1);
		cosfi = cos(fi);
		sinfi = sin(fi);
		costh = 2.0*j/(N-1) - 1;
		sinth = sin(acos(costh));
		
		g = sqrt(gx*sinth*cosfi*gx*sinth*cosfi + gy*sinth*sinfi*gy*sinth*sinfi + gz*costh*gz*costh);
		dH = sqrt(dHx*sinth*cosfi*dHx*sinth*cosfi + dHy*sinth*sinfi*dHy*sinth*sinfi + dHz*costh*dHz*costh);
		Hc = 0.0714477287*freq/g;
        
        for (k=0; k<m; k++)
			*(Y+k) += 0.6366198*(dH*cosph - 2*(*(X+k)-Hc)*sinph)/(4*(*(X+k)-Hc)*(*(X+k)-Hc) + dH*dH);
    }
  }
  /* DERIVATION */
  temp = *(Y+0);
  for (k=1; k<m; k++) {
      *(Y+k-1) = (*(Y+k) - temp)/(*(X+k) - *(X+k-1)); /* n(k-1) = n(k) - n(k-1)*/
      if (*(X+k) - *(X+k-1) > 3*temp1) {
          *(Y+k-1) = *(Y+k-2);           /*!!! A possible source of error, but very unlikely !!!*/
      }
      temp = *(Y+k);
      temp1 = *(X+k) - *(X+k-1);
  }
  *(Y+k-1)=*(Y+k-2); /* last value is the same as one before last */
  
  for (k=0; k<m; k++) *(Y+k) /= N*N;
}



void powdersimuL(double *Y, double *X, mwSize N, double freq, double gx, double gy, double gz, double dHx, double dHy, double dHz, double phase, mwSize m)
{
  unsigned int i,j,k;
  double fi, costh, sinth, cosfi, sinfi, g, dH, Hc, cosph, sinph, maxi=-100000, temp, temp1;
  
  cosph = cos(phase);
  sinph = sin(phase);
	
  for (i=0; i<N; i++) {   /* fi */
	for (j=0; j<N; j++) {   /* th */
		fi = 6.283185*i/(N-1);
		cosfi = cos(fi);
		sinfi = sin(fi);
		costh = 2.0*j/(N-1) - 1;
		sinth = sin(acos(costh));
		
		g = sqrt(gx*sinth*cosfi*gx*sinth*cosfi + gy*sinth*sinfi*gy*sinth*sinfi + gz*costh*gz*costh);
		dH = sqrt(dHx*sinth*cosfi*dHx*sinth*cosfi + dHy*sinth*sinfi*dHy*sinth*sinfi + dHz*costh*dHz*costh);
		Hc = 0.0714477287*freq/g;
        
        for (k=0; k<m; k++)
			*(Y+k) += 0.6366198*(dH*cosph - 2*(*(X+k)-Hc)*sinph)/(4*(*(X+k)-Hc)*(*(X+k)-Hc) + dH*dH);
    }
  }

  
  for (k=0; k<m; k++) *(Y+k) /= N*N;
}





void powdersimudG(double *Y, double *X, mwSize N, double freq, double gx, double gy, double gz, double dHx, double dHy, double dHz, double phase, mwSize m)
{
  unsigned int i,j,k;
  double fi, costh, sinth, cosfi, sinfi, g, dH, Hc, cosph, sinph, maxi=-100000, temp, temp1;
  //FILE *fp;
  
  double *Yhil; //,*Xtab,*Ytab;
  Yhil = mxMalloc(m*sizeof(double)+1);
  
  /*fp=fopen("GaussDispersion.dat", "r");
  if(fp==NULL) 
    mexErrMsgTxt("GaussDispersion.dat not found!");
  
  fscanf(fp,"%d\n",&A);
  fscanf(fp,"%d\n",&B);
  
  Ytab = mxMalloc(A*sizeof(double)+1);
  Xtab = mxMalloc(A*sizeof(double)+1);
  
  for (k=0; k<A; k++) {
      *(Xtab+k) = B*(-0.5+k/(A-1));
      fscanf(fp,"%f\n",(Ytab+k));
      printf("%f\t%f\n",*(Xtab+k),*(Ytab+k));
  }
  
  fclose(fp);*/
  
  for (k=0; k<m; k++) { // set all Yhil to zero
      *(Y+k)=0;
      *(Yhil+k)=0;
  }
  
  cosph = cos(phase);
  sinph = sin(phase);
	
  for (i=0; i<N; i++) {   /* fi */
	for (j=0; j<N; j++) {   /* th */
		fi = 6.283185*i/(N-1);
		cosfi = cos(fi);
		sinfi = sin(fi);
		costh = 2.0*j/(N-1) - 1;
		sinth = sin(acos(costh));
		
		g = sqrt(gx*sinth*cosfi*gx*sinth*cosfi + gy*sinth*sinfi*gy*sinth*sinfi + gz*costh*gz*costh);
		dH = sqrt(dHx*sinth*cosfi*dHx*sinth*cosfi + dHy*sinth*sinfi*dHy*sinth*sinfi + dHz*costh*dHz*costh);
		Hc = 0.0714477287*freq/g;
        //printf("test %f\t%f\t%f\t%f\t%f\n",fi,cosfi,sinfi,costh,sinth);
        //temp = (*(X+50)-Hc)/dH;
        //printf("test %f\t%f\t%f\t%f\n",dH,Hc,temp*temp,0.3989423/dH*exp(-temp*temp));
        for (k=0; k<m; k++) {
            *(Y+k) += 0.3989423/dH*exp(-(*(X+k)-Hc)*(*(X+k)-Hc)/dH/dH/2);
            //*(Y+k) += 0.3989423/dH*exp(-(*(X+k)-Hc)*(*(X+k)-Hc)/dH/dH/2);
        }
    }
  }
  
  /* DERIVATION */
  temp = *(Y+0);
  for (k=1; k<m-1; k++) {
      temp1 = *(Y+k);
      *(Y+k) = (*(Y+k+1) - temp)/(*(X+k+1) - *(X+k-1)); // n(k) = n(k+1) - n(k-1)/
      //if (*(X+k+1) - *(X+k-1) > 6*temp1) {
      //    *(Y+k) = *(Y+k-1);           //!!! A possible source of error, but very unlikely !!!
      //}
      temp = temp1;
      //temp1 = *(X+k) - *(X+k-1);
  }
  *(Y+k-1)=*(Y+k-2); // last value is the same as one before 
  *(Y+0)=*(Y+1); // first value is the same as one after 
  
  
  /* HILBERT TRANSFORM */   
  for (i=0; i<m; i++) { // For every element
      temp=0;
      for (k=1; k<m; k++) { // Integration
          if (i==k) continue;
            temp -= *(Y+k)/(*(X+k)-*(X+i))*(*(X+k) - *(X+k-1));
      }
      *(Yhil+i) = temp/3.1415926536;
  }
  
  /* MIX DISPERSION AND ABSORPTION */
  for (k=0; k<m; k++)  
      *(Y+k) = (cosph*(*(Y+k)) - sinph*(*(Yhil+k)));

  /* NORMALIZE */
  for (k=0; k<m; k++) 
      *(Y+k) /= N*N;
  
  mxFree(Yhil);
}



void powdersimuG(double *Y, double *X, mwSize N, double freq, double gx, double gy, double gz, double dHx, double dHy, double dHz, double phase, mwSize m)
{
  unsigned int i,j,k;
  double fi, costh, sinth, cosfi, sinfi, g, dH, Hc, cosph, sinph, maxi=-100000, temp, temp1;
  //FILE *fp;
  
  double *Yhil; //,*Xtab,*Ytab;
  Yhil = mxMalloc(m*sizeof(double)+1);
  
  
  for (k=0; k<m; k++) { // set all Yhil to zero
      *(Y+k)=0;
      *(Yhil+k)=0;
  }
  
  cosph = cos(phase);
  sinph = sin(phase);
	
  for (i=0; i<N; i++) {   /* fi */
	for (j=0; j<N; j++) {   /* th */
		fi = 6.283185*i/(N-1);
		cosfi = cos(fi);
		sinfi = sin(fi);
		costh = 2.0*j/(N-1) - 1;
		sinth = sin(acos(costh));
		
		g = sqrt(gx*sinth*cosfi*gx*sinth*cosfi + gy*sinth*sinfi*gy*sinth*sinfi + gz*costh*gz*costh);
		dH = sqrt(dHx*sinth*cosfi*dHx*sinth*cosfi + dHy*sinth*sinfi*dHy*sinth*sinfi + dHz*costh*dHz*costh);
		Hc = 0.0714477287*freq/g;
        //printf("test %f\t%f\t%f\t%f\t%f\n",fi,cosfi,sinfi,costh,sinth);
        //temp = (*(X+50)-Hc)/dH;
        //printf("test %f\t%f\t%f\t%f\n",dH,Hc,temp*temp,0.3989423/dH*exp(-temp*temp));
        for (k=0; k<m; k++) {
            *(Y+k) += 0.3989423/dH*exp(-(*(X+k)-Hc)*(*(X+k)-Hc)/dH/dH/2);
            //*(Y+k) += 0.3989423/dH*exp(-(*(X+k)-Hc)*(*(X+k)-Hc)/dH/dH/2);
        }
    }
  }
  
  
  /* HILBERT TRANSFORM */   
  for (i=0; i<m; i++) { // For every element
      temp=0;
      for (k=1; k<m; k++) { // Integration
          if (i==k) continue;
            temp -= *(Y+k)/(*(X+k)-*(X+i))*(*(X+k) - *(X+k-1));
      }
      *(Yhil+i) = temp/3.1415926536;
  }
  
  /* MIX DISPERSION AND ABSORPTION */
  for (k=0; k<m; k++)  
      *(Y+k) = (cosph*(*(Y+k)) - sinph*(*(Yhil+k)));

  /* NORMALIZE */
  for (k=0; k<m; k++) 
      *(Y+k) /= N*N;
  
  mxFree(Yhil);
}



/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  double *Y,*X;
  double freq, gx, gy, gz, dHx, dHy, dHz, phase;
  int N, type;
  mwSize mrows, ncols, m, k;
  
  /*  check for proper number of arguments */
  /* NOTE: You do not need an else statement when using mexErrMsgTxt
     within an if statement, because it will never get to the else
     statement if mexErrMsgTxt is executed. (mexErrMsgTxt breaks you out of
     the MEX-file) */
  if(nrhs!=11) 
    mexErrMsgTxt("Ten inputs required.");
  if(nlhs!=1) 
    mexErrMsgTxt("One output required.");
  
  /* check to make sure the first input argument is a vector */
  if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ||
      mxGetN(prhs[0])*mxGetM(prhs[0])==1 ) {
    mexErrMsgTxt("Input x must be a vector.");
  }

  /*  create a pointer to the input matrix X */
  X = mxGetPr(prhs[0]);
  /*  get the dimensions of the matrix input X */
  mrows = mxGetM(prhs[0]);
  ncols = mxGetN(prhs[0]);
  
  if (mrows>=ncols) 
	m = mrows;
  else
	m = ncols;
	
  /*  get the scalars */
  N = mxGetScalar(prhs[1]); 
  freq = mxGetScalar(prhs[2]);
  gx = mxGetScalar(prhs[3]);
  gy = mxGetScalar(prhs[4]);
  gz = mxGetScalar(prhs[5]);
  dHx = mxGetScalar(prhs[6]);
  dHy = mxGetScalar(prhs[7]);
  dHz = mxGetScalar(prhs[8]);
  phase = mxGetScalar(prhs[9]);
  type = (int)mxGetScalar(prhs[10]);

  /*  set the output pointer to the output matrix */
  plhs[0] = mxCreateDoubleMatrix(mrows,ncols, mxREAL);  
  /*  create a C pointer to a copy of the output matrix */
  Y = mxGetPr(plhs[0]);
  
  /*  call the C subroutine 1 - Lorentz, 2- Gauss*/ 
  switch(type) {
      case 1:
          powdersimudL(Y, X, N, freq, gx, gy, gz, dHx, dHy, dHz, phase, m);
          break;
          
      case 2:
          powdersimudG(Y, X, N, freq, gx, gy, gz, dHx, dHy, dHz, phase, m);
          break;
    
      case 3:
          powdersimuL(Y, X, N, freq, gx, gy, gz, dHx, dHy, dHz, phase, m);
          break;
           
      case 4:
          powdersimuG(Y, X, N, freq, gx, gy, gz, dHx, dHy, dHz, phase, m);
          break;
          
      default:
          printf("Wrong type! Choose 1 for dLor, 2 for dGauss, 3 for Lor, and 4 for Gauss");
          return;
  }                      
                       
}
