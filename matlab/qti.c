/* $Revision: 1.8.6.2 $ */
/*=========================================================
 * convec.c
 * example for illustrating how to use pass complex data 
 * from MATLAB to C and back again
 *
 * convolves  two complex input vectors
 *
 * This is a MEX-file for MATLAB.
 * Copyright 1984-2006 The MathWorks, Inc.
 *=======================================================*/
#include <complex.h>
#include <math.h>
#include "mex.h"

#define min(a,b) ((a)>(b)?(b):(a))
/* computational subroutine */
/*
void convec( double *xr, double *xi, mwSize nx,
             double *yr, double *yi, mwSize ny,
             double *zr, double *zi)
{
    mwSize i,j;
  
    zr[0]=0.0;
    zi[0]=0.0;
    // perform the convolution of the complex vectors 
    for(i=0; i<nx; i++) {
	for(j=0; j<ny; j++) {
	    *(zr+i+j) = *(zr+i+j) + *(xr+i) * *(yr+j) - *(xi+i) * *(yi+j);
	    *(zi+i+j) = *(zi+i+j) + *(xr+i) * *(yi+j) + *(xi+i) * *(yr+j);
	}
    }
}*/
__inline__ bool mxIsInt(const mxArray *pm){
	return (	(!mxIsComplex(pm)) ||
				mxGetN(pm)==1||
				mxGetM(pm)==1||
				mxGetClassID(pm)>=mxINT8_CLASS||
				mxGetClassID(pm)<=mxUINT64_CLASS ) ;
}

/* The gateway routine. */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
	int mxVLoc, mxSVLoc, da, db = 0, len;
	double *vr, *vi, sv[min(da,db)], *svr; 
	complex v[da*db],ia[da*min(da,db)],ib[db*min(da,db)];
	register int i,j;
	/* check for the proper number of arguments */
	switch(nrhs){
		case 3:
		mxVLoc = 2;
		if(!mxIsInt(prhs[0])) {
			mexErrMsgTxt("First input must be an integer if there are three inputs");}
		da = mxGetScalar(prhs[0]);
		if(!mxIsInt(prhs[1])) {
			mexErrMsgTxt("Second input must be an integer if there are three inputs");}
		db = mxGetScalar(prhs[1]);
		break;
		case 2:
		mxVLoc = 1;
		if(!mxIsInt(prhs[0])) {
			mexErrMsgTxt("First input must be an integer if there are two inputs");}
		da = mxGetScalar(prhs[0]);
		break;
		case 1:
		mxVLoc = 0;
		da = db = 2;
		break;
		default:
		mexPrintf("\nnrhs:%d\n", nrhs);
		mexErrMsgTxt("One or three inputs required.");
		break;
	}
	switch(nlhs){
		case 1:
		mxSVLoc = 0;	
		break;
		case 3:
		mxSVLoc = 0;	
		break;
		default:mexErrMsgTxt("One or three outputs required.");
		break;
	}
	/*Check the input is a column vector*/
	if(mxGetN(prhs[mxVLoc]) != 1 ){
		mexErrMsgTxt("Inputs must be a column vector.");}
	/* Check that both inputs are complex*/
	switch(mxGetClassID(prhs[mxVLoc])){
		case mxDOUBLE_CLASS:
		case mxSINGLE_CLASS:
		break;
		case mxINT8_CLASS:
        case mxUINT8_CLASS:
        case mxINT16_CLASS:
        case mxUINT16_CLASS:
        case mxINT32_CLASS:
        case mxUINT32_CLASS:
        case mxINT64_CLASS:
        case mxUINT64_CLASS:
		break;
		default:mexErrMsgTxt("Illegal input type; must be float or int\n");
		break;
	}
	
	/* get the length of the input vector */
	len = mxGetM(prhs[mxVLoc]);
	if(len != da*db){
		if(mxVLoc == 2){
			db = len / da;
		}else{
			mexErrMsgTxt("The length of the array doesn't equal the product of the two dimensions\n");
		}
	}
	
	/* get pointers to the real and imaginary parts of the inputs */
	vr = mxGetPr(prhs[mxVLoc]);
	if(mxIsComplex(prhs[mxVLoc])){
		register int i,j;
		vi = mxGetPi(prhs[mxVLoc]);
		for(i=0;i<da;i++){
			for(j=0;j<db;j++){
				v[i*da+j]=vr[i*da+j]+I*vi[i*da+j];}}
	}else{
		register int i,j;
		for(i=0;i<da;i++){
			for(j=0;j<db;j++){
				v[i*da+j]=vr[i*da+j];}}
	}
	schmidt_decomposition(da,db,v,sv,ia,ib);
	plhs[0] = mxCreateDoubleMatrix(min(da,db), min(da,db), mxREAL);
	svr = mxGetPr(plhs[mxSVLoc]);
	for(i=0;i<min(da,db);i++){
		svr[i*da+i] = sv[i];}
/*
	zr = mxGetPr(plhs[0]);
	zi = mxGetPi(plhs[0]);
*/
	/* create a new array and set the output pointer to it */
/*
	cols = nx + ny - 1;
*/
	
	/* call the C subroutine */
/*	convec(xr, xi, nx, yr, yi, ny, zr, zi);*/
	
	return;
}
