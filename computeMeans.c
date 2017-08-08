/* 
 * computeMeans.c
 * Compile:
 *  mex computeMeans.c
 *or
 *  mex computeMeans.c -DUSE_BLAS -lmwblas
 * (both about the same speed)
 *
 * Usage:
 *      centers = computeMeans( X, k, int32(ind-1) );
 * where
 *      X is p x n, k is an integer, and ind is 1 x n
 *     (in the above example, we have 0-based ind and converted
 *      to the right data format)
 *
 * Stephen.Becker@Colorado.edu, 5/13/2016
 * */

#if defined(__GNUC__) && !(defined(__clang__)) && defined(NEEDS_UCHAR)
#include <uchar.h>
#endif
#include <math.h>
#include "mex.h"

#ifdef USE_BLAS
#include "blas.h"
#endif


#define XMATRIX 0
#define KSCALAR 1
#define INDICES 2
#define MEANS 0

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
     
{ 
    double *means, *X;
    mwSize p, n, i, j, k;
    mwSize *means_size;
    int *indices;
#ifdef USE_BLAS
    ptrdiff_t INCX = 1;
    ptrdiff_t pp;
    double alpha = 1.;
#endif
    if (nrhs == 0 ) {
        mexPrintf("Usage:  centers = computeMeans( X, k, int32(ind-1) );\n");
        return;
    }
    
    /* Check for proper number of arguments */
    if (nrhs != 3) {
        mexErrMsgIdAndTxt( "MATLAB:mexFile:invalidNumInputs",
                "Three input arguments required: X, K, ind.");
    } else if (nlhs > 1) {
        mexErrMsgIdAndTxt( "MATLAB:mexFile:maxlhs",
                "Too many output arguments: should be just MEANS.");
    }
    
    p  = mxGetM(prhs[XMATRIX]);
    n  = mxGetN(prhs[XMATRIX]);
    X  = mxGetPr(prhs[XMATRIX]);
#ifdef USE_BLAS
    pp = p;
#endif
    
    k  = (mwSize) mxGetScalar( prhs[KSCALAR] );
    
    if (mxGetM( prhs[INDICES] ) != 1 ){
        mexErrMsgIdAndTxt( "MATLAB:mexFile:badSize",
                "3rd input should be row vector.");
    }
    if (mxGetN( prhs[INDICES] ) != n ){
        mexErrMsgIdAndTxt( "MATLAB:mexFile:badSize",
                "3rd input should be row vector of length N.");
    }
 
    /* mxClassID mxGetClassID( prhs[INDICES] ); */
    if ( mxGetClassID( prhs[INDICES] ) != mxINT32_CLASS ) {
        mexPrintf("3rd input is of type %s\n", mxGetClassName( prhs[INDICES] ) );
        mexErrMsgIdAndTxt( "MATLAB:mexFile:badType",
                "3rd input should be of type INT32 (and 0-based).");
        
    }
    indices = mxGetData( prhs[INDICES] );
    
    plhs[MEANS] = mxCreateDoubleMatrix( p, k, mxREAL);
    means       = mxGetPr( plhs[MEANS] );
    means_size  = (mwSize*)calloc( k, sizeof( mwSize ) );
    for ( i=0 ; i<n; i++ ) {
        /* DAXPY ( N,  ALPHA, X, INCX, Y, INCY ) , Y <-- ALPHA*X+Y*/
#ifdef USE_BLAS
        daxpy( &pp, &alpha, X+i*p, &INCX, means+indices[i]*p, &INCX );
#else
        for (j=0;j<p;j++)
            means[ indices[i]*p + j ] += X[i*p+j];
#endif
        means_size[ indices[i] ] += 1;
    }
    for (i=0; i<k;i++ )
        for (j=0;j<p;j++)
            means[ i*p + j ] /= (double)means_size[i];

    
}
