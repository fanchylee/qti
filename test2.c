/*******************************************************************************
*  Copyright (C) 2009-2012 Intel Corporation. All Rights Reserved.
*  The information and material ("Material") provided below is owned by Intel
*  Corporation or its suppliers or licensors, and title to such Material remains
*  with Intel Corporation or its suppliers or licensors. The Material contains
*  proprietary information of Intel or its suppliers and licensors. The Material
*  is protected by worldwide copyright laws and treaty provisions. No part of
*  the Material may be copied, reproduced, published, uploaded, posted,
*  transmitted, or distributed in any way without Intel's prior express written
*  permission. No license under any patent, copyright or other intellectual
*  property rights in the Material is granted to or conferred upon you, either
*  expressly, by implication, inducement, estoppel or otherwise. Any license
*  under such intellectual property rights must be express and approved by Intel
*  in writing.
*
********************************************************************************
*/
/*
   LAPACKE_zgesvd Example.
   =======================

   Program computes the singular value decomposition of a general
   rectangular complex matrix A:

   (  5.91, -5.69) (  7.09,  2.72) (  7.78, -4.06) ( -0.79, -7.21)
   ( -3.15, -4.08) ( -1.89,  3.27) (  4.57, -2.07) ( -3.88, -3.30)
   ( -4.89,  4.20) (  4.10, -6.70) (  3.28, -3.84) (  3.84,  1.19)

   Description.
   ============

   The routine computes the singular value decomposition (SVD) of a complex
   m-by-n matrix A, optionally computing the left and/or right singular
   vectors. The SVD is written as

   A = U*SIGMA*VH

   where SIGMA is an m-by-n matrix which is zero except for its min(m,n)
   diagonal elements, U is an m-by-m unitary matrix and VH (V conjugate
   transposed) is an n-by-n unitary matrix. The diagonal elements of SIGMA
   are the singular values of A; they are real and non-negative, and are
   returned in descending order. The first min(m, n) columns of U and V are
   the left and right singular vectors of A.

   Note that the routine returns VH, not V.

   Example Program Results.
   ========================

 LAPACKE_zgesvd (row-major, high-level) Example Program Results

 Singular values
  17.63  11.61   6.78

 Left singular vectors (stored columnwise)
 ( -0.86,  0.00) (  0.40,  0.00) (  0.32,  0.00)
 ( -0.35,  0.13) ( -0.24, -0.21) ( -0.63,  0.60)
 (  0.15,  0.32) (  0.61,  0.61) ( -0.36,  0.10)

 Right singular vectors (stored rowwise)
 ( -0.22,  0.51) ( -0.37, -0.32) ( -0.53,  0.11) (  0.15,  0.38)
 (  0.31,  0.31) (  0.09, -0.57) (  0.18, -0.39) (  0.38, -0.39)
 (  0.53,  0.24) (  0.49,  0.28) ( -0.47, -0.25) ( -0.15,  0.19)
*/
#include <stdlib.h>
#include <stdio.h>
#include <lapacke.h>
#include <math.h>


#define min(a,b) ((a)>(b)?(b):(a))

/* Auxiliary routines prototypes */
extern void print_matrix( char* desc, lapack_int m, lapack_int n, lapack_complex_double* a, lapack_int lda );
extern void print_rmatrix( char* desc, lapack_int m, lapack_int n, double* a, lapack_int lda );

/* Parameters */
#define M 3
#define N 4
#define LDA N
#define LDU M
#define LDVT N

/* Main program */
int main() {
        /* Locals */
		lapack_complex_double v[2*2] = {
			1/sqrt(3), 1/sqrt(3), 0, 1/sqrt(3)};
		double sv[2];
		double superb2[1];
		lapack_complex_double ia[2*2], ib[2*2];
        lapack_int m = M, n = N, lda = LDA, ldu = LDU, ldvt = LDVT, info;
        /* Local arrays */
        double s[M];
        double superb[min(M,N)-1];
        lapack_complex_double u[LDU*M], vt[LDVT*N];
        lapack_complex_double a[LDA*M] = {
			 5.91 +  -5.69i,  7.09 +   2.72i,  7.78 +  -4.06i, -0.79 +  -7.21i,
			-3.15 +  -4.08i, -1.89 +   3.27i,  4.57 +  -2.07i, -3.88 +  -3.30i,
			-4.89 +   4.20i,  4.10 +  -6.70i,  3.28 +  -3.84i,  3.84 +   1.19i
        };
        /* Executable statements */
        printf( "LAPACKE_zgesvd (row-major, high-level) Example Program Results\n" );
        /* Compute SVD */
        info = LAPACKE_zgesvd( LAPACK_ROW_MAJOR, 'A', 'A', m, n, a, lda, s,
         u, ldu, vt, ldvt, superb );
        /* Check for convergence */
        if( info > 0 ) {
                printf( "The algorithm computing SVD failed to converge.\n" );
                exit( 1 );
        }
		info = LAPACKE_zgesvd( LAPACK_ROW_MAJOR, 'A', 'A', 2, 2, v, 2, sv, ia, 2, ib, 2, superb2);
        /* Print singular values */
        print_rmatrix( "Singular values", 1, m, s, 1 );
        /* Print left singular vectors */
        print_matrix( "Left singular vectors (stored columnwise)", m, m, u, ldu );
        /* Print right singular vectors */
        print_matrix( "Right singular vectors (stored rowwise)", m, n, vt, ldvt );

        print_rmatrix( "Singular values", 1, 2, sv, 1 );
        print_matrix( "Left singular vectors (stored columnwise)", 2, 2, ia, 2 );
        print_matrix( "Right singular vectors (stored rowwise)", 2, 2, ib, 2 );
        exit( 0 );
} /* End of LAPACKE_zgesvd Example */

/* Auxiliary routine: printing a matrix */
void print_matrix( char* desc, lapack_int m, lapack_int n, lapack_complex_double* a, lapack_int lda ) {
        lapack_int i, j;
        printf( "\n %s\n", desc );
        for( i = 0; i < m; i++ ) {
                for( j = 0; j < n; j++ )
                        printf( " (%6.2f,%6.2f)", creal(a[i*lda+j]), cimag(a[i*lda+j]));
                printf( "\n" );
        }
}

/* Auxiliary routine: printing a real matrix */
void print_rmatrix( char* desc, lapack_int m, lapack_int n, double* a, lapack_int lda ) {
        lapack_int i, j;
        printf( "\n %s\n", desc );
        for( i = 0; i < m; i++ ) {
                for( j = 0; j < n; j++ ) printf( " %6.2f", a[i*lda+j] );
                printf( "\n" );
        }
}
