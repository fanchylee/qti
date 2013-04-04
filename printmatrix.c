#include <stdio.h>
#include <lapacke.h>
#include <math.h>
#include <complex.h>
/* Auxiliary routine: printing a matrix */
void print_matrix( char* desc, lapack_int m, lapack_int n, lapack_complex_double* a, lapack_int lda ) { 
        lapack_int i, j;
        printf( "\n %s\n", desc );
        for( i = 0; i < m; i++ ) { 
                for( j = 0; j < n; j++ )
                        printf( " (%6.4f,%6.4f)", creal(a[i*lda+j]), cimag(a[i*lda+j]));
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

