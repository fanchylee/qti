#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <lapacke.h>
#include <math.h>
#include <complex.h>

/*
 *compute the Schmidt Decomposition
 *only for real matrices
 */
int sd_real(
	gsl_vector* v,
	/*the amplitute of every states, 
	 *index is the just the state as binary value.
	 *For example, 
	 *		.8*|110> + .6*|001> 
	 *is represented as 
	 *		[0  , 0.6, 0  , 0  , 0  , 0  , 0.8, 0 ]
	 *		---------------------------------------
	 *		000B 001B 010B 011B 100B 101B 110B 111B
	 */
	gsl_vector* sv,
	/*the coefficients of the schmidt decompsition
	 */
	gsl_matrix* ia,
	/*the first Schmidt base 
	 */
	gsl_matrix* ib){
	/*the second Schmidt base 
	 */
	register size_t i,j;
	size_t ia_size = v->size >> 1;
	size_t ib_size = 2 ;
	gsl_matrix* A = gsl_matrix_alloc (ia_size, ib_size);
	gsl_matrix* X = gsl_matrix_alloc (ib_size, ib_size);
	gsl_vector* work = gsl_vector_alloc(ib_size);
	for(i=0; i<ia_size; i++){
		for(j=0; j<ib_size; j++){
			 gsl_matrix_set(A, i, j, gsl_vector_get(v, (i<<1) + j));}}
	gsl_linalg_SV_decomp_mod(A, X, ib, sv, work);
	gsl_matrix_memcpy(ia, A);
	gsl_matrix_transpose(ib);
	gsl_matrix_free (A);
	gsl_matrix_free (X);
	gsl_vector_free (work);
	return 0;
}
long double trace_complex(gsl_matrix_complex* v){
	long double tr = 0;
	register size_t i;
	if(v->size1 != v->size2){
		fprintf(stderr, "not a square matrix");
		exit(EXIT_FAILURE);
	}
	for(i = 0;i < v->size1;i++){
		gsl_complex_add(tr, gsl_matrix_complex_get(v,i,i));
	}
	return tr;
}
int partial_trace_complex(gsl_matrix_complex* v,
	gsl_matrix_complex* r, 
	/*
	 *the result will be stored in r
	 */
	size_t d){
	/*
	 *d is the dimension of the matrix to be traced out
	 */
	size_t step = v->size1 / d ;
	register size_t i = 0, j = 0, k = 0;
	for(i=0;i<step;i++){
		for(j=0;j<step;j++){
			size_t initM = i*step;
			size_t initN = j*step;
			gsl_complex tr;
			GSL_SET_COMPLEX(&tr, 0, 0);
			for(k=0;k<d;k++){
				gsl_complex_add(gsl_matrix_complex_get(v, initM+k, initN+k), tr);
			}
			gsl_matrix_complex_set(r, i, j, tr);
		}
	}
	return 0;
}
int schmidt_decomposition(
	int da, int db,
	/*The dimension of the first and the second
	 *quantum state.
	 */
	double* v,
	/*the amplitute of every states, 
	 *index is the just the state as binary value.
	 *For example, 
	 *		.8*|01> + .6*|11> 
	 *is represented as 
	 *		[0 ,0.8, 0 ,0.6]
	 *		----------------
	 *		00B 01B 10B 11B
	 */
	double* sv,
	/*the coefficients of the schmidt decompsition
	 */
	double* ia,
	/*the first Schmidt base 
	 */
	double* ib){
	/*the second Schmidt base 
	 */
	register int i,j;
	double tempib[db*min(da,db)] ;
	info = LAPACKE_zgesvd( LAPACK_ROW_MAJOR, 'A', 'A', da, db, v, db, sv, ia, min(da,db), tempib, db, superb2);
	/*ib = transpose of tempib
	 */
	for(i=0;i<da;i++){
		for(j=0;j<db;j++){
			ib[j][i] = tempib[i][j];}}
	return 0;
}

int main(){
	gsl_vector* v = gsl_vector_alloc(4);
	gsl_vector* sv = gsl_vector_alloc(2);
	gsl_matrix* ib = gsl_matrix_alloc(2,2);
	gsl_matrix* ia = gsl_matrix_alloc(2,2);
	double complex a[4*3] = {
		 5.91 +  -5.69i,  7.09 +   2.72i,  7.78 +  -4.06i, -0.79 +  -7.21i,
		-3.15 +  -4.08i, -1.89 +   3.27i,  4.57 +  -2.07i, -3.88 +  -3.30i,
		-4.89 +   4.20i,  4.10 +  -6.70i,  3.28 +  -3.84i,  3.84 +   1.19i
	};
	double s[3];
	double superb[3-1];
	lapack_complex_double u[3*3], vt[4*3];
	lapack_int info;
	gsl_matrix_complex* cm = gsl_matrix_complex_calloc(2,2);
	gsl_matrix_complex_set_identity(cm);
	info = LAPACKE_zgesvd(LAPACK_ROW_MAJOR, 'A', 'A', 3, 4, a, 4, s, u, 3, vt, 4, superb);
	if( info > 0 ) {
		printf( "The algorithm computing SVD failed to converge.\n" );
		exit( 1 );
	}
	gsl_vector_set(v, 0, 1/sqrt(3));
	gsl_vector_set(v, 1, 1/sqrt(3));
	gsl_vector_set(v, 2, 1/sqrt(3));
	gsl_vector_set(v, 3, 0);
	sd_real(v,sv,ia,ib);
	printf("ia printed\n");
	gsl_matrix_fprintf(stdout, 	ia, "%g");
	printf("ib printed\n");
	gsl_matrix_fprintf(stdout, 	ib, "%g");
	printf("sv printed\n");
	gsl_vector_fprintf(stdout,  sv, "%g");
	printf("cm printed\n");
	
	gsl_matrix_complex_fprintf(stdout, cm, "%g");
	

	gsl_vector_free(v );
	gsl_vector_free(sv);
	gsl_matrix_free(ib);
	gsl_matrix_free(ia);
	gsl_matrix_complex_free(cm);
}
