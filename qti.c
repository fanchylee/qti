#include <stdio.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <math.h>

/*compute the Schmidt Decomposition
 *recursively
 */
int sd(
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
	register int i,j;
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
int main(){
	gsl_vector* v = gsl_vector_alloc(4);
	gsl_vector* sv = gsl_vector_alloc(2);
	gsl_matrix* ib = gsl_matrix_alloc(2,2);
	gsl_matrix* ia = gsl_matrix_alloc(2,2);
	gsl_vector_set(v, 0, 1/sqrt(3));
	gsl_vector_set(v, 1, 1/sqrt(3));
	gsl_vector_set(v, 2, 1/sqrt(3));
	gsl_vector_set(v, 3, 0);
	sd(v,sv,ia,ib);
	printf("ia printed\n");
	gsl_matrix_fprintf(stdout, 	ia, "%g");
	printf("ib printed\n");
	gsl_matrix_fprintf(stdout, 	ib, "%g");
	printf("sv printed\n");
	gsl_vector_fprintf(stdout,  sv, "%g");
	gsl_vector_free(v );
	gsl_vector_free(sv);
	gsl_matrix_free(ib);
	gsl_matrix_free(ia);
	
}
