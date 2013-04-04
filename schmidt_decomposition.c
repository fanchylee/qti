#include <lapacke.h>
#include <math.h>
#include <complex.h>
#define min(a,b) ((a)>(b)?(b):(a))


int schmidt_decomposition(
	int da, int db,
	/*The dimension of the first and the second
	 *quantum state.
	 */
	lapack_complex_double* v,
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
	lapack_complex_double* ia,
	/*the first Schmidt base 
	 */
	lapack_complex_double* ib){
	/*the second Schmidt base 
	 */
	register int i,j;
	lapack_complex_double ibt[db*min(da,db)];
	double superb[min(da,db)-1];
	if(LAPACKE_zgesvd(LAPACK_ROW_MAJOR, 'A', 'A', da, db, v, db, sv, 
						ia, min(da,db), ibt, db, superb) > 0){
		return(-1);
	}
	/*ib = conjugate transpose of ibt
	 */
	for(i=0;i<da;i++){
		for(j=0;j<db;j++){
			ib[da*j+i] = conj(ibt[db*i+j]);}}
	return 0;
}

