#include <lapacke.h>
#include <math.h>

int main(){
	lapack_complex_double	v[2*2] = {
		1/sqrt(3), 1/sqrt(3), 0, 1/sqrt(3)}, ia[2*2], ib[2*2],
							v1[2*2] = {
		0.5,0.5,0.5,0.5};
	double sv[2];

	schmidt_decomposition(2,2,v1,sv,ia,ib);
	print_rmatrix( "Singular values", 1, 2, sv, 1 );
	print_matrix( "Left singular vectors (stored columnwise)", 2, 2, ia, 2 );
	print_matrix( "Right singular vectors (stored rowwise)", 2, 2, ib, 2 );
}
