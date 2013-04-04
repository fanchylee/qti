qti.mexa64:qti.o ../schmidt_decomposition.o ../printmatrix.o
	mex -fortran $^ -llapacke -llapack -lf77blas -latlas -lm
qti.o:qti.c
	mex -c $< 
clean:
	rm qti.mexa64 qti.o
