%.o:%.c
	gcc -c -o test2.o test2.c
test2:test2.o
	gfortran -o $@ $< -llapacke -llapack -lf77blas -latlas

