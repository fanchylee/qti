#include <lapacke.h>
#include "ruby.h"
#include "complex.h"
#define min(a,b) ((a)>(b)?(b):(a))

extern int schmidt_decomposition(	int da, int db, lapack_complex_double *v,
									double *sv, lapack_complex_double *ia, lapack_complex_double *ib);

VALUE cQti ;
static VALUE sd(VALUE self, VALUE da_v, VALUE db_v, VALUE v_v){
	int da, db;
	long len;
	register int i;
	if(FIXNUM_P(da_v)&&FIXNUM_P(db_v)){
		da = FIX2INT(da_v);
		db = FIX2INT(db_v);
	}else{
		rb_raise(rb_eArgError, "The first two arguments must be integers");
	}
	Check_Type(v_v, T_ARRAY);
	len = RARRAY_LEN(v_v) ;
	if(len != da * db){
		rb_raise(rb_eArgError, "Length of the array must equal the product of the first two arguments");}
	complex v[da*db],ia[da*min(da,db)],ib[db*min(da,db)];
	double sv[min(da,db)];
	for(i=0;i<len;i++){
		VALUE ary_i = rb_ary_entry(v_v, i);
		switch(TYPE(ary_i)){
			case(T_COMPLEX):
			v[i]=NUM2DBL(RCOMPLEX(ary_i)->real) 
				+NUM2DBL(RCOMPLEX(ary_i)->imag)*I;
			break;
			case(T_FLOAT):
			case(T_FIXNUM):
			case(T_RATIONAL):
			v[i] = NUM2DBL(ary_i);
			break;
			default:
			rb_raise(rb_eArgError, "Contains elements of illegal type within the array");
			break;
		}
	}
	schmidt_decomposition(da,db,v,sv,ia,ib);
	return rb_float_new(sv[0]);
}
void Init_qti(){
	cQti = rb_define_class("Qti", rb_cObject);
	rb_define_singleton_method(cQti, "sd", sd, 3);
}

