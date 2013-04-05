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
	register int i,j;
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
	lapack_complex_double v[da*db],ia[da*min(da,db)],ib[db*min(da,db)];
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
	VALUE retary, rethash[min(da,db)], va[min(da,db)], vb[min(da,db)];
	for(j=0;j<min(da,db);j++){
		rethash[j] = rb_hash_new();
		va[j] = rb_ary_new2(da);
		for(i=0;i<da;i++){
			rb_ary_store(va[j], i, rb_Complex(rb_float_new(creal(ia[i*da+j])), rb_float_new(cimag(ia[i*da+j]))));}
		rb_hash_aset(rethash[j], ID2SYM(rb_intern("ia")), va[j]);
		vb[j] = rb_ary_new2(db);
		for(i=0;i<db;i++){
			rb_ary_store(vb[j], i, rb_Complex(rb_float_new(creal(ib[i*db+j])), rb_float_new(cimag(ib[i*db+j]))));}
		rb_hash_aset(rethash[j], ID2SYM(rb_intern("ib")), vb[j]);
		rb_hash_aset(rethash[j], ID2SYM(rb_intern("sv")), rb_float_new(sv[j]));
	}
	retary = rb_ary_new4(min(da,db), rethash);
	
	return retary;
}
void Init_qti(){
	cQti = rb_define_class("Qti", rb_cObject);
	rb_define_singleton_method(cQti, "sd", sd, 3);
}

