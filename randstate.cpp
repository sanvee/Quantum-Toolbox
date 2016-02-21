#include <math.h>
#include <cstdio>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include "qr_complex.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>




void sample_random_normal_matrix_cplx(gsl_matrix_complex * cmat){
  size_t m = cmat->size1;
  size_t n = cmat->size2;
  
  for( int i=0 ;i<m;i++){
    for( int j=0 ;j<n;j++){
      gsl_matrix_complex_set(cmat,i,j,gsl_complex_rect(gsl_ran_gaussian(r_g,1),gsl_ran_gaussian(r_g,1)));
    }
  }
  
}
void sample_random_normal_matrix_real(gsl_matrix * mat){
  size_t m = mat->size1;
  size_t n = mat->size2;
  
  for( int i=0 ;i<m;i++){
    for( int j=0 ;j<n;j++){
      gsl_matrix_set(mat,i,j,gsl_ran_gaussian(r_g,1));
    }
  }
  
}
void sample_uniform_unitary_matrix_cplx(gsl_matrix_complex * cmat){

  gsl_vector_complex * tau = gsl_vector_complex_alloc(cmat->size1);
  gsl_matrix_complex * Q = gsl_matrix_complex_alloc(cmat->size1,cmat->size2);
  gsl_matrix_complex * R = gsl_matrix_complex_alloc(cmat->size1,cmat->size2);
  
  sample_random_normal_matrix_cplx(Q);
  gsl_linalg_complex_QR_decomp(Q, tau);
  gsl_linalg_complex_QR_unpack (Q,tau,cmat,R);
}
void sample_uniform_unitary_matrix_real(gsl_matrix* mat){

  gsl_vector * tau = gsl_vector_alloc(mat->size1);
  gsl_matrix * Q = gsl_matrix_alloc(mat->size1,mat->size2);
  gsl_matrix * R = gsl_matrix_alloc(mat->size1,mat->size2);
  
  sample_random_normal_matrix_real(Q);
  gsl_linalg_QR_decomp(Q, tau);
  gsl_linalg_QR_unpack (Q,tau,mat,R);
}
void sample_uniform_2_sphere_rad_sqrt_r(gsl_complex * z){ 
  
  /* This function samples Uniformly a complex number from a 2-sphere of radius sqrt(r), r stored in the real part
   *of the argument*/ 
  
  double r = sqrt(GSL_REAL (*z));
  
  GSL_SET_COMPLEX (z, gsl_ran_exponential(r_g,1), gsl_ran_exponential(r_g,1));
  
  r *= 1/sqrt(GSL_REAL (*z)*GSL_REAL (*z)+GSL_IMAG (*z)*GSL_IMAG (*z));
  
  GSL_SET_COMPLEX (z,r*GSL_REAL (*z),r*GSL_IMAG (*z));  
}
void sample_uniform_n_simplex_cplx (gsl_vector_complex * carr){
  
  gsl_vector_view real_part = gsl_vector_complex_real(carr); 
  gsl_vector * arr = &real_part.vector;
  size_t b = arr->stride;
  size_t n = arr->size;
  // we generate an Uniform Ordered array with Sum arr[i] = 1
  // 1) generate n exp distributed x
  for (int i=0;i<n;++i){
    arr->data[i*b]=gsl_ran_exponential(r_g,1);
  } 
   //2) implicit reorder by x[i]<---x[i-1]+(x[i]/total)
  double buff = gsl_blas_dasum(arr);
  arr->data[0]=(arr->data[0]/buff);
  for (int i=1;i<n-1;++i){  
    arr->data[i*b]=arr->data[(i-1)*b]+(arr->data[i*b]/buff);
  }
  arr->data[(n-1)*b]=1;
  //take the difference x[i]<---X[i]-x[i-1] ---> This is now an n Dimensional uniform distribution over the
  //n-simplex 
  buff = arr->data[0];
  for (int i=0;i<n-1;++i){  
    arr->data[i*b]=arr->data[(i+1)*b]-arr->data[i*b];
  }
  arr->data[(n-1)*b]=buff; 
}
void sample_uniform_n_simplex_real (gsl_vector * arr){
   
  size_t b = arr->stride;
  size_t n = arr->size;
  // we generate an Uniform Ordered array with Sum arr[i] = 1
  // 1) generate n exp distributed x
  for (int i=0;i<n;++i){
    arr->data[i*b]=gsl_ran_exponential(r_g,1);
  } 
   //2) implicit reorder by x[i]<---x[i-1]+(x[i]/total)
  double buff = gsl_blas_dasum(arr);
  arr->data[0]=(arr->data[0]/buff);
  for (int i=1;i<n-1;++i){  
    arr->data[i*b]=arr->data[(i-1)*b]+(arr->data[i*b]/buff);
  }
  arr->data[(n-1)*b]=1;
  //take the difference x[i]<---X[i]-x[i-1] ---> This is now an n Dimensional uniform distribution over the
  //n-simplex 
  buff = arr->data[0];
  for (int i=0;i<n-1;++i){  
    arr->data[i*b]=arr->data[(i+1)*b]-arr->data[i*b];
  }
  arr->data[(n-1)*b]=buff; 
}
void sample_pure_n_qudit_ket_real (gsl_vector * arr){
  
    sample_uniform_n_simplex_real (arr);
  size_t n = arr->size;
  size_t b= arr->stride;
  
  for (int i = 0;i < n; ++i){
    arr->data[i*b]=sqrt(arr->data[i*b]);
  }
}
void sample_pure_n_qudit_ket_cplx (gsl_vector_complex * carr){
    sample_uniform_n_simplex_cplx (carr);
  size_t n = carr->size;
  for ( int i = 0 ; i < n ; ++i){
  sample_uniform_2_sphere_rad_sqrt_r(gsl_vector_complex_ptr(carr,i));
  }
}
void sample_uniform_2x2_unitary_matrix_cplx(gsl_matrix_complex * cmat){
  
  
  gsl_vector_complex_view row1 = gsl_matrix_complex_row(cmat,0);
  gsl_vector_complex * carr = &row1.vector;
  sample_pure_n_qudit_ket_cplx(carr);
  gsl_matrix_complex_set (cmat,1,0,gsl_complex_negative (gsl_complex_conjugate(gsl_matrix_complex_get(cmat,0,1))));
  gsl_matrix_complex_set (cmat,1,1,                     (gsl_complex_conjugate(gsl_matrix_complex_get(cmat,0,0))));   
}
void sample_uniform_2x2_unitary_matrix_real(gsl_matrix * mat){
  /* Samples a and b from unit circle so that a^2+b^2=1*/
  mat->data[0] = gsl_ran_gaussian(r_g,1);
  mat->data[1] = gsl_ran_gaussian(r_g,1);
  double r = 1/sqrt(mat->data[0]* mat->data[0]+mat->data[1]* mat->data[1]);
  mat->data[0]*= r;
  mat->data[1]*= r;
  
  mat->data[mat->tda]=-mat->data[1];
  mat->data[mat->tda+1]=mat->data[0];
}


