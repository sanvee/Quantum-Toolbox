#ifndef __QUIT_TOOLBOX__
#define __QUIT_TOOLBOX__

#include <iostream>
#include <fstream>
#include <math.h>
#include <algorithm>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_vector.h>
#include <vector>
#include <gsl/gsl_complex.h>
#include <complex>
#include <gsl/gsl_complex_math.h>
#include "qr_complex.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#define HAVE_LAPACK_CONFIG_H
#define LAPACK_COMPLEX_CPP
#include <openblas/lapacke_config.h>
#include <openblas/lapacke_utils.h>
#include <openblas/lapacke.h>

extern std::complex< double > _0;
extern std::complex< double > _1;
extern std::complex< double > _neg_1;
extern std::complex< double > _i;
extern std::complex< double > _neg_i;

extern gsl_rng * _R_G;

gsl_rng * initialise_random_number_gernerators();

const double CUT_OFF_AS_ZERO = 1e-15;


bool compare_real_part ( std::complex< double > a,
                         std::complex< double > b );

void sort_vec_asc ( std::complex< double > * eigen,
                    int size );

void sort_vec_asc ( double * eigen,
                    int size );

void cutoff ( std::complex< double > * in );

void matrix_cutoff_neg ( std::complex< double > * in,
                         int size );

void vector_cutoff_neg ( std::complex< double > * in,
                         int size );

void vector_cutoff_neg ( double *in,
                         int size );


/*---------------------------------------------------------------------------------
 * GENERATORS (level 0)
 *------------------------------------------------------------------------------ -*/

void matrix_initialize_unity ( std::complex< double > * cmat,
                               int size );

void matrix_initialize_zero ( std::complex< double > * cmat,
                              int size );

void sphere_3_rand ( double vec[3] );

void sample_from_complex_sphere_n_rand ( std::complex< double > * vec, int size );

void sample_uniform_2_sphere_rad_sqrt_r ( gsl_complex * z );

void sample_uniform_2_unit_sphere ( double arr[2] );

void sample_random_ginibre_matrix ( std::complex< double > * cmat,
                                    int size );

void sample_random_ginibre_matrix ( double * mat,
                                    int n );

void sample_uniform_n_simplex ( std::complex< double >* carr,
                                int size );

void sample_uniform_n_simplex2 ( std::complex< double > * carr,
                                 int size );

//---------------------------------------------------------------------------------

void sample_uniform_hilbert_schmidt_density_matrix ( std::complex< double > * out,
        std::complex< double > * buffer,
        int size );

void sample_density_matrix_hs ( std::complex< double > * out,
                                int size );

void sample_density_matrix_up ( std::complex< double >* out,
                                std::complex< double >* buff,
                                std::complex< double >* tau,
                                int size );

void sample_density_matrix_up ( std::complex< double >* out,
                                int size );

void sample_separable_density_matrix_hs ( std::complex< double > * cmat,
        int * sub_spaces,
        int tot_nr_subspaces,
        int size );

void sample_separable_density_matrix_up ( std::complex< double > * cmat,
        int * sub_spaces,
        int tot_nr_subspaces,
        int size );

void sample_pure_density_matrix ( std::complex< double > * out,
                                  int size );

void sample_separable_pure_density_matrix ( std::complex< double > * cmat,
        int size,
        int * sub_spaces,
        int tot_nr_subspaces );


void sample_uniform_2x2_unitary_matrix ( std::complex< double > * cmat );


void sample_uniform_2x2_unitary_matrix ( gsl_matrix * mat );


int sample_unitary_matrix ( std::complex< double > * cmat,
                            std::complex< double > * tau_v,
                            std::complex< double > * temp_v,
                            int size );

int sample_unitary_matrix ( std::complex< double > * cmat,
                            int size );

void sample_diag_state ( std::complex< double > * cmat,
                         int size );



void sample_pure_n_qudit_ket ( gsl_vector * arr );

void sample_pure_n_qudit_ket ( gsl_vector_complex * carr );

template < typename T >
void vector_show (T * vec,int n );


void vector_write_to_file ( double * vec,
			    double res,
			    const char * s,
			    int n );

void vector_show ( std::complex< double > * cvec,
                   int n );

void matrix_show ( std::complex<double> * cmat,
                   int n );

void matrix_show_real_part ( std::complex<double> * cmat,
                             int n );


void matrix_show ( const double* mat, int n );
/*---------------------------------------------------------------------------------
 * CALCULATIONS (level 0)
 *------------------------------------------------------------------------------ -*/
void matrix_transpose ( double * in,
                        double * out,int dim );

void matrix_transpose ( std::complex< double > * in,
                        std::complex< double > * out,
                        int dim );

void matrix_transpose ( std::complex< double > * in,
                        int dim );

void matrix_hermitian_conjugate ( std::complex< double > * in,
                                  std::complex< double > * out
                                  ,int dim );
void matrix_hermitian_conjugate ( std::complex< double > * in ,
                                  int dim );

void transpose_block_mxm_at_offset_in_nxn ( std::complex< double > * in ,
        std::complex< double > * out,
        int dim,
        int offset_i,
        int offset_j,
        int dim_sub );

void transpose_block_mxm_at_offset_in_nxn ( double  * in,
        double  * out,
        int dim,
        int offset_i,
        int offset_j,
        int dim_sub );

void matrix_get_column ( int col, const std::complex< double >*const cmat, std::complex< double >* out, int size );
void matrix_get_partial_column (
    int col,
    int col_lenght,
    std::complex < double >*cmat,
    std::complex < double >*out,
    int size );
    void
matrix_set_partial_column (
    int col,
    int col_lenght,
    std::complex < double >*cvec,
    std::complex < double >*out,
    int size );


void matrix_partial_trace ( const std::complex< double >* in, std::complex< double >* out, int* trace_vector, int nr_subspaces, int size_tot );



void flip_block_kxk_in_block_mxm_in_nxn ( std::complex <double> * in,
        std::complex <double> *out,
        int dim_in, int offset_i,
        int offset_j,
        int k,
        int m );

void flip_block_kxk_in_block_mxm_in_nxn ( double * in,
        double *out,
        int dim_in,
        int offset_i,
        int offset_j,
        int k,
        int m );

void state_partial_transpose_nxn ( std::complex< double > * in,
                                   std::complex< double > * out,
                                   int dim_in ,
                                   int * transpose_vector,
                                   int tot_nr_subspaces );

void matrix_partial_cross_transpose ( std::complex< double > * in ,
                                      std::complex< double > * out,
                                      int dim_a,
                                      int dim_b,
                                      int size );

void matrix_trace_nxn ( std::complex< double > * in,
                        int dim_in,
                        std::complex< double > * result );

std::complex< double > matrix_trace_nxn ( std::complex< double > * in,
        int dim_in );

void vector_dyadic_nxn ( std::complex< double > * in1,
                         std::complex< double > * in2,
                         std::complex< double > * out,
                         int size );

void vector_dyadic_self ( std::complex< double > * in1,
                          std::complex< double > * out,
                          int size );
void vector_dyadic_self2 (
    std::complex < double >*in1,
    std::complex < double >*out,
    int size_v,
    int size_m);

void vector_dyadic_self ( double * in1,
                          double * out,
                          int size );

void vector_tensor_prod ( std::complex< double > * in1,
                          int dim_in1,
                          std::complex< double > * in2,
                          int dim_in2,
                          std::complex< double > * out,
                          int dim_out );

void matrix_log_nat_hermitian ( const std::complex< double > * in ,
                                std::complex< double > * in_copy ,
                                std::complex< double > * temp,
                                std::complex< double > * out,
                                double * eigen_temp,
				double * eigen_out,
                                int size );

void matrix_log_2_hermitian ( const std::complex< double >* in, std::complex< double >* in_copy, std::complex< double >* temp, std::complex< double >* out, double* eigen_temp, double* eigen_out, int size );

void matrix_log_nat_hermitian ( std::complex< double > * in ,
                                std::complex< double > * out,
                                int size );

void matrix_log_2_hermitian ( std::complex< double > * in ,
                              std::complex< double > * out,
                              int size );




void matrix_square_root ( std::complex< double > * in,
                          std::complex< double > * in_copy,
                          std::complex< double > * temp,
                          std::complex< double > * out,
                          double * eigen_dvec,
                          int size );

void matrix_square_root ( std::complex< double > * in,
                          std::complex< double > * out,
                          int size );

/*---------------------------------------------------------------------------------
 * CALCULATIONS (level 1)
 *------------------------------------------------------------------------------ -*/

void matrix_tensor_prod ( std::complex< double > * in1,
                          int dim_in1,
                          std::complex< double > * in2,
                          int dim_in2,
                          std::complex< double > * out,
                          int dim_out );

void matrix_tensor_prod ( double * in1,
                          int dim_in1,
                          double * in2,
                          int dim_in2,
                          double * out,
                          int dim_out );

void matrix_add ( std::complex< double > * cmat1,
                  std::complex< double > * cmat2,
                  int size );

void matrix_sub ( std::complex< double > * cmat1,
                  std::complex< double > * cmat2,
                  int size );

void matrix_add_to_first ( std::complex< double >* cmat1,
                           std::complex< double >* cmat2,
                           double factor,
                           int size );

void matrix_normalize ( std::complex< double > * cmat,
                        int size );

void vector_add ( std::complex< double > * cvec1,
                  double factor1,
                  std::complex< double > * cvec2,
                  double factor2,
                  std::complex< double > * out,
                  int size );


void vector_add_to_first ( std::complex< double > * cvec1,
                           std::complex< double > * cvec2,
                           double factor,
                           int size );

void matrix_mlt ( const std::complex< double > * cmat1,
                  const std::complex< double > * cmat2,
                  std::complex< double > * out,
                  int size );

void matrix_scalar_mult ( std::complex< double > * cmat,
                          double factor,
                          int size );

void matrix_scalar_mult ( std::complex< double > * cmat,
                          std::complex< double > factor,
                          int size );

void matrix_copy ( const std::complex< double >*const cmat,
		   std::complex< double >* copy,
		   int size );

template <typename T>
void
vector_copy ( const T * vector,
	      T * copy,
	      int size );

void
vector_scalar_mult (
    std::complex < double >*cvec,
    double factor,
    int size );
void
vector_scalar_mult (
    int * cvec,
    double factor,
    int size );

void vector_scalar_prod ( std::complex< double > * cvec1,
                          std::complex< double > * cvec2,
                          std::complex< double > * out,
                          int size );

std::complex< double >  vector_norm ( std::complex< double > * cvec1,
                                      int size );

void vector_norm ( std::complex< double > * cvec1,
                   std::complex< double >  * out,
                   int size );

std::complex< double > vector_coef_sum ( std::complex< double > * cvec,
        int size );

void matrix_vector_mult ( std::complex< double > * cmat,
                          std::complex< double > * cvec,
                          std::complex< double > * out,
                          int size );

/*---------------------------------------------------------------------------------
 * CALCULATIONS (level 2)
 *------------------------------------------------------------------------------ -*/
void matrix_eigenvalues ( const std::complex< double >*const cmat, std::complex< double >* eigen, int size );

void matrix_eigenvalues ( const std::complex< double >*const cmat, std::complex< double >* copy, std::complex< double >* eigen, int size );

/*---------------------------------------------------------------------------------
 * CALCULATIONS (level Quantum)
 *------------------------------------------------------------------------------ -*/

double state_hs_norm ( std::complex< double >* in,
                       int size );

double state_purity ( std::complex< double >* cmat,
                      int size );


void state_get_projector ( std::complex< double > * unitary,
                           int nr,
                           std::complex< double > * out,
                           int size );

double state_v_n_entropy ( const std::complex< double >* cmat, std::complex< double >* copy, std::complex< double >* eigen, int size );

double state_v_n_entropy ( const std::complex< double >*const cmat, int size );




double state_rel_entropy ( const std::complex< double >* in1,     // in1 log in1 - in1 log in2 
                           const std::complex< double >* in2,
                           int size );

double state_rel_entropy ( const std::complex< double >* in1,
			   const std::complex< double >* in2,
			   std::complex< double >* in1_c,
			   std::complex< double >* in2_c,
			   std::complex< double >* temp_cmat_1,
			   std::complex< double >* temp_cmat_2, 
			   std::complex< double >* temp_cmat_3,
			   std::complex< double >* cvec,
			   double* eigen_1,
			   double* eigen_2,
			   double* eigen_temp,
			   int size );

double state_negativity ( std::complex< double > * cmat,
                          std::complex< double > * cmat_buff,
                          std::complex< double > * buff,
                          int * subspaces,
                          int size );

double state_negativity ( std::complex< double > * cmat,
                          int * subspaces,
                          int size );

void state_decohere_on_subspace_ran ( std::complex< double >* in,
                                  std::complex< double >* out,
                                  int dim_in, int* dec_vector,
                                  int tot_nr_subspaces );

void state_decohere_on_subspace_ran ( std::complex< double >*const in,
				      std::complex< double >* out,
				      int dim_in,
				      int* dec_vector,
				      int tot_nr_subspaces,
				      std::complex< double >* temp1,
				      std::complex< double >* temp2,
				      std::complex< double >* total_op,
				      std::complex< double >* random_unitary,
				      std::complex< double >* basis_vector );

void state_decohere_on_subspace ( const std::complex< double >*const in,
				  std::complex< double >* out,
				  const std::complex< double > * const unitary,
				  int dim_in,
				  int* dec_vector,
				  int tot_nr_subspaces );

void state_decohere_on_subspace ( const std::complex< double >* const in,
				  std::complex< double >* out,
				  int dim_in,
				  int * dec_vector,
				  int tot_nr_subspaces,
				  std::complex< double >* temp1,
				  std::complex< double >* temp2,
				  std::complex< double >* total_op,
				  const std::complex< double > * const unitary,
				  std::complex< double >* basis_vector );


double state_quantum_commutance ( std::complex< double > * in,
                                  int dim_a,
                                  int dim_b,
                                  int dim_max,
                                  int dim_mat_max,
                                  int size,
                                  std::complex< double > * temp_cmat_1,
                                  std::complex< double > * temp_cmat_2,
                                  std::complex< double > * temp_cmat_3,
                                  double   * temp_dvec,
                                  std::complex< double > * temp_cvec,
                                  std::complex< double > * X );

double state_quantum_commutance ( std::complex< double > * in,
                                  int dim_a,
                                  int dim_b,
                                  int size );
double
state_Mutual_information_I_A_B ( const std::complex< double > *  const rho_AB,
                        int size_A,
                        int size_B,
                        int size_AB );





/*-----------------------------------------------------------------------------------
 * ----------------------------------------------------------------------------------
 */

#endif  // __QUIT_TOOLBOX__
