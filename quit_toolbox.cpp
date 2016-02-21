#include "quit_toolbox.h"

std::complex< double > _0 = std::complex< double > ( 0.0,0.0 );
std::complex< double > _1 = std::complex< double > ( 1.0,0.0 );
std::complex< double > _neg_1 = std::complex< double > ( -1.0,0.0 );
std::complex< double > _i = std::complex< double > ( 0.0,1.0 );
std::complex< double > _neg_i = std::complex< double > ( 0.0,-1.0 );


gsl_rng * _R_G = initialise_random_number_gernerators();



gsl_rng * initialise_random_number_gernerators (
)
    {

    /* Declare generator and generator_type variables */
    const gsl_rng_type *T;
    /* Setup variables for the generator */
    gsl_rng_env_setup ();
    /* Specify the type of generator we will use */
    T = gsl_rng_taus;
    /* Getting a seed from the Linux System */
    std::ifstream urandom_FILE;
    urandom_FILE.open ( "/dev/urandom", std::ios::in | std::ios::binary );
    if ( urandom_FILE.is_open () )
        {

        urandom_FILE.read ( reinterpret_cast < char *> ( &gsl_rng_default_seed ),
                            sizeof ( gsl_rng_default_seed ) );

        }

    urandom_FILE.close ();
    /* Creating the number generator with specified options */
    _R_G = gsl_rng_alloc ( T );

    return _R_G;

    }



void
matrix_eigenvalues (
    const std::complex < double > * const cmat,
    std::complex < double >*copy,
    std::complex < double >*eigen,
    int size )
    {

    matrix_copy ( cmat, copy, size );

    LAPACKE_zgehrd ( LAPACK_ROW_MAJOR, size, 1, size, copy, size, eigen );
    LAPACKE_zhseqr ( LAPACK_ROW_MAJOR, 'E', 'N', size, 1, size, copy, size,
                     eigen, 0, size );

    }

void
matrix_eigenvalues (
    const std::complex < double > * const cmat,
    std::complex < double >*eigen,
    int size )
    {

    std::complex < double >*copy;

    copy = new std::complex < double >[size * size];

    matrix_copy ( cmat, copy, size );

    LAPACKE_zgehrd ( LAPACK_ROW_MAJOR, size, 1, size, copy, size, eigen );
    LAPACKE_zhseqr ( LAPACK_ROW_MAJOR, 'E', 'N', size, 1, size, copy, size,
                     eigen, 0, size );

    delete[]copy;

    }

void
sphere_3_rand (
    double vec[3] )
    {

    vec[0] = gsl_ran_gaussian ( _R_G, 1 );
    vec[1] = gsl_ran_gaussian ( _R_G, 1 );
    vec[2] = gsl_ran_gaussian ( _R_G, 1 );
    double r = 1 / sqrt ( vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2] );
    vec[0] = r * vec[0];
    vec[1] = r * vec[1];
    vec[2] = r * vec[2];

    }


void
sample_from_complex_sphere_n_rand ( std::complex< double > * vec, int size )
    {
    double norm ( 0 );

    for ( int i = 0 ; i < size ; ++i )
        {
        vec[i] = std::complex< double > ( gsl_ran_gaussian_ziggurat ( _R_G,1 ),gsl_ran_gaussian_ziggurat ( _R_G,1 ) );
        norm += vec[i].real() * vec[i].real() + vec[i].imag() *vec[i].imag();
        }

    norm = 1.0 / ( sqrt ( norm ) );
    for ( int i = 0 ; i < size ; ++i )
        {
        vec[i] = norm * vec[i];
        }

    }


void
sample_random_ginibre_matrix ( std::complex < double > * cmat,int size )

    {

    for ( int i = 0; i < size; ++i )
        {

        for ( int j = 0; j < size; ++j )
            {

            cmat[j + size * i] = std::complex < double > ( gsl_ran_gaussian_ziggurat ( _R_G,1 ), gsl_ran_gaussian_ziggurat ( _R_G,1 ) );

            }

        }

    }

void
sample_density_matrix_hs (
    std::complex < double >*out,
    int size )
    {

    std::complex < double >*buffer;
    buffer = new std::complex < double >[size * size];

    sample_random_ginibre_matrix ( buffer, size );

    for ( int i = 0; i < size; ++i )
        {

        for ( int j = 0; j < size; ++j )
            {

            out[j + size * i] = std::complex < double > (
                                    0,
                                    0 );

            for ( int k = 0; k < size; ++k )
                {

                out[j + size * i] +=
                    buffer[k + i * size] * std::conj ( buffer[k + j * size] );

                }

            }

        }

    buffer[0] = matrix_trace_nxn ( out, size );
    buffer[0] = 1.0 / buffer[0];
    matrix_scalar_mult ( out, buffer[0], size );

    delete[]buffer;

    }

void
sample_density_matrix_up (
    std::complex < double >*out,
    std::complex < double >*buff,
    std::complex < double >*tau,
    int size )
    {

    sample_unitary_matrix ( buff, tau, out, size );

    sample_uniform_n_simplex ( tau, size );

    for ( int i = 0; i < size; ++i )
        {

        for ( int j = 0; j < size; ++j )
            {

            out[j + size * i] = std::complex < double > (
                                    0,
                                    0 );

            for ( int k = 0; k < size; ++k )
                {

                out[j + size * i] +=
                    buff[k + i * size] * tau[k] * std::conj ( buff[k + j * size] );
                }

            }

        }

    }

void
sample_density_matrix_up (
    std::complex < double >*out,
    int size )
    {

    std::complex < double >*buff;
    buff = new std::complex < double >[size * size];

    std::complex < double >*tau;
    tau = new std::complex < double >[size];

    sample_density_matrix_up ( out, buff, tau, size );


    delete[]buff;
    delete[]tau;

    }

int
sample_unitary_matrix (
    std::complex < double >*cmat,
    std::complex < double >*tau_v,
    std::complex < double >*temp_v,
    int size )
    {

    sample_random_ginibre_matrix ( cmat, size );

    int info1 =
        LAPACKE_zgeqrf ( LAPACK_ROW_MAJOR, size, size, cmat, size, tau_v );

    for ( int i = 0; i < size; ++i )
        {

        temp_v[i] = cmat[i * ( size + 1 )] / abs ( cmat[i * ( size + 1 )] );

        }

    int info2 =
        LAPACKE_zungqr ( LAPACK_ROW_MAJOR, size, size, size, cmat, size, tau_v );

    for ( int i = 0; i < size; ++i )
        {

        for ( int j = 0; j < size; ++j )
            {

            cmat[i * size + j] = cmat[i * size + j] * temp_v[i];

            }

        }

    return info1 + info2;

    }

int
sample_unitary_matrix (
    std::complex < double >*cmat,
    int size )
    {

    std::complex < double >*tau;
    tau = new std::complex < double >[size];

    std::complex < double >*temp;
    temp = new std::complex < double >[size];

    sample_random_ginibre_matrix ( cmat, size );

    int info1 = LAPACKE_zgeqrf ( LAPACK_ROW_MAJOR, size, size, cmat, size, tau );

    for ( int i = 0; i < size; ++i )
        {

        temp[i] = cmat[i * ( size + 1 )] / abs ( cmat[i * ( size + 1 )] );

        }

    int info2 = LAPACKE_zungqr ( LAPACK_ROW_MAJOR, size, size, size, cmat, size, tau );

    for ( int i = 0; i < size; ++i )
        {

        for ( int j = 0; j < size; ++j )
            {

            cmat[i * size + j] = cmat[i * size + j] * temp[i];

            }

        }

    delete[]tau;
    delete[]temp;

    return info1 + info2;

    }

void
sample_separable_density_matrix_hs (
    std::complex < double >*cmat,
    int *sub_spaces,
    int tot_nr_subspaces,
    int size )
    {

    int dim_check = 1;
    int curr_dim = 1;
    int max_dim = 1;

    for ( int i = 0; i < tot_nr_subspaces; ++i )
        {
        dim_check *= sub_spaces[i];
        if ( max_dim < sub_spaces[i] )
            {
            max_dim = sub_spaces[i];
            }
        }

    if ( dim_check != size )
        {
        std::cout <<
                  "ERROR IN FUNCTION sample_separable_density_matrix: INCOMPATIBLE DIMENSIONS";
        return;
        }

    std::complex < double >*local_unitary,
        *local_unitary_hc, *local_density, *temp, *temp_global;

    temp_global = new std::complex < double >[dim_check * dim_check];
    temp_global[0] = std::complex < double > ( 1,0 );

    local_unitary = new std::complex < double >[max_dim * max_dim];
    local_unitary_hc = new std::complex < double >[max_dim * max_dim];
    local_density = new std::complex < double >[max_dim * max_dim];
    temp = new std::complex < double >[max_dim * max_dim];

    for ( int i = tot_nr_subspaces - 1; i >= 0; --i )
        {

        sample_density_matrix_hs ( local_density, sub_spaces[i] );
        sample_unitary_matrix ( local_unitary, sub_spaces[i] );
        matrix_hermitian_conjugate ( local_unitary, local_unitary_hc,sub_spaces[i] );
        matrix_mlt ( local_density, local_unitary_hc, temp, sub_spaces[i] );
        matrix_mlt ( local_unitary, temp, local_density, sub_spaces[i] );
        matrix_tensor_prod ( local_density, sub_spaces[i], temp_global, curr_dim,cmat, curr_dim * sub_spaces[i] );

        curr_dim *= sub_spaces[i];
        matrix_copy ( cmat, temp_global, curr_dim );

        }

    delete[]temp_global;
    delete[]local_unitary;
    delete[]local_unitary_hc;
    delete[]local_density;
    delete[]temp;

    }

void
sample_separable_density_matrix_up (
    std::complex < double >*cmat,
    int *sub_spaces,
    int tot_nr_subspaces,
    int size )
    {

    int dim_check = 1;
    int curr_dim = 1;
    int max_dim = 1;

    for ( int i = 0; i < tot_nr_subspaces; ++i )
        {
        dim_check *= sub_spaces[i];
        if ( max_dim < sub_spaces[i] )
            {
            max_dim = sub_spaces[i];
            }
        }

    if ( dim_check != size )
        {
        std::cout <<
                  "ERROR IN FUNCTION sample_separable_density_matrix: INCOMPATIBLE DIMENSIONS";
        return;
        }

    std::complex < double >*local_unitary,
        *local_unitary_hc, *local_density, *temp, *temp_global;

    temp_global = new std::complex < double >[dim_check * dim_check];
    temp_global[0] = std::complex < double > (
                         1,
                         0 );

    local_unitary = new std::complex < double >[max_dim * max_dim];
    local_unitary_hc = new std::complex < double >[max_dim * max_dim];
    local_density = new std::complex < double >[max_dim * max_dim];
    temp = new std::complex < double >[max_dim * max_dim];

    for ( int i = tot_nr_subspaces - 1; i >= 0; --i )
        {

        sample_density_matrix_up ( local_density, sub_spaces[i] );
        sample_unitary_matrix ( local_unitary, sub_spaces[i] );

        matrix_hermitian_conjugate ( local_unitary, local_unitary_hc, sub_spaces[i] );

        matrix_mlt ( local_density, local_unitary_hc, temp, sub_spaces[i] );
        matrix_mlt ( local_unitary, temp, local_density, sub_spaces[i] );
        matrix_tensor_prod ( local_density, sub_spaces[i], temp_global, curr_dim, cmat, curr_dim * sub_spaces[i] );

        curr_dim *= sub_spaces[i];
        matrix_copy ( cmat, temp_global, curr_dim );

        }

    delete[]temp_global;
    delete[]local_unitary;
    delete[]local_unitary_hc;
    delete[]local_density;
    delete[]temp;

    }

void
sample_separable_pure_density_matrix (
    std::complex < double >*cmat,
    int size,
    int *sub_spaces,
    int tot_nr_subspaces )
    {

    int dim_check = 1;
    int curr_dim = 1;
    int max_dim = 1;

    for ( int i = 0; i < tot_nr_subspaces; ++i )
        {
        dim_check *= sub_spaces[i];
        if ( max_dim < sub_spaces[i] )
            {
            max_dim = sub_spaces[i];
            }
        }

    if ( dim_check != size )
        {
        std::cout <<
                  "ERROR IN FUNCTION sample_separable_density_matrix: INCOMPATIBLE DIMENSIONS";
        return;
        }

    std::complex < double >*local_unitary,
        *local_unitary_hc, *local_density, *temp, *temp_global;

    temp_global = new std::complex < double >[dim_check * dim_check];
    temp_global[0] = std::complex < double > (
                         1,
                         0 );

    local_unitary = new std::complex < double >[max_dim * max_dim];
    local_unitary_hc = new std::complex < double >[max_dim * max_dim];
    local_density = new std::complex < double >[max_dim * max_dim];
    temp = new std::complex < double >[max_dim * max_dim];

    for ( int i = tot_nr_subspaces - 1; i >= 0; --i )
        {

        sample_pure_density_matrix ( local_density, sub_spaces[i] );
        sample_unitary_matrix ( local_unitary, sub_spaces[i] );

        matrix_hermitian_conjugate ( local_unitary, local_unitary_hc, sub_spaces[i] );

        matrix_mlt ( local_density, local_unitary_hc, temp, sub_spaces[i] );
        matrix_mlt ( local_unitary, temp, local_density, sub_spaces[i] );
        matrix_tensor_prod ( local_density, sub_spaces[i], temp_global, curr_dim, cmat, curr_dim * sub_spaces[i] );

        curr_dim *= sub_spaces[i];
        matrix_copy ( cmat, temp_global, curr_dim );

        }

    delete[]temp_global;
    delete[]local_unitary;
    delete[]local_unitary_hc;
    delete[]local_density;
    delete[]temp;

    }

void
sample_pure_density_matrix (
    std::complex < double >*out,
    int size )
    {

    std::complex < double >*cvec;
    cvec = new std::complex < double >[size];

    sample_unitary_matrix ( out, size );
    matrix_get_column ( 1, out, cvec, size );
    vector_dyadic_self ( cvec, out, size );

    delete[]cvec;

    }

void
sample_uniform_2_sphere_rad_sqrt_r (
    gsl_complex * z )
    {

    /* This function samples Uniformly a complex number from a 2-sphere of radius sqrt(r), r stored in the real part
     *of the argument*/
    double r = sqrt ( GSL_REAL ( *z ) );
    GSL_SET_COMPLEX ( z, gsl_ran_exponential ( _R_G, 1 ),
                      gsl_ran_exponential ( _R_G, 1 ) );
    r *=1 / sqrt ( GSL_REAL ( *z ) * GSL_REAL ( *z ) + GSL_IMAG ( *z ) * GSL_IMAG ( *z ) );

    GSL_SET_COMPLEX ( z, r * GSL_REAL ( *z ), r * GSL_IMAG ( *z ) );

    }

void
sample_uniform_n_simplex (
    std::complex < double >*carr,
    int size )
    {

    double buff = 0;

    for ( int i = 0; i < size; ++i )
        {

        carr[i] = std::complex < double > (
                      gsl_ran_exponential ( _R_G,
                                            1 ),
                      0 );
        buff += carr[i].real ();

        }

    buff = 1.0 / buff;

    for ( int i = 0; i < size; ++i )
        {

        carr[i] = carr[i] * buff;

        }

    }

void
sample_uniform_n_simplex2 (
    std::complex < double >*carr,
    int size )
    {

    // we generate an Uniform Ordered array with Sum arr[i] = 1
    // 1) generate n exp distributed x

    carr[0] = std::complex < double > (
                  0,
                  0 );

    for ( int i = 1; i < size; ++i )
        {

        carr[i] = std::complex < double > (
                      gsl_rng_uniform ( _R_G ),
                      0 );
        }

    sort_vec_asc ( carr, size );

    for ( int i = 0; i < size - 1; ++i )
        {

        carr[i] = ( carr[i + 1] - carr[i] );

        }

    carr[size - 1] = ( 1.0, 0.0 ) - carr[size - 1];

    }

void
matrix_show (
    std::complex < double >*cmat,
    int n )
    {

    for ( int i = 0; i < n; ++i )
        {

        for ( int j = 0; j < n; ++j )
            {

            std::cout << cmat[j + n * i];

            }

        std::cout << std::endl;

        }

    std::cout << std::endl << std::endl;

    }

void
matrix_show_real_part (
    std::complex < double >*cmat,
    int n )
    {

    std::cout.precision ( 2 );

    for ( int i = 0; i < n; ++i )
        {

        for ( int j = 0; j < n; ++j )
            {

            std::cout << cmat[j + n * i].real () << "  ";

            }

        std::cout << std::endl;

        }

    std::cout << std::endl << std::endl;

    }

void
matrix_show (
    const double *mat,
    int n )
    {

    for ( int i = 0; i < n; ++i )
        {

        for ( int j = 0; j < n; ++j )
            {

            std::cout << mat[j + n * i] << " ";

            }

        std::cout << std::endl;

        }

    std::cout << std::endl << std::endl;

    }

template < typename T >
void vector_show (
    T * vec,
    int n )
    {

    for ( int i = 0; i < n; ++i )
        {

        std::cout << vec[i] << " ";

        }

    std::cout << std::endl << std::endl;

    }
    
   
    

void
vector_write_to_file (
    double * vec,
    double res,
    const char * s,
    int n )
    {

    std::ofstream T;
    T.open ( s, std::fstream::app );

    for ( int i = 0; i < n; ++i )
        {

        T << vec[i] << "\t";

        }

    T << res << "\t";

    T << std::endl ;
    T.close();

    }



void
vector_show (
    std::complex < double >*cvec,
    int n )
    {

    for ( int i = 0; i < n; ++i )
        {

        std::cout << cvec[i] << " ";

        }

    std::cout << std::endl << std::endl;

    }

void
vector_dyadic_nxn (
    std::complex < double >*in1,
    std::complex < double >*in2,
    std::complex < double >*out,
    int size )
    {

    for ( int i = 0; i < size; ++i )
        {

        for ( int j = 0; j < size; ++j )
            {

            out[j + i * size] = in1[i] * std::conj ( in2[j] );

            }

        }

    }

void
dyadic_nxn (
    double *in1,
    double *in2,
    double *out,
    int size )
    {

    for ( int i = 0; i < size; ++i )
        {

        for ( int j = 0; j < size; ++j )
            {

            out[j + i * size] = in1[i] * in2[j];

            }

        }

    }

void
vector_dyadic_self (
    std::complex < double >*in1,
    std::complex < double >*out,
    int size )
    {

    for ( int i = 0; i < size; ++i )
        {

        for ( int j = 0; j < size; ++j )
            {

            out[j + i * size] = in1[i] * std::conj ( in1[j] );

            }

        }

    }
void
vector_dyadic_self2 (
    std::complex < double >*in1,
    std::complex < double >*out,
    int size_v,
    int size_m )
    {
    for ( int i = 0; i < size_v; ++i )
        {

        for ( int j = 0; j < size_v; ++j )
            {

            out[j + i * size_m] = in1[i] * std::conj ( in1[j] );

            }

        }

    }


void
vector_dyadic_self (
    double *in1,
    double *out,
    int dim_in1 )
    {

    for ( int i = 0; i < dim_in1; ++i )
        {

        for ( int j = 0; j < dim_in1; ++j )
            {

            out[j + i * dim_in1] = in1[i] * in1[j];

            }

        }

    }

void
matrix_transpose (
    std::complex < double >*in,
    std::complex < double >*out,
    int dim )
    {

    for ( int i = 0; i < dim; ++i )
        {

        for ( int j = 0; j < dim; ++j )
            {

            out[i + dim * j] = in[j + dim * i];
            }
        }
    }

void
matrix_transpose (
    std::complex < double >*in,
    int dim )
    {

    std::complex < double >temp;

    for ( int i = 0; i < dim; ++i )
        {

        for ( int j = i + 1; j < dim; ++j )
            {

            temp = in[i + dim * j];
            in[i + dim * j] = in[j + dim * i];
            in[j + dim * i] = temp;
            }
        }
    }

void
matrix_hermitian_conjugate (
    std::complex < double >*in,
    std::complex < double >*out,
    int dim )
    {

    for ( int i = 0; i < dim; ++i )
        {

        for ( int j = 0; j < dim; ++j )
            {

            out[i + dim * j] = std::conj ( in[j + dim * i] );

            }

        }

    }

void
matrix_hermitian_conjugate (
    std::complex < double >*in,
    int dim )
    {

    matrix_transpose ( in, dim );

    for ( int i = 0; i < dim; ++i )
        {

        for ( int j = 0; j < dim; ++j )
            {

            in[i + dim * j] = std::conj ( in[i + dim * j] );

            }

        }

    }

void
matrix_transpose (
    double *in,
    double *out,
    int dim )
    {

    for ( int i = 0; i < dim; ++i )
        {

        for ( int j = 0; j < dim; ++j )
            {

            out[i + dim * j] = in[j + dim * i];

            }

        }

    }

void
matrix_tensor_prod (
    std::complex < double >*in1,
    int dim_in1,
    std::complex < double >*in2,
    int dim_in2,
    std::complex < double >*out,
    int dim_out )
    {

    int index_out = dim_in1 * dim_in2;
    int m = 0;
    if ( index_out == dim_out )
        {

        for ( int i = 0; i < dim_in1; ++i )
            {

            for ( int k = 0; k < dim_in2; ++k )
                {

                for ( int j = 0; j < dim_in1; ++j )
                    {

                    for ( int l = 0; l < dim_in2; ++l )
                        {

                        out[m] = in1[i * dim_in1 + j] * in2[k * dim_in2 + l];
                        ++m;

                        }

                    }

                }

            }

        }

    else
        {

        std::cout <<
                  "ERROR in function: matrix_tensor_prod: INCOMPATIBLE DIMENSIONS";

        }

    }

void
matrix_tensor_prod (
    double *in1,
    int dim_in1,
    double *in2,
    int dim_in2,
    double *out,
    int dim_out )
    {

    int index_out = dim_in1 * dim_in2;
    int m = 0;
    if ( index_out == dim_out )
        {

        for ( int i = 0; i < dim_in1; ++i )
            {

            for ( int k = 0; k < dim_in2; ++k )
                {

                for ( int j = 0; j < dim_in1; ++j )
                    {

                    for ( int l = 0; l < dim_in2; ++l )
                        {

                        out[m] = in1[i * dim_in1 + j] * in2[k * dim_in2 + l];
                        m++;

                        }

                    }

                }

            }

        }

    else
        {

        std::cout << "ERROR IN TENSOR PRODUCT. INCOMPATIBLE DIMENSIONS";

        }

    }

void
transpose_block_mxm_at_offset_in_nxn (
    std::complex < double >*in,
    std::complex < double >*out,
    int dim,
    int offset_i,
    int offset_j,
    int dim_sub )
    {

    for ( int i = 0; i < dim_sub; ++i )
        {

        for ( int j = 0; j < dim_sub; ++j )
            {

            out[ ( offset_i + j ) * dim + offset_j +
                 i].real ( in[offset_j + j + dim * ( i + offset_i )].real () );
            out[ ( offset_i + j ) * dim + offset_j +
                 i].imag ( in[offset_j + j + dim * ( i + offset_i )].imag () );

            }

        }

    }

void
transpose_block_mxm_at_offset_in_nxn (
    double *in,
    double *out,
    int dim,
    int offset_i,
    int offset_j,
    int dim_sub )
    {

    for ( int i = 0; i < dim_sub; ++i )
        {

        for ( int j = 0; j < dim_sub; ++j )
            {

            out[ ( offset_i + j ) * dim + offset_j + i] =
                in[offset_j + j + dim * ( i + offset_i )];

            }

        }

    }

void
flip_block_kxk_in_block_mxm_in_nxn (
    double *in,
    double *out,
    int dim_in,
    int offset_i,
    int offset_j,
    int k,
    int m )
    {

    if ( m % k )
        {

        int stride = m / k;
        for ( int i = 0; i < stride; ++i )
            {

            for ( int j = 0; j < stride; ++j )
                {

                for ( int i_index = 0; i_index < k; ++i_index )
                    {

                    for ( int j_index = 0; j_index < k; ++j_index )
                        {

                        out[offset_j + i * k + j_index +
                            ( offset_i + j * k + i_index ) * dim_in] =
                                in[offset_j + j * k + j_index +
                                   ( offset_i + i * k + i_index ) * dim_in];

                        }

                    }

                }

            }

        }

    else
        {

        std::cout <<
                  " Non compatible dimensions: m dimension must be mutiple of k.";

        }

    }

void
flip_block_kxk_in_block_mxm_in_nxn (
    std::complex < double >*in,
    std::complex < double >*out,
    int dim_in,
    int offset_i,
    int offset_j,
    int k,
    int m )
    {

    if ( m % k == 0 )
        {

        int stride = m / k;
        for ( int i = 0; i < stride; ++i )
            {

            for ( int j = 0; j < stride; ++j )
                {

                for ( int i_index = 0; i_index < k; ++i_index )
                    {

                    for ( int j_index = 0; j_index < k; ++j_index )
                        {

                        //std::cout << " in funktion out=" <<  offset_j+i*k+j_index+(j*k+i_index)*dim_in+offset_i << std::endl;
                        //std::cout << " in funktion in =" <<  offset_j+j*k+j_index+(i*k+i_index)*dim_in+offset_i << std::endl;

                        out[offset_j + i * k + j_index +
                            ( j * k + i_index ) * dim_in + offset_i] =
                                in[offset_j + j * k + j_index +
                                   ( i * k + i_index ) * dim_in + offset_i];

                        }

                    }

                }

            }

        }

    else
        {

        std::cout <<
                  "Error in Subrutine flip_block_kxk_in_block_mxm_in_nxn: Non compatible dimensions: m dimension must be mutiple of k.";

        }

    }

void
state_partial_transpose_nxn (
    std::complex < double >*in,
    std::complex < double >*out,
    int dim_in,
    int *transpose_vector,
    int tot_nr_subspaces )
    {

    int sub_mat_size = dim_in;

    int block_size = dim_in;

    int steps = 0;

    for ( int i = 0; i < tot_nr_subspaces; ++i )
        {

        sub_mat_size = block_size;

        block_size /= abs ( transpose_vector[i] );

        steps = dim_in / sub_mat_size;

        //std::cout << "steps =" << (int) steps << std::endl;
        //std::cout << "sub_mat_size=" << sub_mat_size << std::endl;
        //std::cout << "block_size =" << block_size<< std::endl;

        if ( transpose_vector[i] < 0 )
            {

            for ( int j = 0; j < steps; ++j )
                {

                for ( int k = 0; k < steps; ++k )
                    {

                    // std::cout << "offset i =" << j*dim_in*sub_mat_size << std::endl ;

                    //  std::cout << "offset j =" << k*sub_mat_size << std::endl ;

                    flip_block_kxk_in_block_mxm_in_nxn ( in, out, dim_in,
                                                         j * dim_in *
                                                         sub_mat_size,
                                                         k * sub_mat_size,
                                                         block_size,
                                                         sub_mat_size );

                    }

                }

            }

        }

    }

std::complex < double >
matrix_trace_nxn (
    std::complex < double >*in,
    int dim_in )
    {

    std::complex < double >result;
    for ( int i = 0; i < dim_in; ++i )
        {

        result += in[ ( dim_in + 1 ) * i];

        }

    return result;

    }


void
state_decohere_on_subspace_ran (
    std::complex < double >*in,
    std::complex < double >*out,
    int dim_in,
    int *dec_vector,
    int tot_nr_subspaces )
    {

    int curr_dim = 1;
    int abs_dim = 0;
    int tot_dim = 1;
    int dim_dec = 0;		// <------------------------------------------------------
    int witness = 0;

    for ( int i = 0; i < tot_nr_subspaces; ++i )
        {
        tot_dim *= abs ( dec_vector[i] );
        if ( dec_vector[i] < 0 )
            {
            dim_dec = -dec_vector[i];
            ++witness;
            }
        }

    if ( witness > 1 )
        {
        std::cout <<
                  "Error in function Decorhere_on_subspace: only one subspace to decohere allowed";
        }

    std::complex < double >*total_op,
        *random_unitary, *temp1, *temp2, *basis_vector;

    total_op = new std::complex < double >[tot_dim * tot_dim];
    temp1 = new std::complex < double >[tot_dim * tot_dim];
    temp2 = new std::complex < double >[tot_dim * tot_dim];
    basis_vector = new std::complex < double >[tot_dim];
    random_unitary = new std::complex < double >[dim_dec * dim_dec];

    sample_unitary_matrix ( random_unitary, dim_dec );
    //matrix_initialize_unity( random_unitary ,dim_dec);
    if ( ( dim_in = !tot_dim ) )
        {

        std::cout <<
                  "Error in function decohere on subspace: dimension missmatch";

        return;
        }

    matrix_initialize_zero ( out, tot_dim );

    for ( int i = 0; i < dim_dec; ++i )
        {

        temp1[0] = std::complex < double > ( 1.0 , 0 );
        curr_dim = 1;

        for ( int j = tot_nr_subspaces - 1; j >= 0; --j )
            {

            abs_dim = abs ( dec_vector[j] );

            if ( dec_vector[j] < 0 )
                {

                matrix_get_column ( i, random_unitary, basis_vector, abs_dim );
                vector_dyadic_self ( basis_vector, temp2, abs_dim );
                //show_matrix(temp2,2);
                }

            else
                {

                matrix_initialize_unity ( temp2, abs_dim );
                //show_matrix(temp2,2);

                }

            matrix_tensor_prod ( temp2, abs_dim, temp1, curr_dim, total_op,
                                 curr_dim * abs_dim );
            //std::cout << curr_dim << std::endl;
            curr_dim *= abs_dim;
            //show_matrix(total_op,curr_dim);
            matrix_copy ( total_op, temp1, curr_dim );
            //std::cout << curr_dim << std::endl;
            //show_matrix(total_op,curr_dim);

            }

        matrix_mlt ( total_op, in, temp1, tot_dim );
        //show_matrix(temp1,tot_dim);
        matrix_mlt ( temp1, total_op, temp2, tot_dim );
        //show_matrix(temp2,tot_dim);
        matrix_add_to_first ( out, temp2, 1.0, tot_dim );
        //show_matrix(out,tot_dim);

        }

    delete[]temp1;
    delete[]temp2;
    delete[]total_op;
    delete[]random_unitary;
    delete[]basis_vector;

    }

void
state_decohere_on_subspace_ran (
    std::complex < double >*const in,
    std::complex < double >*out,
    int dim_in,
    int *dec_vector,
    int tot_nr_subspaces,
    std::complex < double >*temp1,
    std::complex < double >*temp2,
    std::complex < double >*total_op,
    std::complex < double >*random_unitary,
    std::complex < double >*basis_vector )
    {

    int curr_dim = 1;
    int abs_dim = 0;
    int tot_dim = 1;
    int dim_dec = 0;		// <------------------------------------------------------
    int witness = 0;
    matrix_initialize_zero ( out, tot_dim );

    for ( int i = 0; i < tot_nr_subspaces; ++i )
        {
        tot_dim *= abs ( dec_vector[i] );
        if ( dec_vector[i] < 0 )
            {
            dim_dec = -dec_vector[i];
            ++witness;
            }
        }
    if ( witness > 1 )
        {
        std::cout <<
                  "Error in function Decorhere_on_subspace: only one subspace to decohere allowed";
        }
    if ( ( dim_in = !tot_dim ) )
        {
        std::cout <<
                  "Error in function decohere on subspace: dimension missmatch";
        return;
        }

    sample_unitary_matrix ( random_unitary, basis_vector, temp1, dim_dec );

    //---------------------------------------------------------------------------------------------------------------------

    for ( int i = 0; i < dim_dec; ++i )
        {

        temp1[0] = std::complex < double > (
                       1.0,
                       0 );
        curr_dim = 1;

        for ( int j = tot_nr_subspaces - 1; j >= 0; --j )
            {

            abs_dim = abs ( dec_vector[j] );

            if ( dec_vector[j] < 0 )
                {

                matrix_get_column ( i, random_unitary, basis_vector, abs_dim );
                vector_dyadic_self ( basis_vector, temp2, abs_dim );
                }

            else
                {

                matrix_initialize_unity ( temp2, abs_dim );
                //show_matrix(temp2,2);

                }

            matrix_tensor_prod ( temp2, abs_dim, temp1, curr_dim, total_op,
                                 curr_dim * abs_dim );
            //std::cout << curr_dim << std::endl;
            curr_dim *= abs_dim;
            //show_matrix(total_op,curr_dim);
            matrix_copy ( total_op, temp1, curr_dim );
            //std::cout << curr_dim << std::endl;
            //matrix_show ( total_op, curr_dim );

            }

        matrix_mlt ( total_op, in, temp1, tot_dim );
        //show_matrix(temp1,tot_dim);
        matrix_mlt ( temp1, total_op, temp2, tot_dim );
        //show_matrix(temp2,tot_dim);
        matrix_add_to_first ( out, temp2, 1.0, tot_dim );
        //show_matrix(out,tot_dim);

        }

    }

void
state_decohere_on_subspace (
    const std::complex < double >*const in,
    std::complex < double >*out,
    const std::complex < double >* const unitary,
    int dim_in,
    int * dec_v,
    int tot_nr_subspaces )
    {

    int curr_dim = 1;
    int abs_dim = 0;
    int tot_dim = 1;
    int dim_dec = 0;		// <------------------------------------------------------
    int witness = 0;
   
    for ( int i = 0; i < tot_nr_subspaces; ++i )
        {
        tot_dim *= abs ( dec_v[i] );
        if ( dec_v[i] < 0 )
            {
            dim_dec = -dec_v[i];
            ++witness;
            }
        }

    if ( witness > 1 )
        {
        std::cout <<
                  "Error in function Decorhere_on_subspace: only one subspace to decohere allowed";
        }

    std::complex < double >*total_op,
        *temp1, *temp2, *basis_vector;

    total_op = new std::complex < double >[tot_dim * tot_dim];
    temp1 = new std::complex < double >[tot_dim * tot_dim];
    temp2 = new std::complex < double >[tot_dim * tot_dim];
    basis_vector = new std::complex < double >[tot_dim];

    if ( ( dim_in = !tot_dim ) )
        {

        std::cout <<
                  "Error in function decohere on subspace: dimension missmatch";

        return;
        }

    matrix_initialize_zero ( out, tot_dim );

    for ( int i = 0; i < dim_dec; ++i )
        {

        temp1[0] = std::complex < double > ( 1.0 , 0 );
        curr_dim = 1;

        for ( int j = tot_nr_subspaces - 1; j >= 0; --j )
            {

            abs_dim = abs ( dec_v[j] );

            if ( dec_v[j] < 0 )
                {

                matrix_get_column ( i, unitary, basis_vector, abs_dim );
                vector_dyadic_self ( basis_vector, temp2, abs_dim );
                //show_matrix(temp2,2);
                }

            else
                {

                matrix_initialize_unity ( temp2, abs_dim );
                //show_matrix(temp2,2);

                }

            matrix_tensor_prod ( temp2, abs_dim, temp1, curr_dim, total_op,
                                 curr_dim * abs_dim );
            //std::cout << curr_dim << std::endl;
            curr_dim *= abs_dim;
            //show_matrix(total_op,curr_dim);
            matrix_copy ( total_op, temp1, curr_dim );
            //std::cout << curr_dim << std::endl;
            //show_matrix(total_op,curr_dim);

            }

        matrix_mlt ( total_op, in, temp1, tot_dim );
        //show_matrix(temp1,tot_dim);
        matrix_mlt ( temp1, total_op, temp2, tot_dim );
        //show_matrix(temp2,tot_dim);
        matrix_add_to_first ( out, temp2, 1.0, tot_dim );
        //show_matrix(out,tot_dim);

        }

    delete[]temp1;
    delete[]temp2;
    delete[]total_op;
    delete[]basis_vector;

    }

void
state_decohere_on_subspace (
    const std::complex < double >*in,
    std::complex < double >*out,
    int dim_in,
    int *dec_vector,
    int tot_nr_subspaces,
    std::complex < double >*temp1,
    std::complex < double >*temp2,
    std::complex < double >*total_op,
    std::complex < double >*unitary,
    std::complex < double >*basis_vector )
    {

    int curr_dim = 1;
    int abs_dim = 0;
    int tot_dim = 1;
    int dim_dec = 0;		// <------------------------------------------------------
    int witness = 0;
    matrix_initialize_zero ( out, tot_dim );

    for ( int i = 0; i < tot_nr_subspaces; ++i )
        {
        tot_dim *= abs ( dec_vector[i] );
        if ( dec_vector[i] < 0 )
            {
            dim_dec = -dec_vector[i];
            ++witness;
            }
        }
    if ( witness > 1 )
        {
        std::cout <<
                  "Error in function Decorhere_on_subspace: only one subspace to decohere allowed";
        }
    if ( ( dim_in = !tot_dim ) )
        {
        std::cout <<
                  "Error in function decohere on subspace: dimension missmatch";
        return;
        }

    //---------------------------------------------------------------------------------------------------------------------

    for ( int i = 0; i < dim_dec; ++i )
        {

        temp1[0] = _1;
        curr_dim = 1;

        for ( int j = tot_nr_subspaces - 1; j >= 0; --j )
            {

            abs_dim = abs ( dec_vector[j] );

            if ( dec_vector[j] < 0 )
                {

                matrix_get_column ( i, unitary, basis_vector, abs_dim );
                vector_dyadic_self ( basis_vector, temp2, abs_dim );
                }

            else
                {

                matrix_initialize_unity ( temp2, abs_dim );
                //show_matrix(temp2,2);

                }

            matrix_tensor_prod ( temp2, abs_dim, temp1, curr_dim, total_op,
                                 curr_dim * abs_dim );
            //std::cout << curr_dim << std::endl;
            curr_dim *= abs_dim;
            //show_matrix(total_op,curr_dim);
            matrix_copy ( total_op, temp1, curr_dim );
            //std::cout << curr_dim << std::endl;
            //matrix_show ( total_op, curr_dim );

            }

        matrix_mlt ( total_op, in, temp1, tot_dim );
        //show_matrix(temp1,tot_dim);
        matrix_mlt ( temp1, total_op, temp2, tot_dim );
        //show_matrix(temp2,tot_dim);
        matrix_add_to_first ( out, temp2, 1.0, tot_dim );
        //show_matrix(out,tot_dim);

        }

    }


void
matrix_get_column (
    int col,
    const std::complex < double >* const cmat,
    std::complex < double >*out,
    int size )
    {

    for ( int i = 0; i < size; ++i )
        {

        out[i] = cmat[i * size + col];

        }

    }

void
matrix_get_partial_column (
    int col,
    int col_lenght,
    std::complex < double >*cmat,
    std::complex < double >*out,
    int size )
    {

    for ( int i = 0; i < col_lenght; ++i )
        {
        out[i] = cmat[i * size + col];
        }

    }
void
matrix_set_partial_column (
    int col,
    int col_lenght,
    std::complex < double >*cvec,
    std::complex < double >*out,
    int size )
    {

    for ( int i = 0; i < col_lenght; ++i )
        {
	  out[i * size + col]=cvec[i];
        }

    }
    

void
matrix_normalize (
    std::complex < double >*cmat,
    int size )
    {

    std::complex < double > norm = 1.0 / ( matrix_trace_nxn ( cmat, size ) );

    for ( int i = 0; i < size; ++i )
        {

        for ( int j = 0; j < size; ++j )
            {

            cmat[i * size + j] = cmat[i * size + j] * norm;

            }

        }

    }

void
matrix_add (
    std::complex < double >*cmat1,
    std::complex < double >*cmat2,
    int size )
    {

    for ( int i = 0; i < size; ++i )
        {

        for ( int j = 0; j < size; ++j )
            {

            cmat1[i * size + j] = cmat1[i * size + j] + cmat2[i * size + j];

            }

        }

    }

void
matrix_sub (
    std::complex < double >*cmat1,
    std::complex < double >*cmat2,
    int size )
    {

    for ( int i = 0; i < size; ++i )
        {

        for ( int j = 0; j < size; ++j )
            {

            cmat1[i * size + j] = cmat1[i * size + j] - cmat2[i * size + j];

            }

        }

    }

void
matrix_add_to_first (
    std::complex < double >*cmat1,
    std::complex < double >*cmat2,
    double factor,
    int size )
    {

    for ( int i = 0; i < size; ++i )
        {

        for ( int j = 0; j < size; ++j )
            {

            cmat1[i * size + j] += cmat2[i * size + j] * factor;

            }

        }

    }

void
vector_add (
    std::complex < double >*cvec1,
    double factor1,
    std::complex < double >*cvec2,
    double factor2,
    std::complex < double >*out,
    int size )
    {

    for ( int i = 0; i < size; ++i )
        {

        out[i] = cvec1[i] * factor1 + cvec2[i] * factor2;

        }

    }

void
vector_add_to_first (
    std::complex < double >*cvec1,
    std::complex < double >*cvec2,
    double factor,
    int size )
    {

    for ( int i = 0; i < size; ++i )
        {

        cvec1[i] += cvec2[i] * factor;

        }

    }

void
matrix_mlt (
    const std::complex < double >*cmat1,
    const std::complex < double >*cmat2,
    std::complex < double >*out,
    int size )
    {

    for ( int i = 0; i < size; ++i )
        {

        for ( int j = 0; j < size; ++j )
            {

            out[i * size + j] = std::complex < double > (
                                    0,
                                    0 );

            for ( int k = 0; k < size; ++k )
                {

                out[i * size + j] += cmat1[i * size + k] * cmat2[k * size + j];

                }

            }

        }

    }

void
matrix_scalar_mult (
    std::complex < double >*cmat,
    double factor,
    int size )
    {

    for ( int i = 0; i < size; ++i )
        {

        for ( int j = 0; j < size; ++j )
            {

            cmat[i * size + j] = cmat[i * size + j] * factor;

            }

        }

    }

void
matrix_scalar_mult (
    std::complex < double >*cmat,
    std::complex < double >factor,
    int size )
    {

    for ( int i = 0; i < size; ++i )
        {

        for ( int j = 0; j < size; ++j )
            {

            cmat[i * size + j] = cmat[i * size + j] * factor;

            }

        }

    }

void
matrix_copy (
    const std::complex < double > * const cmat,
    std::complex < double >*copy,
    int size )
    {

    for ( int i = 0; i < size; ++i )
        {

        for ( int j = 0; j < size; ++j )
            {

            copy[i * size + j] = cmat[i * size + j];

            }

        }

    }

template <typename T>
void
vector_copy (
    const T * vector,
    T * copy,
    int size )
    {

    for ( int i = 0; i < size; ++i )
        {
        copy[i] = vector[i];
        }

    }

    

void
vector_scalar_mult (
    std::complex < double >*cvec,
    double factor,
    int size )
    {

    for ( int i = 0; i < size; ++i )
        {

        cvec[i] = cvec[i] *  factor;

        }

    }

void
vector_scalar_mult (
    int*cvec,
    double factor,
    int size )
    {

    for ( int i = 0; i < size; ++i )
        {

        cvec[i] = cvec[i] *  factor;

        }

    }


void
vector_scalar_prod (
    std::complex < double >*cvec1,
    std::complex < double >*cvec2,
    std::complex < double >*out,
    int size )
    {

    *out = 0;
    for ( int i = 0; i < size; ++i )
        {

        *out += ( cvec1[i] * std::conj ( cvec2[i] ) );

        }

    }

void
matrix_vector_mult (
    std::complex < double >*cmat,
    std::complex < double >*cvec,
    std::complex < double >*out,
    int size )
    {

    for ( int i = 0; i < size; ++i )
        {

        out[i] = 0;
        for ( int j = 0; j < size; ++j )
            {

            out[i] += cmat[i * size + j] * cvec[j];

            }

        }

    }

void
vector_norm (
    std::complex < double >*cvec1,
    std::complex < double >*out,
    int size )
    {

    *out = 0;
    for ( int i = 0; i < size; ++i )
        {

        *out += ( cvec1[i] * std::conj ( cvec1[i] ) );

        }

    *out = sqrt ( *out );

    }

std::complex < double >
vector_norm (
    std::complex < double >*cvec1,
    int size )
    {

    std::complex < double >result;
    for ( int i = 0; i < size; ++i )
        {

        result += ( cvec1[i] * std::conj ( cvec1[i] ) );

        }

    return result;

    }

void
state_get_projector (
    std::complex < double >*unitary,
    int nr,
    std::complex < double >*out,
    int size )
    {

    for ( int i = 0; i < size; ++i )
        {

        for ( int j = 0; j < size; ++j )
            {

            out[j + i * size] =
                unitary[i * size + nr] * std::conj ( unitary[j * size + nr] );

            }

        }

    }

double
state_v_n_entropy (
    const std::complex < double > * const cmat,
    int size )
    {

    double result = 0;

    std::complex < double >*eigen;

    eigen = new std::complex < double >[size];

    matrix_eigenvalues ( cmat, eigen, size );

    vector_cutoff_neg ( eigen, size );

    for ( int i = 0; i < size; ++i )
        {

        if ( eigen[i].real () > 0 )
            {

            result += eigen[i].real () * log2 ( eigen[i].real () );

            }

        else if ( eigen[i].real () < 0 )
            {

            std::cout << " ERROR: DENSITY MATRIX HAS NEGATIVE EIGENVALUES ";
            return -1;

            }

        }
    delete[]eigen;
    return -result;

    }

double
state_v_n_entropy (
    const std::complex < double > * const cmat,
    std::complex < double >*copy,
    std::complex < double >*eigen, // eigenvalues vector
    int size )
    {

    double result = 0;

    matrix_eigenvalues ( cmat, copy, eigen, size );
    for ( int i = 0; i < size; ++i )
        {

        if ( eigen[i].real () > 0 )
            {

            result += eigen[i].real () * log2 ( eigen[i].real () );

            }

        else if ( eigen[i].real () < 0 )
            {

            std::cout << " ERROR: DENSITY MATRIX HAS NEGATIVE EIGENVALUES ";
            return -1;

            }

        }

    return -result;

    }
    
    

double
state_rel_entropy (
    const std::complex< double > * in1,
    const std::complex< double > * in2,
    std::complex< double > * in1_c,
    std::complex< double > * in2_c,
    std::complex< double > * temp_cmat_1,
    std::complex< double > * temp_cmat_2,
    std::complex< double > * temp_cmat_3,
    std::complex< double > * cvec,
    double  * eigen_1,
    double  * eigen_2,
    double  * eigen_temp,
    const int size )
    {

    matrix_log_2_hermitian ( in1,in1_c,temp_cmat_1,in2_c,eigen_temp,eigen_1,size );

    vector_cutoff_neg ( eigen_1,size );

    matrix_mlt ( in1,in2_c,temp_cmat_2,size ); //<-------------- temp_cmat_2 = Rho*log(Rho).

    matrix_copy ( in2,in2_c,size );

    LAPACKE_zheevd ( LAPACK_ROW_MAJOR,'V','U',size,in2_c,size,eigen_2 );

    matrix_hermitian_conjugate ( in2_c,in1_c,size ); //"Store EigenMat^H Sigma in inc_1"

    vector_cutoff_neg ( eigen_2,size );

    matrix_mlt ( in1,in2_c,temp_cmat_3,size ); // store Rho * EigenVecMat Sigma in temp_cmat_3

    for ( int i = 0; i < size; ++i )
        {
        if ( eigen_2[i] < 0.0 )
            {
            std::cout << "Relative entropy error: negative eigenvalue" << std::endl;
            return -1;
            }
        else if ( eigen_2[i] == 0 )
            {

            matrix_get_column ( i, temp_cmat_3, cvec, size );

            vector_cutoff_neg ( cvec,size );


            if ( vector_coef_sum ( cvec ,size ).real() != 0 )
                {

                return std::numeric_limits< double >::infinity();

                }

            }
        else if ( eigen_2[i] > 0.0 )

            {

            eigen_2[i]= log2 ( eigen_2[i] );

            }

        }

    for ( int i = 0; i < size; i++ )

        {

        for ( int j = 0; j < size; j++ )

            {

            temp_cmat_1[j * size + i] = temp_cmat_3[j * size + i] * eigen_2[i];

            }

        }

    matrix_mlt ( temp_cmat_1, in1_c, temp_cmat_3, size ); //temp_cmat_3 has now Rho*log(sigma)
    matrix_sub ( temp_cmat_2,temp_cmat_3,size );

    return matrix_trace_nxn ( temp_cmat_2 , size ).real();

    }


double
state_rel_entropy (
    const std::complex< double > * in1,
    const std::complex< double > * in2,
    const int size )
    {
    std::complex< double > * in1_c        = new std::complex< double >[size*size];
    std::complex< double > * in2_c        = new std::complex< double >[size*size];
    std::complex< double > * temp_cmat_1  = new std::complex< double >[size*size];
    std::complex< double > * temp_cmat_2  = new std::complex< double >[size*size];
    std::complex< double > * temp_cmat_3  = new std::complex< double >[size*size];
    std::complex< double > * cvec  = new std::complex< double >[size];
    double  * eigen_1                     = new   double [size];
    double  * eigen_2                     = new   double [size];
    double  * eigen_temp                  = new   double [size];

    double rel_ent = state_rel_entropy ( in1, in2, in1_c, in2_c, temp_cmat_1, temp_cmat_2, temp_cmat_3, cvec ,eigen_1, eigen_2, eigen_temp, size );

    delete[] in1_c;
    delete[] in2_c;
    delete[] temp_cmat_1;
    delete[] temp_cmat_2;
    delete[] temp_cmat_3;
    delete[] cvec;
    delete[] eigen_1;
    delete[] eigen_2;
    delete[] eigen_temp;

    return rel_ent;
    }


double
state_negativity (
    std::complex < double >*cmat,
    std::complex < double >*cmat_buff,
    std::complex < double >*buff,
    int * subspaces,
    int size )
    {

    state_partial_transpose_nxn ( cmat,cmat_buff,size, subspaces, 2 );

    double result = 0.0;
    matrix_eigenvalues ( cmat_buff, buff, size );

    for ( int i = 0; i < size; ++i )
        {

        if ( buff[i].real () < 0.0 )
            {

            result += -buff[i].real ();

            }

        }

    return result;

    }

double
state_negativity (
    std::complex < double >*cmat,
    int * subspaces,
    int size )
    {
    double result;

    std::complex < double >*cmat_buff;
    cmat_buff = new std::complex < double >[size*size];

    std::complex < double >*eigen;
    eigen = new std::complex < double >[size];

    result = state_negativity ( cmat, cmat_buff, eigen, subspaces, size );

    delete[]eigen;

    return result;

    }

void
vector_tensor_prod (
    std::complex < double >*in1,
    int dim_in1,
    std::complex < double >*in2,
    int dim_in2,
    std::complex < double >*out,
    int dim_out )
    {

    if ( dim_in1 * dim_in2 != dim_out )
        {

        std::cout << "Error in funtion VECTOR_TENSOR_PROD: dimension Mismatch";
        return;

        }

    for ( int i = 0; i < dim_in1; ++i )
        {

        for ( int j = 0; j < dim_in2; ++j )
            {

            out[i * dim_in2 + j] = in1[i] * in2[j];

            }

        }

    }

double
state_hs_norm (
    std::complex < double >*in,
    int size )
    {

    std::complex < double >norm = std::complex < double > (
                                      0, 0 );

    for ( int i = 0; i < size; ++i )
        {

        for ( int j = 0; j < size; ++j )
            {

            norm += in[j + i * size] * std::conj ( in[j + i * size] );

            }

        }

    return sqrt ( norm.real () );

    }

bool
compare_real_part (
    std::complex < double >a,
    std::complex < double >b )
    {

    return ( a.real () < b.real () );

    }

void
sort_vec_asc (
    std::complex < double >*eigen,
    int size )
    {

    std::sort ( eigen, eigen + size, compare_real_part );

    }

void
sort_vec_asc (
    double *eigen,
    int size )
    {

    std::sort ( eigen, eigen + size );

    }

double
state_purity (
    std::complex < double >*cmat,
    int size )
    {

    std::complex < double >purity = std::complex < double > (
                                        0, 0 );

    for ( int i = 0; i < size; ++i )
        {

        for ( int j = 0; j < size; ++j )
            {

            purity += cmat[j + i * size] * cmat[i + j * size];
            }
        }

    return purity.real ();
    }

std::complex < double >
vector_coef_sum (
    std::complex < double >*cvec,
    int size )
    {

    std::complex < double >coef_sum = std::complex < double > (
                                          0.0, 0.0 );

    for ( int i = 0; i < size; ++i )
        {

        coef_sum += cvec[i];
        }

    return coef_sum;

    }

void
sample_diag_state (
    std::complex < double >*cmat,
    int size )
    {

    std::complex < double >*diag_v;
    diag_v = new std::complex < double >[size];

    sample_uniform_n_simplex ( diag_v, size );

    for ( int i = 0; i < size * size; ++i )
        {
        cmat[i] = std::complex < double > (
                      0.0,
                      0.0 );
        }
    for ( int i = 0; i < size; ++i )
        {
        cmat[i * ( size + 1 )] = diag_v[i];
        }

    delete[]diag_v;
    }

void
matrix_initialize_unity (
    std::complex < double >*cmat,
    int size )
    {

    for ( int i = 0; i < size * size; ++i )
        {
        cmat[i] = std::complex < double > (
                      0.0,
                      0.0 );
        }
    for ( int i = 0; i < size; ++i )
        {
        cmat[i * ( size + 1 )] = std::complex < double > (
                                     1.0,
                                     0.0 );
        }
    }

void
matrix_initialize_zero (
    std::complex < double >*cmat,
    int size )
    {

    for ( int i = 0; i < size * size; ++i )
        {

        cmat[i] = std::complex < double > ( 0.0,0.0 );
        }

    }

void
matrix_cutoff_neg (
    std::complex < double >*in,
    int size )
    {

    for ( int i = 0; i < size*size; ++i )
        {


        if ( ( in[i].real () > -CUT_OFF_AS_ZERO ) && ( in[i].real () < 0.0 ) )
            {
            in[i].real ( 0.0 );
            }

        if ( ( in[i].imag () > -CUT_OFF_AS_ZERO ) && ( in[i].imag () < 0.0 ) )
            {
            in[i].imag ( 0.0 );
            }
        }
    }

void
vector_cutoff_neg (
    std::complex < double >*in,
    int size )
    {

    for ( int i = 0; i < size; ++i )
        {

        if ( ( in[i].real () > -CUT_OFF_AS_ZERO ) && ( in[i].real () < 0.0 ) )
            {
            in[i].real ( +0.0 );
            }

        if ( ( in[i].imag () > -CUT_OFF_AS_ZERO ) && ( in[i].imag () < 0.0 ) )
            {
            in[i].imag ( +0.0 );
            }

        }

    }

void
vector_cutoff_neg (
    double *in,
    int size )
    {

    for ( int i = 0; i < size; ++i )
        {

        if ( ( in[i] > -CUT_OFF_AS_ZERO ) && ( in[i] < 0.0 ) )
            {
            in[i] = 0;
            }

        }

    }


void
matrix_svd_full (
    std::complex < double >*in,
    std::complex < double >*U,
    double *sv,
    std::complex < double >*V_H,
    int size )
    {

    LAPACKE_zgesdd ( LAPACK_ROW_MAJOR, 'A', size, size, in, size, sv, U, size,
                     V_H, size );

    }

void
matrix_partial_cross_transpose (
    std::complex < double >*in,
    std::complex < double >*out,
    int dim_a,
    int dim_b,
    int size )
    {

    int index_out = dim_a * dim_b;

    if ( size == index_out )
        {

        for ( int i = 0; i < dim_a; ++i )
            {

            for ( int j = 0; j < dim_a; ++j )
                {

                for ( int k = 0; k < dim_b; ++k )
                    {

                    for ( int l = 0; l < dim_b; ++l )
                        {

                        out[ ( i * dim_b + j ) * size + k * dim_b + l] =
                            in[ ( i * dim_b + k ) * size + j * dim_b + l];
                        //std::cout << (i*dim_b+k)*size+j*dim_b+l << '\t' << (i*dim_b+j)*size+k*dim_b+l << std::endl;

                        }

                    }

                }

            }

        }

    else
        {

        std::cout <<
                  "ERROR in function: matrix_partial_cross_transpose: INCOMPATIBLE DIMENSIONS";

        }

    }

void
matrix_square_root (
    std::complex < double >*in,
    std::complex < double >*in_copy,
    std::complex < double >*temp,
    std::complex < double >*out,
    double *eigen_dvec,
    int size )
    {

    matrix_copy ( in, in_copy, size );

    LAPACKE_zheevd ( LAPACK_ROW_MAJOR, 'V', 'U', size, in_copy, size,
                     eigen_dvec );

    vector_cutoff_neg ( eigen_dvec, size );


    for ( int i = 0; i < size; i++ )
        {

        eigen_dvec[i] = sqrt ( eigen_dvec[i] );

        }
    // Rücktrafo : S*D*S^-1:
    // S*D

    for ( int i = 0; i < size; i++ )
        {
        for ( int j = 0; j < size; j++ )
            {
            temp[j * size + i] = in_copy[j * size + i] * eigen_dvec[i];
            }
        }

    matrix_hermitian_conjugate ( in_copy, size );

    matrix_mlt ( temp, in_copy, out, size );

    }


void
matrix_log_nat_hermitian ( const std::complex< double > * in ,
                           std::complex< double > * in_copy ,
                           std::complex< double > * temp,
                           std::complex< double > * out,
                           double * eigen_temp,
                           double * eigen_out,
                           int size )
    {


    matrix_copy ( in,in_copy,size );

    LAPACKE_zheevd ( LAPACK_ROW_MAJOR,'V','U',size,in_copy,size,eigen_temp );

    vector_cutoff_neg ( eigen_temp, size );
    vector_copy ( eigen_temp,eigen_out,size );

    for ( int i=0; i<size; i++ )
        {
        if ( eigen_temp[i] > 0 )
            {
            eigen_temp[i] = log ( eigen_temp[i] );
            }
        else if ( eigen_temp[i] < 0 )
            {
            std::cout <<  "Matrix log Error: Negative Eigenvalues" <<  std::endl;
            return;
            }
        else if ( eigen_temp == 0 )
            {
            std::cout <<  "Infinite Matrix logaritm" <<  std::endl;
            return;
            }
        }

    // Rücktrafo : S*D*S^-1:
    // S*D

    for ( int i=0; i<size; i++ )
        {
        for ( int j=0; j<size; j++ )
            {
            temp[j*size+i]=in_copy[j*size+i]*eigen_temp[i];
            }
        }

    matrix_hermitian_conjugate ( in_copy,size );

    matrix_mlt ( temp,in_copy,out,size );


    }

void
matrix_log_nat_hermitian ( std::complex< double > * in ,
                           std::complex< double > * out,
                           int size )
    {

    std::complex< double > * in_copy    = new std::complex< double >[size * size];
    std::complex< double > * temp       = new std::complex< double >[size * size];
    double                 * eigen_dvec = new double [size];
    double                 * eigen_out  = new double [size];

    matrix_log_nat_hermitian ( in , in_copy , temp , out, eigen_dvec ,eigen_out , size );

    delete[] in_copy;
    delete[] temp;
    delete[] eigen_dvec;
    delete[] eigen_out;

    }


void
matrix_log_2_hermitian ( const std::complex< double > * in ,
                         std::complex< double > * in_copy ,
                         std::complex< double > * temp,
                         std::complex< double > * out,
                         double * eigen_temp,
                         double * eigen_out,
                         int size )
    {


    matrix_copy ( in,in_copy,size );

    LAPACKE_zheevd ( LAPACK_ROW_MAJOR,'V','U',size,in_copy,size,eigen_temp );

    vector_cutoff_neg ( eigen_temp, size );
    vector_copy ( eigen_temp,eigen_out,size );

    for ( int i=0; i<size; i++ )
        {
        if ( eigen_temp[i] > 0 )
            {
            eigen_temp[i] = log2 ( eigen_temp[i] );
            }
        else if ( eigen_temp[i] < 0 )
            {
            std::cout <<  "Matrix log Error: Negative Eigenvalues" <<  std::endl;
            return;
            }
        else if ( eigen_temp == 0 )
            {
            std::cout <<  "Infinite Matrix logaritm" <<  std::endl;
            return;
            }
        }

    // Rücktrafo : S*D*S^-1:
    // S*D

    for ( int i=0; i<size; i++ )
        {
        for ( int j=0; j<size; j++ )
            {
            temp[j*size+i]=in_copy[j*size+i]*eigen_temp[i];
            }
        }

    matrix_hermitian_conjugate ( in_copy,size );

    matrix_mlt ( temp,in_copy,out,size );


    }

void
matrix_log_2_hermitian ( std::complex< double > * in ,
                         std::complex< double > * out,
                         int size )
    {

    std::complex< double > * in_copy    = new std::complex< double >[size * size];
    std::complex< double > * temp       = new std::complex< double >[size * size];
    double                 * eigen_dvec = new double [size];
    double                 * eigen_out = new double [size];

    matrix_log_2_hermitian ( in , in_copy , temp , out, eigen_dvec,eigen_out , size );

    delete[] in_copy;
    delete[] temp;
    delete[] eigen_dvec;
    delete[] eigen_out;

    }

void
matrix_square_root_hermitian (
    std::complex < double >*in,
    std::complex < double >*out,
    int size )
    {

    std::complex < double >*in_copy;
    std::complex < double >*temp;
    double *eigen_dvec;

    in_copy = new std::complex < double >[size * size];
    temp = new std::complex < double >[size * size];
    eigen_dvec = new double[size];

    matrix_square_root ( in, in_copy, temp, out, eigen_dvec, size );

    delete[]in_copy;
    delete[]temp;
    delete[]eigen_dvec;

    }

double
state_quantum_commutance (
    std::complex < double >*in,
    int dim_a,
    int dim_b,
    int dim_max,
    int dim_max_mat,
    int size,
    std::complex < double >*temp_cmat_1,
    std::complex < double >*temp_cmat_2,
    std::complex < double >*temp_cmat_3,
    double *temp_dvec,
    std::complex < double >*temp_cvec,
    std::complex < double >*X )
    {
    matrix_partial_cross_transpose ( in, temp_cmat_1, dim_a, dim_b, size );	// calculate R

    matrix_hermitian_conjugate ( temp_cmat_1, temp_cmat_2, dim_max_mat );	// calculate R^H

    matrix_mlt ( temp_cmat_1, temp_cmat_2, temp_cmat_3, dim_max_mat );	// calculate R*R^H

    matrix_square_root ( temp_cmat_3, temp_cmat_1, temp_cmat_2, X, temp_dvec, dim_max_mat );	// calculate X = sqrt(R*R^H)

    matrix_partial_cross_transpose ( X, temp_cmat_1, dim_max, dim_max, dim_max_mat );	// calculate X^<

    matrix_mlt ( temp_cmat_1, temp_cmat_1, temp_cmat_2, dim_max_mat );	//calculate X^< * X^<

    matrix_partial_cross_transpose ( temp_cmat_2, temp_cmat_3, dim_max, dim_max,
                                     dim_max_mat );

    matrix_mlt ( temp_cmat_1, X, temp_cmat_2, dim_max_mat );

    matrix_sub ( temp_cmat_3, temp_cmat_2, dim_max_mat );


    return matrix_trace_nxn ( temp_cmat_3, dim_max_mat ).real ();
    }

double
state_quantum_commutance (
    std::complex < double >*in,
    int dim_a,
    int dim_b,
    int size )
    {
    int dim_max;

    if ( dim_a > dim_b )
        {
        dim_max = dim_a;
        }
    else
        {
        dim_max = dim_b;
        }

    int dim_mat_max = dim_max * dim_max;
    int dim_mat_max_global = dim_max * dim_max * dim_max * dim_max;

    std::complex < double >*temp_cmat_1 =
        new std::complex < double >[dim_mat_max_global];
    std::complex < double >*temp_cmat_2 =
        new std::complex < double >[dim_mat_max_global];
    std::complex < double >*temp_cmat_3 =
        new std::complex < double >[dim_mat_max_global];
    double *temp_dvec =
        new double[dim_mat_max];
    std::complex < double >*temp_cvec =
        new std::complex < double >[dim_mat_max];
    std::complex < double >*X = new std::complex < double >[dim_mat_max_global];

    double result =
        state_quantum_commutance ( in, dim_a, dim_b, dim_max, dim_mat_max, size,
                                   temp_cmat_1, temp_cmat_2, temp_cmat_3,
                                   temp_dvec, temp_cvec, X );

    delete[]temp_cmat_1;
    delete[]temp_cmat_2;
    delete[]temp_cmat_3;
    delete[]temp_dvec;
    delete[]temp_cvec;
    delete[]X;

    return result;
    }


double
state_Mutual_information_I_A_B ( const std::complex< double > *  const rho_AB,
                                 int size_A,
                                 int size_B,
                                 int size_AB )
    {
    double S_A ( 0 ),S_B ( 0 ),S_AB ( 0 ),S_A_w_B ( 0 );

    std::complex< double > * rho_A = new std::complex< double > [size_A*size_A];
    std::complex< double > * rho_B = new std::complex< double > [size_B*size_B];

    int tr_out_a[2] = {-size_A,size_B};
    int tr_out_b[2] = {size_A,-size_B};


    matrix_partial_trace ( rho_AB,rho_A,tr_out_a,2,size_AB );
    matrix_partial_trace ( rho_AB,rho_B,tr_out_b,2,size_AB );

    S_A= state_v_n_entropy ( rho_A,size_A );
    S_B= state_v_n_entropy ( rho_B,size_B );
    S_AB= state_v_n_entropy ( rho_AB,size_AB );

    delete[] rho_A;
    delete[] rho_B;

    return S_A+S_B-S_AB;

    }


double
state_conditional_entropy ( const std::complex< double > * const rho_AB,
                            const std::complex< double > * const von_neumann_set,
                            int size_A,
                            int size_B,
                            int size_AB ) //decoheres AXB on subspace B
    {

    std::complex< double > * rho_dec = new std::complex< double >[size_AB*size_AB];

    int dec_vec[2]= {size_A,-size_B};

    state_decohere_on_subspace ( rho_AB,rho_dec,von_neumann_set,size_AB,dec_vec,2 );

    double result = state_v_n_entropy ( rho_dec,size_AB );

    delete[] rho_dec;

    return result;


    }



void matrix_partial_trace ( const std::complex< double > * in,
                            std::complex< double > * out,
                            int * trace_vector,
                            int nr_subspaces,
                            int size_tot )
    {

    int temp_size = size_tot;
    int block_size = size_tot;
    int sub_block_size = size_tot;
    int out_size=size_tot;
    int nr_blocks = 1;
    int ind = 0;

    std::complex< double > temp;
    std::complex< double > * buffer = new std::complex< double >[size_tot*size_tot];
    std::complex< double > * buffer2 = new std::complex< double >[size_tot*size_tot];

    matrix_copy ( in,buffer,size_tot );

    for ( int subspace = 0 ; subspace < nr_subspaces ; ++subspace )
        {
        sub_block_size /= trace_vector[subspace];

        if ( trace_vector[subspace] <0 )
            {
            sub_block_size=-sub_block_size;
            out_size /=-trace_vector[subspace];
            nr_blocks=temp_size/block_size;

            for ( int block_j = 0 ; block_j < nr_blocks ; ++block_j )
                {
                for ( int block_i = 0 ; block_i < nr_blocks ; ++block_i )
                    {
                    for ( int in_block_j = 0 ; in_block_j < sub_block_size ; ++in_block_j )
                        {
                        for ( int in_block_i = 0 ; in_block_i < sub_block_size ; ++in_block_i )
                            {
                            temp=std::complex< double > ( 0,0 );
                            ind = in_block_i+block_i*block_size+ ( in_block_j +  block_j*block_size ) *temp_size;
                            for ( int sum_i=0; sum_i< -trace_vector[subspace]; ++sum_i )
                                {
                                temp +=buffer[ind + sum_i*sub_block_size* ( temp_size+1 )];
                                }
                            buffer2[in_block_i+block_i*sub_block_size+ ( in_block_j + block_j*sub_block_size ) *out_size ] = temp;
                            }
                        }
                    }
                }
            matrix_copy ( buffer2,buffer,out_size );
            temp_size = out_size;
            block_size /= -trace_vector[subspace];
            }
        else
            {
            block_size /= trace_vector[subspace];
            }
        }
    matrix_copy ( buffer2,out,out_size );
    delete[] buffer;
    delete[] buffer2;
    }


    
    
    
    
    
    template void vector_copy <int> (const int * vector,
				     int *copy,
				     int size );
    
    template void vector_copy <std::complex< double >> (const std::complex< double > * vector,
							std::complex< double > *copy,
							int size );
    
    template void vector_show <int> (int * vec ,
				     int n);
    template void vector_show <double> (double * vec ,
				     int n);
    
    
    
    