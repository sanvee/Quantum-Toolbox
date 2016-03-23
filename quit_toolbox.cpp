#include "quit_toolbox.h"

std::complex< double > _0 = std::complex< double > ( 0.0,0.0 );
std::complex< double > _1 = std::complex< double > ( 1.0,0.0 );
std::complex< double > _neg_1 = std::complex< double > ( -1.0,0.0 );
std::complex< double > _i = std::complex< double > ( 0.0,1.0 );
std::complex< double > _neg_i = std::complex< double > ( 0.0,-1.0 );
double _1_[2] = {1,0};
double _0_[2] = {0,0};


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

    copy = new std::complex < double >[size * size]();

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

    for ( int i = 0; i < size*size; ++i )
        {
        cmat[i] = std::complex < double > ( gsl_ran_gaussian_ziggurat ( _R_G,1 ), gsl_ran_gaussian_ziggurat ( _R_G,1 ) );
        }

    }

void
sample_density_matrix_hs ( std::complex < double >*out,
                           int size )
    {
// generate a ginibre matrix X ; do U=X*(X^H) ; do U=U/tr(U); return U

    std::complex < double > * buffer;
    buffer = new std::complex < double >[size * size];

    sample_random_ginibre_matrix ( buffer, size );

    cblas_zgemm ( CblasRowMajor,CblasNoTrans,CblasConjTrans,size,size,size,_1_,buffer,size,buffer,size,_0_,out,size );

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

            out[j + size * i] = std::complex < double > ( 0,0 );

            for ( int k = 0; k < size; ++k )
                {
                out[j + size * i] += buff[k + i * size] * tau[k] * std::conj ( buff[k + j * size] );
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

        temp_v[i] = cmat[i * ( size + 1 )] / std::abs ( cmat[i * ( size + 1 )] );

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

        temp[i] = cmat[i * ( size + 1 )] / std::abs ( cmat[i * ( size + 1 )] );

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
sample_pure_density_matrix (
    std::complex < double >*out,
    int size )
    {

    std::complex < double >*cvec;
    cvec = new std::complex < double >[size];

    sample_unitary_matrix ( out, size );
    matrix_get_column ( 0, out, cvec, size );
    vector_dyadic_self ( cvec, out, size );

    delete[]cvec;

    }

void sample_pure_density_matrix2 ( std::complex < double >*out,
                                   int size )
    {
    std::complex <double> * temp_cmat1 = new std::complex<double>[size*size];

    sample_random_ginibre_matrix ( temp_cmat1, size );

    for ( int i = 0 ; i < size ; ++i )
        {
        for ( int j = 0; j < size ; ++j )
            {
            for ( int inner = 0; inner < size; ++inner )
                {
                out[i*size+j] = temp_cmat1[i*size+inner] * std::conj ( temp_cmat1[j*size+inner] );
                }
            }
        }

    double t = matrix_trace_nxn ( out,size ).real();
    t= 1/t;

    matrix_scalar_mult ( out,t,size );
    delete[]temp_cmat1;

    }
void sample_separable_state ( std::complex <double> * out,int nr_mixtures ,int * subspaces , int nr_subspaces,int dim_tot ) // generates mixture of productstates with nr_para terms.
    {

    std::complex <double> * coefficients = new std::complex <double>[nr_mixtures]();
    std::complex <double> * temp_cmat1 = new std::complex<double>[dim_tot*dim_tot]();
    std::complex <double> * temp_cmat2 = new std::complex<double>[dim_tot*dim_tot]();
    std::complex <double> * temp_cmat3 = new std::complex<double>[dim_tot*dim_tot]();
    int curr_dim;
    abs ( subspaces[nr_subspaces-1] );
    matrix_initialize_zero ( out,dim_tot );

    sample_uniform_n_simplex ( coefficients,nr_mixtures );

    for ( int i = 0 ; i < nr_mixtures; ++i )
        {

        curr_dim = abs ( subspaces[nr_subspaces-1] );
        sample_pure_density_matrix ( temp_cmat1,curr_dim );

        for ( int j = nr_subspaces-2; j>=0 ; --j )
            {

            sample_pure_density_matrix ( temp_cmat2,abs ( subspaces[j] ) );
            matrix_tensor_prod ( temp_cmat2,abs ( subspaces[j] ),temp_cmat1,curr_dim,temp_cmat3,abs ( curr_dim*subspaces[j] ) );
            curr_dim *=  abs ( subspaces[j] );
            matrix_copy ( temp_cmat3,temp_cmat1,curr_dim );
            }
        matrix_add_to_first ( out,temp_cmat1,coefficients[i],dim_tot );

        }

    delete[]coefficients;
    delete[]temp_cmat1;
    delete[]temp_cmat2;
    delete[]temp_cmat3;

    }

void
sample_uniform_n_simplex ( std::complex < double >*carr,
                           int size )
    {

    double buff = 0;

    for ( int i = 0; i < size; ++i )
        {

        carr[i] = std::complex < double > ( gsl_ran_exponential ( _R_G, 1 ),0 );
        buff += carr[i].real ();

        }

    buff = 1.0 / buff;

    for ( int i = 0; i < size; ++i )
        {

        carr[i] = carr[i] * buff;

        }
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
void
matrix_show (
    T * mat,
    int dim_i, int dim_j )
    {

    for ( int i = 0; i < dim_i; ++i )
        {

        for ( int j = 0; j < dim_j; ++j )
            {

            std::cout << mat[j + dim_j * i] << " ";

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

void matrix_hermitian_conjugate (
    std::complex < double >*in,
    std::complex < double >*out,
    int dim_i, int dim_j )
    {

    for ( int i = 0; i < dim_i; ++i )
        {

        for ( int j = 0; j < dim_j; ++j )
            {

            out[j*dim_i+i] = std::conj ( in[i*dim_j+j] );

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

            out[ ( offset_i + j ) * dim + offset_j +i].real ( in[offset_j + j + dim * ( i + offset_i )].real () );
            out[ ( offset_i + j ) * dim + offset_j +i].imag ( in[offset_j + j + dim * ( i + offset_i )].imag () );

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

    std::complex < double > * temp = new std::complex< double > [dim_in * dim_in];

    matrix_copy ( in, temp, dim_in );

    for ( int i = 0; i < tot_nr_subspaces; ++i )
        {

        sub_mat_size = block_size;

        block_size /= abs ( transpose_vector[i] );

        steps = dim_in / sub_mat_size;

        if ( transpose_vector[i] < 0 )
            {
            for ( int j = 0; j < steps; ++j )
                {

                for ( int k = 0; k < steps; ++k )
                    {

                    // std::cout << "offset i =" << j*dim_in*sub_mat_size << std::endl ;

                    //  std::cout << "offset j =" << k*sub_mat_size << std::endl ;

                    flip_block_kxk_in_block_mxm_in_nxn ( temp, out, dim_in,
                                                         j * dim_in *
                                                         sub_mat_size,
                                                         k * sub_mat_size,
                                                         block_size,
                                                         sub_mat_size );

                    }

                }
            matrix_copy ( out, temp,dim_in );
            }

        }
    delete[] temp;
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
        matrix_add_to_first ( out, temp2, std::complex<double> ( 1,0 ), tot_dim );
        //show_matrix(out,tot_dim);

        }

    delete[]temp1;
    delete[]temp2;
    delete[]total_op;
    delete[]random_unitary;
    delete[]basis_vector;

    }


void
state_decohere_on_subspace ( // this function calculates RHO_A/B
    const std::complex < double >*const in,
    std::complex < double >*out,
    const std::complex < double >* const unitary,
    int dim_in,
    int * dec_v, // vector that gives the subspace on witch to decohere: {2,2,-2} will apply 1x1xP ...
    int tot_nr_subspaces )
    {

    int curr_dim = 1;
    int abs_dim = 0;
    int tot_dim = 1;
    int dim_dec = 0;		// <-----------dimension of the space to decohere
    int witness = 0;

    for ( int i = 0; i < tot_nr_subspaces; ++i )   // <--- checking that only one subspace is decohered
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

    for ( int i = 0; i < dim_dec; ++i )               // <---- now constructing the projectors.
        {

        temp1[0] = std::complex < double > ( 1.0 , 0 );
        curr_dim = 1;

        for ( int j = tot_nr_subspaces - 1; j >= 0; --j )    // starting in reverse oder because of the tensor product.
            {

            abs_dim = abs ( dec_v[j] );

            if ( dec_v[j] < 0 )
                {

                matrix_get_column ( i, unitary, basis_vector, abs_dim );  // takes the i th column of the unitary matrix and construct a Projektor
                vector_dyadic_self ( basis_vector, temp2, abs_dim );
                //show_matrix(temp2,2);
                }

            else
                {

                matrix_initialize_unity ( temp2, abs_dim );  // in the case of the invariant subspace the unity matrix is taken
                //show_matrix(temp2,2);

                }

            matrix_tensor_prod ( temp2, abs_dim, temp1, curr_dim, total_op,curr_dim * abs_dim ); // construct tensorproduct and proceed the same with next subspace.
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
        matrix_add_to_first ( out, temp2,std::complex<double> ( 1.0,0.0 ), tot_dim ); //  here we build the sum    P x rho x P
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
        matrix_add_to_first ( out, temp2,std::complex<double> ( 1,0 ), tot_dim );
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
    double factor,
    std::complex < double >*out,
    int size )
    {
    for ( int i = 0; i < size*size ; ++i )
        {
        out[i] = cmat1[i] + factor * cmat2[i];
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
    std::complex < double > factor,
    int size )
    {

    for ( int i = 0; i < size*size; ++i )
        {
        cmat1[i] += cmat2[i] * factor;
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

void matrix_mlt2 ( const std::complex< double >* cmat1,
                   const std::complex< double >* cmat2,
                   std::complex< double >* out, int size )
    {

    for ( int i = 0; i < size; ++i )
        {

        for ( int j = 0; j < size; ++j )
            {

            out[i * size + j] = std::complex < double > ( 0,0 );

            for ( int k = 0; k < size; ++k )
                {

                out[i * size + j] += cmat1[i * size + k] * cmat2[k * size + j];

                }

            }

        }

    }

void matrix_mlt ( const std::complex< double >* cmat1,
                  const std::complex< double >* cmat2,
                  std::complex< double >* out, int size )
    {

    cblas_zgemm ( CblasRowMajor,CblasNoTrans,CblasNoTrans,size,size,size,_1_,cmat1,size,cmat2,size,_0_,out,size );

    }


template < typename T >
void matrix_multiplication ( T * in_1, int dim_i_1, int dim_j_1, T * in_2, int dim_i_2 , int dim_j_2 ,T * out )
    {
    for ( int j_2 = 0 ; j_2 < dim_j_2; ++j_2 )
        {
        for ( int i_1 = 0 ; i_1 < dim_i_1; ++i_1 )
            {
            out[ j_2+i_1*dim_j_2]=T();
            for ( int inner = 0 ; inner < dim_j_1; ++inner )
                {
                out[ j_2+i_1*dim_j_2]+= in_1[i_1*dim_j_1+inner]*in_2[j_2+inner*dim_j_2];
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

    for ( int i = 0; i < size*size; ++i )
        {
        cmat[i] = cmat[i] * factor;
        }

    }

void
matrix_scalar_mult (
    std::complex < double >*cmat,
    std::complex < double >*out,
    std::complex < double > factor,
    int size )
    {

    for ( int i = 0; i < size*size; ++i )
        {
        out[i] = cmat[i] * factor;
        }

    }

void
matrix_scalar_mult (
    std::complex < double >*cmat,
    std::complex < double >factor,
    int size )
    {

    for ( int i = 0; i < size*size; ++i )
        {
        cmat[i] = cmat[i] * factor;
        }
    }

void
matrix_copy (
    const std::complex < double > * const cmat,
    std::complex < double >*copy,
    int size )
    {

    for ( int i = 0; i < size*size; ++i )
        {
        copy[i] = cmat[i];
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

            std::cout << " ERROR: DENSITY MATRIX HAS NEGATIVE EIGENVALUES: " << eigen[i] << std::endl;
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
    int nr_of_subspaces,// 2 by 3 would be {2,3}
    int size )
    {
    state_partial_transpose_nxn ( cmat,cmat_buff,size, subspaces, nr_of_subspaces );

    //matrix_show(cmat_buff,4);

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
    int nr_of_subspaces,
    int size )
    {
    double result;

    std::complex < double >*cmat_buff;
    cmat_buff = new std::complex < double >[size*size];

    std::complex < double >*eigen;
    eigen = new std::complex < double >[size];

    result = state_negativity ( cmat, cmat_buff, eigen, subspaces, nr_of_subspaces, size );

    delete[]eigen;
    delete[] cmat_buff;

    return result;

    }

double
state_logarithmic_negativity (
    std::complex < double >*cmat,
    int * subspaces,
    int nr_of_subspaces,
    int size )
    {
    double result;

    std::complex < double >*cmat_buff1 = new std::complex < double >[size*size];
    std::complex < double >*cmat_buff2 = new std::complex < double >[size*size];

    state_partial_transpose_nxn ( cmat,cmat_buff1,size, subspaces, nr_of_subspaces );

    cblas_zgemm ( CblasRowMajor,CblasConjTrans,CblasNoTrans,size,size,size,_1_,cmat_buff1,size,cmat_buff1,size,_0_,cmat_buff2,size );

    matrix_square_root ( cmat_buff2,cmat_buff1,size );

    result = log2 ( matrix_trace_nxn ( cmat_buff1,size ).real() );

    delete[] cmat_buff1;
    delete[] cmat_buff2;
    return result;


    }


double
state_hs_norm_herm (
    std::complex < double >*in,
    int size )
    {

    std::complex < double > result ( 0 );

    for ( int i = 0; i < size*size; ++i )
        {
        result +=   in[i]*std::conj ( in[i] );
        }
    return sqrt ( result.real() );
    }

double state_hilbert_schmidt_distance (
    std::complex < double > * in1,
    std::complex < double > * in2,
    int size )
    {

    std::complex< double > * temp1 = new std::complex< double >[size*size];
    std::complex< double > * temp2 = new std::complex< double >[size*size];

    matrix_add ( in1,in2,-1.0,temp1,size );
    matrix_mlt ( temp1,temp1,temp2,size );

    double result = sqrt ( matrix_trace_nxn ( temp2,size ).real() );
    delete[] temp1;
    delete[] temp2;
    return result;

    };



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

    std::complex < double > * temp = new std::complex < double >[size*size];

    std::complex < double >purity = std::complex < double > ( 0, 0 );

    matrix_mlt ( cmat,cmat,temp,size );

    purity = matrix_trace_nxn ( temp,size );

    delete[] temp;
    return purity.real();

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

                        out[ ( i * dim_a + j ) * dim_b*dim_b + k * dim_b + l] =
                            in[ ( i * dim_b + k ) * size        + j * dim_b + l];

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
    //vector_show(eigen_dvec,size);

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

    cblas_zgemm ( CblasRowMajor,CblasNoTrans,CblasConjTrans,size,size,size,_1_,temp,size,in_copy,size,_0_,out,size );

    // matrix_hermitian_conjugate ( in_copy, size );

    // matrix_mlt ( temp, in_copy, out, size );

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
            std::cout <<  "Error Matrix_log_nat_herm: Negative Eigenvalues" <<  std::endl;
            vector_show ( eigen_temp,size );
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
            std::cout <<  "Matrix matrix_log_2_hermitian: Negative Eigenvalues" <<  std::endl;
            vector_show ( eigen_out,size );
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
matrix_square_root (
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
    int size,
    std::complex < double >*temp_cmat_1,
    std::complex < double >*temp_cmat_2,
    std::complex < double >*temp_cmat_3,
    double *temp_dvec,
    std::complex < double >*X )
    {
    int dim_a2 = dim_a*dim_a;

    matrix_partial_cross_transpose ( in, temp_cmat_1, dim_a, dim_b, size );	// calculate R has dimension

    matrix_hermitian_conjugate ( temp_cmat_1, temp_cmat_2, dim_a2, dim_b*dim_b );	// calculate R^H <-------ADAPT TO RECTANGULAR

    matrix_multiplication ( temp_cmat_1,dim_a2,dim_b*dim_b, temp_cmat_2,dim_b*dim_b, dim_a2, temp_cmat_3 );	// calculate R*R^H  <-------ADAPT TO RECTANGULAR

    matrix_square_root ( temp_cmat_3, temp_cmat_1, temp_cmat_2, X, temp_dvec,dim_a2 );	// calculate X = sqrt(R*R^H)

    matrix_partial_cross_transpose ( X, temp_cmat_1, dim_a, dim_a, dim_a2 );	// calculate X^<

    matrix_mlt ( temp_cmat_1, temp_cmat_1, temp_cmat_2,dim_a2 );	//calculate X^< * X^<

    matrix_partial_cross_transpose ( temp_cmat_2, temp_cmat_3, dim_a, dim_a,
                                     dim_a2 );

    matrix_mlt ( temp_cmat_1, X, temp_cmat_2, dim_a2 );

    matrix_sub ( temp_cmat_3, temp_cmat_2, dim_a2 );


    return matrix_trace_nxn ( temp_cmat_3,dim_a2 ).real ();
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

    std::complex < double >*temp_cmat_1 = new std::complex < double >[dim_mat_max_global];
    std::complex < double >*temp_cmat_2 = new std::complex < double >[dim_mat_max_global];
    std::complex < double >*temp_cmat_3 = new std::complex < double >[dim_mat_max_global];

    double *temp_dvec = new double[dim_mat_max];
    std::complex < double >*X = new std::complex < double >[dim_mat_max_global];

    double result =
        state_quantum_commutance ( in, dim_a, dim_b,size,
                                   temp_cmat_1, temp_cmat_2, temp_cmat_3,
                                   temp_dvec,X );

    delete[]temp_cmat_1;
    delete[]temp_cmat_2;
    delete[]temp_cmat_3;
    delete[]temp_dvec;
    delete[]X;

    return result;
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
    std::complex< double > * buffer = new std::complex< double >[size_tot*size_tot]();
    std::complex< double > * buffer2 = new std::complex< double >[size_tot*size_tot]();

    matrix_copy ( in,buffer,size_tot );

    for ( int subspace = 0 ; subspace < nr_subspaces ; ++subspace )
        {
        sub_block_size /= trace_vector[subspace];

        if ( trace_vector[subspace] < 0 )
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
                            temp = std::complex< double > ( 0.0,0.0 );
                            ind  = in_block_i+block_i*block_size+ ( in_block_j +  block_j*block_size ) *temp_size;
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

double state_relative_entropy_of_coherence ( std::complex< double > * in, int size )
    {
    std::complex< double > * in_diag = new std::complex< double >[size*size];
    matrix_initialize_zero ( in_diag,size );

    for ( int i = 0; i < size ; ++i )
        {
        in_diag[i* ( size+1 )]=in[i* ( size+1 )];
        }

    double result = state_v_n_entropy ( in_diag,size )-state_v_n_entropy ( in,size );
    delete[] in_diag;
    return result;
    }

void state_trivial_expand ( std::complex< double > * in, std::complex< double > * out,int klein, int gross )

    {
    matrix_initialize_zero ( out,gross );

    for ( int i = 0; i < klein; ++i )
        {
        for ( int j = 0; j < klein; ++j )
            {

            out[i+j*gross]=in[i+j*klein];

            }

        }
    }
double trace_distance ( const std::complex< double > * in1 , std::complex< double > * in2, int size )
    {
    std::complex< double > * temp  = new std::complex< double > [size*size];
    std::complex< double > * eigen = new std::complex< double > [size];
    for ( int i = 0 ; i < size*size; ++i )
        {
        temp[i] = in1[i]-in2[i];
        }
    matrix_eigenvalues ( temp,eigen,size );

    double result ( 0 );

    for ( int i = 0 ; i < size; ++i )
        {
        result += std::abs ( eigen[i].real() );
        }
    delete[] temp;
    delete[] eigen;

    return 0.5*result;


    }

////////-----------------------------------------------------schrott !!! muss noch verallgemeinert werden ...
void
state_bit_flip (
    const std::complex < double >* in,
    std::complex < double >*out,
    double p,
    int size )
    {

    std::complex< double > * temp1 = new std::complex< double >[8*8];
    std::complex< double > * temp2 = new std::complex< double >[8*8];

    std::complex< double > E_0 [2*2] = { _1*sqrt ( 1-p ), _0,
                                         _0, _1*sqrt ( 1-p )
                                       };

    std::complex< double > E_1 [2*2] = {_0,_1*sqrt ( p ),
                                        _1*sqrt ( p ),_0
                                       };

    std::complex< double >  R [4*4];

    matrix_initialize_unity ( R,4 );
//matrix_show(E_0,2);
//matrix_show(E_1,2);

    matrix_tensor_prod ( R,4,E_0,2,temp1,8 );
    //matrix_show(temp1,8);

    cblas_zgemm ( CblasRowMajor,CblasNoTrans,CblasNoTrans,size,size,size,_1_,temp1,size,in,size,_0_,temp2,size );

    cblas_zgemm ( CblasRowMajor,CblasNoTrans,CblasConjTrans,size,size,size,_1_,temp2,size,temp1,size,_0_,out,size );

    matrix_tensor_prod ( R,4,E_1,2,temp1,8 );
    //matrix_show(temp1,8);

    cblas_zgemm ( CblasRowMajor,CblasNoTrans,CblasNoTrans,size,size,size,_1_,temp1,size,in,size,_0_,temp2,size );

    cblas_zgemm ( CblasRowMajor,CblasNoTrans,CblasConjTrans,size,size,size,_1_,temp2,size,temp1,size,_1_,out,size );

    //matrix_show(out,8);

    }

void
state_amp_damp (
    const std::complex < double >* in,
    std::complex < double >*out,
    double p,
    int size )
    {

    std::complex< double > * temp1 = new std::complex< double >[8*8];
    std::complex< double > * temp2 = new std::complex< double >[8*8];

    std::complex< double > E_0 [2*2] = { _0, _1* sqrt ( p ),
                                         _0,        _0
                                       };

    std::complex< double > E_1 [2*2] = {_1,_0,
                                        _0,_1*sqrt ( 1-p )
                                       };

    std::complex< double >  R [4*4];

    matrix_initialize_unity ( R,4 );
//matrix_show(E_0,2);
//matrix_show(E_1,2);

    matrix_tensor_prod ( R,4,E_0,2,temp1,8 );
    //matrix_show(temp1,8);

    cblas_zgemm ( CblasRowMajor,CblasNoTrans,CblasNoTrans,size,size,size,_1_,temp1,size,in,size,_0_,temp2,size );

    cblas_zgemm ( CblasRowMajor,CblasNoTrans,CblasConjTrans,size,size,size,_1_,temp2,size,temp1,size,_0_,out,size );

    matrix_tensor_prod ( R,4,E_1,2,temp1,8 );
    //matrix_show(temp1,8);

    cblas_zgemm ( CblasRowMajor,CblasNoTrans,CblasNoTrans,size,size,size,_1_,temp1,size,in,size,_0_,temp2,size );

    cblas_zgemm ( CblasRowMajor,CblasNoTrans,CblasConjTrans,size,size,size,_1_,temp2,size,temp1,size,_1_,out,size );

    //matrix_show(out,8);

    }
void
state_amp_damp_bob (
    const std::complex < double >* in,
    std::complex < double >*out,
    double p,
    int size )
    {
    std::complex< double > * temp = new std::complex< double >[8*8];
    std::complex< double > * temp1 = new std::complex< double >[8*8];
    std::complex< double > * temp2 = new std::complex< double >[8*8];

    std::complex< double > E_0 [2*2] = { _0, _1* sqrt ( p ),
                                         _0,        _0
                                       };

    std::complex< double > E_1 [2*2] = {_1,_0,
                                        _0,_1*sqrt ( 1-p )
                                       };

    std::complex< double >  R [2*2];

    matrix_initialize_unity ( R,2 );
//matrix_show(E_0,2);
//matrix_show(E_1,2);

    matrix_tensor_prod ( R,2,E_0,2,temp,4 );
    matrix_tensor_prod ( temp,4,R,2,temp1,8 );
    //matrix_show(temp1,8);

    cblas_zgemm ( CblasRowMajor,CblasNoTrans,CblasNoTrans,size,size,size,_1_,temp1,size,in,size,_0_,temp2,size );

    cblas_zgemm ( CblasRowMajor,CblasNoTrans,CblasConjTrans,size,size,size,_1_,temp2,size,temp1,size,_0_,out,size );

    matrix_tensor_prod ( R,2,E_1,2,temp,4 );
    matrix_tensor_prod ( temp,4,R,2,temp1,8 );

    cblas_zgemm ( CblasRowMajor,CblasNoTrans,CblasNoTrans,size,size,size,_1_,temp1,size,in,size,_0_,temp2,size );

    cblas_zgemm ( CblasRowMajor,CblasNoTrans,CblasConjTrans,size,size,size,_1_,temp2,size,temp1,size,_1_,out,size );

    //matrix_show(out,8);
    delete[] temp;
    delete[]temp1;
    delete[] temp2;

    }

void
state_phase_damping (
    const std::complex < double >* in,
    std::complex < double >*out,
    double p,
    int size )
    {

    std::complex< double > * temp1 = new std::complex< double >[8*8];
    std::complex< double > * temp2 = new std::complex< double >[8*8];

    std::complex< double > E_0 [2*2] = { _1*sqrt ( p/2.0 ), _0,
                                         _0, -_1*sqrt ( p/2.0 )
                                       };

    std::complex< double > E_1 [2*2] = {_1*sqrt ( 1- ( p/2.0 ) ),_0,
                                        _0,_1*sqrt ( 1- ( p/2.0 ) )
                                       };

    std::complex< double >  R [4*4];

    matrix_initialize_unity ( R,4 );
//matrix_show(E_0,2);
//matrix_show(E_1,2);

    matrix_tensor_prod ( R,4,E_0,2,temp1,8 );
    //matrix_show(temp1,8);

    cblas_zgemm ( CblasRowMajor,CblasNoTrans,CblasNoTrans,size,size,size,_1_,temp1,size,in,size,_0_,temp2,size );

    cblas_zgemm ( CblasRowMajor,CblasNoTrans,CblasConjTrans,size,size,size,_1_,temp2,size,temp1,size,_0_,out,size );

    matrix_tensor_prod ( R,4,E_1,2,temp1,8 );
    //matrix_show(temp1,8);

    cblas_zgemm ( CblasRowMajor,CblasNoTrans,CblasNoTrans,size,size,size,_1_,temp1,size,in,size,_0_,temp2,size );

    cblas_zgemm ( CblasRowMajor,CblasNoTrans,CblasConjTrans,size,size,size,_1_,temp2,size,temp1,size,_1_,out,size );

    //matrix_show(out,8);
    delete[] temp1;
    delete[] temp2;

    }

void
state_depolarizing (
    const std::complex < double >* in,
    std::complex < double >*out,
    double p,
    int size )
    {

    std::complex< double > * temp1 = new std::complex< double >[8*8];
    std::complex< double > * temp2 = new std::complex< double >[8*8];

    std::complex< double > E_0 [2*2] = { _1*sqrt ( 1- ( 3.0*p/4.0 ) ),                    _0,
                                         _0, _1*sqrt ( 1- ( 3.0*p/4.0 ) )
                                       };

    std::complex< double > E_1 [2*2] = {_0,                   _1*sqrt ( p/4.0 ),
                                        _1*sqrt ( p/4.0 )       ,_0
                                       };

    std::complex< double > E_2 [2*2] = { _0,                 -_i*sqrt ( p/4.0 ),
                                         _i*sqrt ( p/4.0 )   ,_0
                                       };

    std::complex< double > E_3 [2*2] = {_1*sqrt ( p/4.0 ),_0,
                                        _0,           -_1*sqrt ( p/4.0 )
                                       };


    std::complex< double >  R [4*4];

    matrix_initialize_unity ( R,4 );
//matrix_show(E_0,2);
//matrix_show(E_1,2);

    matrix_tensor_prod ( R,4,E_0,2,temp1,8 );
    //matrix_show(temp1,8);

    cblas_zgemm ( CblasRowMajor,CblasNoTrans,CblasNoTrans,size,size,size,_1_,temp1,size,in,size,_0_,temp2,size );

    cblas_zgemm ( CblasRowMajor,CblasNoTrans,CblasConjTrans,size,size,size,_1_,temp2,size,temp1,size,_0_,out,size );

    matrix_tensor_prod ( R,4,E_1,2,temp1,8 );
    //matrix_show(temp1,8);

    cblas_zgemm ( CblasRowMajor,CblasNoTrans,CblasNoTrans,size,size,size,_1_,temp1,size,in,size,_0_,temp2,size );

    cblas_zgemm ( CblasRowMajor,CblasNoTrans,CblasConjTrans,size,size,size,_1_,temp2,size,temp1,size,_1_,out,size );

    matrix_tensor_prod ( R,4,E_2,2,temp1,8 );
    //matrix_show(temp1,8);

    cblas_zgemm ( CblasRowMajor,CblasNoTrans,CblasNoTrans,size,size,size,_1_,temp1,size,in,size,_0_,temp2,size );

    cblas_zgemm ( CblasRowMajor,CblasNoTrans,CblasConjTrans,size,size,size,_1_,temp2,size,temp1,size,_1_,out,size );

    matrix_tensor_prod ( R,4,E_3,2,temp1,8 );
    //matrix_show(temp1,8);

    cblas_zgemm ( CblasRowMajor,CblasNoTrans,CblasNoTrans,size,size,size,_1_,temp1,size,in,size,_0_,temp2,size );

    cblas_zgemm ( CblasRowMajor,CblasNoTrans,CblasConjTrans,size,size,size,_1_,temp2,size,temp1,size,_1_,out,size );

    //matrix_show(out,8);
    delete[] temp1;
    delete[] temp2;

    }
    
    //------------------------template declaration------------------

template void vector_copy <int> ( const int * vector,
                                  int *copy,
                                  int size );

template void vector_copy <std::complex< double >> ( const std::complex< double > * vector,
        std::complex< double > *copy,
        int size );

template void vector_show <int> ( int * vec ,
                                  int n );
template void vector_show <double> ( double * vec ,
                                     int n );

template void matrix_multiplication <std::complex<double>> ( std::complex <double> * in_1,
        int dim_i_1,
        int dim_j_1,
        std::complex <double> * in_2,
        int dim_i_2,
        int dim_j_2,
        std::complex <double> * out );

template void matrix_show <std::complex< double >> ( std::complex< double > * mat, int dim_i, int dim_j );





