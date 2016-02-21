#include "quit_jarlskog.h"


//CLASS functions

void jarlskog_matrix::calculate_repesentation ( void )
    {

    representation[0]=eigenvalues[0];

    for ( int n = 1; n < size; ++n )
        {
        j_expand ( n,n+1 );
        create_V_n_n ( n ); // n = current size - 1
        matrix_hermitian_conjugate ( temp_V,temp_V_dag,n+1 );
        representation [n * ( n+2 )] = eigenvalues[n];
        matrix_mlt ( temp_V_dag,representation,temp_representation,n+1 );
        matrix_mlt ( temp_representation,temp_V,representation,n+1 );
        }

    }

void jarlskog_matrix::create_V_n_n ( int n )
    {

    matrix_initialize_zero ( temp_V,n+1 );
    matrix_get_partial_column ( n-1,n,vectors_normalized,temp_cvec,size-1 );
    vector_dyadic_self2 ( temp_cvec,temp_V,n,n+1 );

    std::complex< double > eins = std::complex< double > ( 1,0 );
    std::complex< double > m_c_eins = std::complex< double > ( -1,-1 );
    c = cos ( norms[n-1] );
    s = sin ( norms[n-1] );

    matrix_scalar_mult ( temp_V, - ( 1-c ),n+1 );

    for ( int i = 0; i < n; ++i )
        {
        temp_V[i* ( n+2 )]= eins+temp_V[i* ( n+2 )];
        }

    vector_scalar_mult ( temp_cvec,s,n );

    for ( int i = 0; i < n; ++i )
        {
        temp_V[n+i* ( n+1 )] = temp_cvec[i];
        temp_V[i+n* ( n+1 )] = -std::conj ( temp_cvec[i] );
        }

    temp_V[n * ( n+2 )]= c;

    }


void jarlskog_matrix::set_eigen ( std::complex< double >* eigen_cvec )
    {

    for ( int i = 0 ; i < size ; ++i )
        {
        eigenvalues[i] = eigen_cvec[i];
        }

    }
void jarlskog_matrix::set_eigen_random ()
    {
    sample_uniform_n_simplex ( eigenvalues,size );
    }


void jarlskog_matrix::set_vectors_random ( void )
    {

    for ( int n = 0; n < size-1; ++n )
        {
        double test ( 0.0 );
        sample_from_complex_sphere_n_rand ( temp_cvec,n+1 );
        for ( int i = 0 ; i < n+1 ; ++i )
            {
            vectors_normalized[n+i* ( size-1 )]=temp_cvec[i];
            }
        }
    }

void jarlskog_matrix::set_norms ( double  * norms )
    {
    for ( int i = 0; i < size-1; ++i )
        {
        this->norms[i] =  norms[i];
        }

    }

void jarlskog_matrix::set_norms_random()
    {
    for ( int i = 0 ; i < size-1; ++i )
        {
        norms[i] = gsl_rng_uniform ( _R_G );
        }

    }
void jarlskog_matrix::set_random()
    {
    this->set_eigen_random();
    this->set_norms_random();
    this->set_vectors_random();
    }


void jarlskog_matrix::j_expand ( int klein, int gross )
    {
    matrix_initialize_zero ( temp_representation,gross );
    for ( int i = 0; i < klein; ++i )
        {
        for ( int j = 0; j < klein; ++j )
            {

            temp_representation[i+j*gross]=representation[i+j*klein];
            }

        }

    for ( int i = 0; i < gross*gross; ++i )
        {
        representation[i]=temp_representation[i];
        }

    }

void jarlskog_matrix::minimize_rel_entropy ( std::complex< double >* ref, double tolerance )
    {

    std::complex< double > * in1_c        = new std::complex< double >[size*size];
    std::complex< double > * in2_c        = new std::complex< double >[size*size];
    std::complex< double > * temp_cmat_1  = new std::complex< double >[size*size];
    std::complex< double > * temp_cmat_2  = new std::complex< double >[size*size];
    std::complex< double > * temp_cmat_3  = new std::complex< double >[size*size];
    std::complex< double > * cvec         = new std::complex< double >[size];
    double  * eigen_1                     = new   double [size];
    double  * eigen_2                     = new   double [size];
    double  * eigen_temp                  = new   double [size];
    double * grad = new double[size]();


    double x_0 ( 0 );
    double f_x_0 ( 0 );
    double y_0 ( 1 );
    double f_y_0 = ( 0 );



    double delta ( 0 );


    double step_lenght ( 0.000001 );
    double last_entr ( 0 );
    double min_entr ( 0 );
    double p_e_s;
    double * eigenvalues_plus_delta        = new double [size];

    do
        {
        this->calculate_repesentation();

        last_entr = state_rel_entropy ( ref,representation,in1_c,in2_c,temp_cmat_1,temp_cmat_2,temp_cmat_3,cvec,eigen_1,eigen_2,eigen_temp,size );

        vector_copy ( eigenvalues,cvec,size );

        for ( int sub_vec = 1 ; sub_vec < size ; ++sub_vec )
            {
            // Konstruiert die Parameter Darstellung der einzelnen Simplex-vectoren. Die Parameter respekieren Implizit die nebenbedingung Sum x_i = 1

            p_e_s=partial_eigen_sum( 0,sub_vec );       // im simplex ist der i subvektor der vektor (X0,...,0,0,Xi,0,0,...,0). Alle vektoren spannen R^n-1 auf

            // intersection search --------------------------------------------------------------------------------- hier werden die Parameter x_0=a-sum und x_i=1-a-sum minimiert.



            eigenvalues[0]= x_0* ( 1-p_e_s );
            eigenvalues[sub_vec]= ( 1-x_0 ) * ( 1-p_e_s );





            eigenvalues_plus_delta [0] = y_0* ( 1-p_e_s );
            eigenvalues_plus_delta [sub_vec] = ( 1- y_0 ) * ( 1-p_e_s );



















            }

        }
    while ( last_entr - min_entr > tolerance ); // nach durchlaufen aller vectoren wird die abbruchbedingung gesetstet




    }

double jarlskog_matrix::partial_eigen_sum ( int a, int b )
    {
    double result;
    for ( int i = 0; i < size ; ++i )
        {
        if ( (i =! a) && (i =! b) )
            {
            result+= this->eigenvalues[i].real();
            }

        }
    return result;
    }


void parametrized_biseparable_state::minimize_relative_entropy_of_entanglement ( std::complex< double >* ref, double tolerance )
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

    for ( int term = 0; term < nr_factors ; ++term )
        {




        }


    delete[] in1_c;
    delete[] in2_c;
    delete[] temp_cmat_1;
    delete[] temp_cmat_2;
    delete[] temp_cmat_3;
    delete[] cvec;
    delete[] eigen_1;
    delete[] eigen_2;
    delete[] eigen_temp;
    }

void parametrized_biseparable_state::fill_random()
    {
    for ( int i = 0; i < nr_factors; ++i )
        {
        this->substates_A[i].set_random();
        this->substates_B[i].set_random();
        }
    }

















