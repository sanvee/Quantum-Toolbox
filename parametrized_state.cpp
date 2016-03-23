#include"parametrized_state.h"

void param_state::calculate_simplex ( void )
    {
    for ( int i = 0; i < dim ; ++i )
        {
        this->eigen[i] = 1.0 ;
        }

    double s ( 1.0 ),c ( 1.0 ),n ( 1.0 );

    for ( int i = 0; i < dim-1 ; ++i )
        {

        c = cos ( simplex_angles[i] );
        s = sin ( simplex_angles[i] );
        n=1/ ( c+s );

        eigen[i] *= n*c;

        for ( int k = i+1 ; k < dim ; ++k )
            {

            eigen[k]*=n*s;

            }
        }
    }

void param_state::calculate_representation ( void )
    {
    
    set_paramters();
    calculate_simplex();  
    UN.calculate_representation();
    matrix_initialize_zero ( representation,dim );
    double a=0;
    for ( int i = 0 ; i < dim ; ++i )
        {
        representation[i* ( dim+1 )]=eigen[i];
        }

    //matrix_show ( representation,dim );

    double  _1_[2] = {1,0};
    double  _0_[2] = {0,0};

  matrix_initialize_zero(temp1,dim);
    
    cblas_zgemm ( CblasRowMajor,CblasNoTrans,CblasNoTrans,dim,dim,dim,_1_,UN.representation,dim,representation,dim,_0_,temp1,dim );
    cblas_zgemm ( CblasRowMajor,CblasNoTrans,CblasConjTrans,dim,dim,dim,_1_,temp1,dim,UN.representation,dim,_0_,representation,dim );
    
    }

void param_state::set_paramters ( void )
    {
      
    for ( int i = 0; i < dim * dim ; ++i )
        {
        UN.parameters[i] = parameters[i];
        }
    for ( int i = dim * dim, a = 0; i < dim * dim + dim -1 ; ++i, ++a )
        {
        simplex_angles[a] = parameters[i];
        }
    }


void param_state::set_limits()
    {
      
       for ( int i = 0 ; i < dim ; ++i )
        {
        par_max[i]= 2.0*M_PI; //<-- phases of unitary
        par_min[i] = 0;
        }
        
       for ( int i = dim; i < dim + dim - 1; ++i)
        {
        par_max[i]= M_PI/2.0; // <-- norms of the vectors.
        par_min[i] = 0;
        }
        
        for ( int i = dim + dim -1 ; i < dim * dim ; ++i )
        {
        par_max[i] = 2.0*M_PI; //<-----------------angles of the C^n
        par_min[i]=0;
        }
     
    for ( int i = dim * dim ; i < nr_parameters; ++i )
        {
        par_max[i] = M_PI/2.0; // simplex_angles[a] = parameters[i];
        par_min[i]=0;
        }
    }

void param_state_multipartite::set_limits()
    {
    int offset = 0;
    for ( int i = 0; i < nr_subspaces; ++i )
        {
        for ( int j = 0; j < states[i].nr_parameters; ++j )
            {
            par_max[j+offset] = states[i].par_max[j];//states[i].parameters[j] = parameters[j+offset];
            par_min[j+offset] = 0;
            }
        offset += this->states[i].nr_parameters;
        states[i].set_paramters();
        }
    }

void param_state_multipartite::set_paramters ( void )
    {
    int offset = 0;
    for ( int i = 0; i < nr_subspaces; ++i )
        {
        for ( int j = 0; j < states[i].nr_parameters; ++j )
            {
            states[i].parameters[j] = parameters[j+offset];
            }
        offset += this->states[i].nr_parameters;
        states[i].set_paramters();
        }
    }

void param_state_multipartite::calculate_representation ( void ) //dicht
    {
     
    set_paramters();
    for ( int i = 0 ; i < nr_subspaces ; ++i )
        {
        states[i].calculate_representation();
        }

    int curr_dim=subspaces[nr_subspaces-1];

    matrix_copy ( states[nr_subspaces-1].representation,temp1,subspaces[nr_subspaces-1] );

    for ( int i = nr_subspaces-2 ; i >= 0 ; --i )
        {

        matrix_tensor_prod ( states[i].representation,subspaces[i],temp1,curr_dim ,representation,subspaces[i]*curr_dim );
        curr_dim *= subspaces[i];
        matrix_copy ( representation,temp1,curr_dim );
        }
    }

void param_state::initialize ( int n )

    {
    dim = n ;
    nr_parameters = n*n+n-1;
    parameters = new double[nr_parameters] {};
    representation = new std::complex< double >[n*n] {};
    temp1 = new std::complex< double >[n*n] {};
    simplex_angles = new double [n-1] {};
    eigen = new std::complex< double >[n] {};
    par_max = new double [nr_parameters] {};
    par_min = new double [nr_parameters] {};
    UN.initialize ( n );
    set_limits();

    }

void param_state_multipartite::initialize ( int n, int* subs )
    {

    nr_subspaces = n ;
    subspaces = new int[n] {};
    dim = 1;
    for ( int i = 0; i < n ; i++ )
        {
        nr_parameters += subs[i]*subs[i]+abs ( subs[i] )-1;
        dim *= abs ( subs[i] );
        subspaces[i] = abs ( subs[i] );
        }


    states = new param_state[n] {};
    for ( int i = 0; i < n ; i++ )
        {
        states[i].initialize ( abs ( subs[i] ) );
        }

    parameters = new double [nr_parameters] {};
    par_max = new double [nr_parameters] {};
    par_min = new double [nr_parameters] {};
    representation = new std::complex< double > [dim*dim] {};
    temp1 = new std::complex< double > [dim*dim] {};
    buffer = new std::complex< double > [dim*dim] {};
    set_limits();

    }

void param_state_multipartite_mixed::initialize ( int n, int* subs, int d )
    {
    nr_mixtures = d;
    nr_subspaces = n ;
    subspaces = new int[n] {};
    dim = 1;
    nr_parameters=0;
    for ( int i = 0; i < n ; i++ )
        {
        nr_parameters += subs[i]*subs[i]+abs ( subs[i] )-1;
        dim *= abs ( subs[i] );
        subspaces[i] = abs ( subs[i] );
        }

    nr_parameters *=nr_mixtures;
    nr_parameters += d-1;

    product_states = new param_state_multipartite[d] {};

    for ( int i = 0; i < d ; i++ )
        {
        product_states[i].initialize ( n , subs );
        }

    parameters = new double [nr_parameters] {};
    par_max = new double [nr_parameters] {};
    par_min = new double [nr_parameters] {};
    representation = new std::complex< double > [dim*dim] {};
    temp1 = new std::complex< double > [dim*dim] {};
    buffer = new std::complex< double > [dim*dim] {};
    simplectic_prob = new double [d] {};
    theta_simplex = new double [d-1] {};
    set_limits();

    }
void param_state_multipartite_mixed::set_limits()
    {
        {
        int offset ( 0 );
        for ( int i = 0; i < nr_mixtures; ++i )
            {
            for ( int j = 0; j < product_states[i].nr_parameters; ++j )
                {
                par_max[j+offset]=product_states[i].par_max[j]; //product_states[i].parameters[j]=parameters[j+offset];
                par_min[j+offset]=0;

                }
            offset += product_states[i].nr_parameters;
            }
        for ( int j = offset; j < nr_parameters ; j++ )
            {
            par_max[j]=M_PI/2.0;//theta_simplex[j]=parameters[j];
            par_min[j]=0;
            }
        }
    }

void param_state_multipartite_mixed::set_paramters ( void )
    {
    int offset = 0;

    for ( int i = 0; i < nr_mixtures; ++i )
        {
        for ( int j = 0; j < product_states[i].nr_parameters; ++j )
            {
            product_states[i].parameters[j] = parameters[j+offset];
            }
        offset += product_states[i].nr_parameters;
        }

    for ( int j = offset , a=0; j < nr_parameters ; j++,a++ )
        {
        theta_simplex[a]=parameters[j];
        }
    }

void param_state_multipartite_mixed::calculate_simplex ()
    {

    for ( int i = 0; i < nr_mixtures ; ++i )
        {
        simplectic_prob[i] = 1.0 ;
        }

    double s ( 1.0 ),c ( 1.0 ),n ( 1.0 );

    for ( int i = 0; i < nr_mixtures-1 ; ++i )
        {

        c = cos ( theta_simplex[i] );
        s = sin ( theta_simplex[i] );
        n=1/ ( c+s );
        simplectic_prob[i] *=n ;
        simplectic_prob[i] *=c ;

        for ( int k = i+1 ; k < nr_mixtures ; ++k )
            {
            simplectic_prob[k] *= n*s;
            }
        }
    }

void param_state_multipartite_mixed::calculate_representation ( void )
    {

    set_paramters();

    calculate_simplex();

    for ( int i = 0 ; i < nr_mixtures ; ++ i )
        {
        product_states[i].calculate_representation();
        }
    matrix_initialize_zero ( representation,dim );

    for ( int i = 0; i < nr_mixtures ; ++i )
        {
        matrix_add_to_first ( representation, product_states[i].representation,simplectic_prob[i],dim );
        }
    }
















