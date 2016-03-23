#include "opt-func.h"

using namespace std;

void convert_col_arr ( const col_vec col, double * mat )
    {

    for ( int i=0; i< col.size() ; ++i )

        {
        mat[i] = col ( 0,i );
        }

    }

void convert_col_arr ( double * mat, col_vec col )
    {

    for ( int i=0; i< col.size() ; ++i )

        {
        col ( 0,i ) = mat[i];
        std::cout << col ( 0,i ) << " ";
        }

    }


void convert_col_arr ( const col_vec col, std::complex < double > * mat )
    {

    for ( int i=0; i< col.size() ; ++i )

        {

        mat[i].real ( col ( 0,i ) );
        }

    }


//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


opt_conditional_entropy_qubit_ancilla::opt_conditional_entropy_qubit_ancilla (
	const std::complex < double > * state,
	std::complex < double > * _ancilla,								       
        unitary_jarlskog &mesurement,
        int * deco_vec,
        int nr_subspaces,
        int dim_state ) :
    von_neumann_set ( mesurement )
    {
      
    
    decohered_state = new std::complex< double > [dim_state*dim_state]();
    rho_no_ancilla = new std::complex< double >[ ( dim_state/2 ) * ( dim_state/2 )]();
    decohered_ancilla = new std::complex< double >[2*2]();
    
    dec_vec = deco_vec;
    dec_anc = new int[1]();
    dec_anc[0]=2;
    rho = const_cast <std::complex< double >*> ( state);
    ancilla = _ancilla;
    //dec_vec  =  new int[nr_subspaces]() ;

    size = dim_state;
    size_no_ancilla = dim_state/2;
    nr_sub = nr_subspaces;
    //vector_copy ( dec_vector, dec_vec, nr_sub );

    };




double opt_conditional_entropy_qubit_ancilla::operator() ( const col_vec p ) const
    {

    convert_col_arr ( p,von_neumann_set.parameters );
    //vector_show<double> ( von_neumann_set.parameters,3 );

    von_neumann_set.calculate_representation();
   //matrix_show(von_neumann_set.representation,2);
    state_decohere_on_subspace ( rho,decohered_state,von_neumann_set.representation,size,dec_vec,nr_sub );

    matrix_partial_trace ( decohered_state,rho_no_ancilla,dec_vec,nr_sub, size );
     
    vector_scalar_mult ( dec_vec,-1,nr_sub );
    matrix_partial_trace ( decohered_state,decohered_ancilla,dec_vec, nr_sub, size );                               // trace out ancilla
    vector_scalar_mult ( dec_vec,-1,nr_sub );
    double result = state_v_n_entropy ( decohered_state, size) - state_v_n_entropy( decohered_ancilla,2 );
    return result;
    }

//-----------------opt relative entropy-------------------------------------------------------------------------


opt_relative_entropy_qubit_ancilla::opt_relative_entropy_qubit_ancilla ( const std::complex < double > * state,
        unitary_jarlskog &mesurement,
        int * dec_vector,
        int nr_subspaces,
        int dim_state ) :
    von_neumann_set ( mesurement )
    {

    rho = ( new std::complex< double >[dim_state*dim_state]() );
    decohered_state = ( new std::complex< double >[dim_state*dim_state]() );
    dec_vec  = ( new int[nr_subspaces]() );

    size = dim_state;
    matrix_copy ( state, rho, dim_state );

    nr_sub = nr_subspaces;
    vector_copy ( dec_vector, dec_vec, nr_sub );
    

    };




double opt_relative_entropy_qubit_ancilla::operator() ( const col_vec p ) const
    {

    double result;

    convert_col_arr ( p,von_neumann_set.parameters );

    // vector_show<double> ( von_neumann_set.parameters,3 );
    von_neumann_set.set_paramters();
    von_neumann_set.calculate_representation();
    //  matrix_show(von_neumann_set.representation,2);

    state_decohere_on_subspace ( rho,decohered_state,von_neumann_set.representation,size,dec_vec,nr_sub );


    result = state_rel_entropy ( rho,decohered_state,size );

    //std::cout << result << std::endl;
    return result;

    }

//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

opt_trace_distance_qubit_ancilla::opt_trace_distance_qubit_ancilla ( const std::complex< double >* state,
        unitary_jarlskog& mesurement,
        int* _subspaces,
        int nr_subspaces,
        int dim_state ) :
    von_neumann_set ( mesurement )
    {

    rho = ( new std::complex< double >[dim_state*dim_state]() );
    subspaces  = ( new int[nr_subspaces]() );
    decohered_state = ( new std::complex< double >[dim_state*dim_state]() );

    size = dim_state;
    matrix_copy ( state, rho, dim_state );

    nr_sub = nr_subspaces;
    vector_copy ( _subspaces, subspaces, nr_sub );

    };

double opt_trace_distance_qubit_ancilla::operator() ( const col_vec p ) const
    {


    convert_col_arr ( p,von_neumann_set.parameters );
    von_neumann_set.set_paramters();
    von_neumann_set.calculate_representation();

    state_decohere_on_subspace ( rho,decohered_state,von_neumann_set.representation,size,subspaces,nr_sub );
    return trace_distance ( rho,decohered_state,size );


    }
//-----------------------------------------------------------------------------------------------

opt_relative_entropy_of_entanglement::opt_relative_entropy_of_entanglement (
    std::complex< double > * state,
    param_state_multipartite_mixed& sigma ) :
    separable_state ( sigma )
    {
    rho=state;
    };




double opt_relative_entropy_of_entanglement::operator() ( const col_vec p ) const
    {

    double result;

    convert_col_arr ( p,separable_state.parameters );

    separable_state.calculate_representation();

    return state_rel_entropy ( const_cast < std::complex< double >*> ( rho ), separable_state.representation, separable_state.dim );



    }

//-----------------------------------------------------------------------------------------------------------




//---------------------------------------------------------------------------------------------------------

opt_hilbert_schmidt_distance::opt_hilbert_schmidt_distance ( const std::complex< double >* state,
        unitary_jarlskog& mesurement,
        int* dec_vector,
        int nr_subspaces,
        int dim_state ) :
    von_neumann_set ( mesurement )
    {

    rho = ( new std::complex< double >[dim_state*dim_state]() );
    decohered_state = ( new std::complex< double >[dim_state*dim_state]() );

    dec_vec  = ( new int[nr_subspaces]() );

    size = dim_state;
    matrix_copy ( state, rho, dim_state );

    nr_sub = nr_subspaces;
    vector_copy ( dec_vector, dec_vec, nr_sub );

    };

double opt_hilbert_schmidt_distance::operator() ( const col_vec p ) const
    {

    convert_col_arr ( p,von_neumann_set.parameters );

    // vector_show<double> ( von_neumann_set.parameters,3 );
    von_neumann_set.set_paramters();
    von_neumann_set.calculate_representation();
    //  matrix_show(von_neumann_set.representation,2);

    state_decohere_on_subspace ( rho,decohered_state,von_neumann_set.representation,size,dec_vec,nr_sub );

    return state_hilbert_schmidt_distance ( rho,decohered_state,size );

    //std::cout << result << std::endl;

    }

//

double state_quantum_discord_qubit_ancilla ( const std::complex< double > * const rho_tot,
        unitary_jarlskog &minimal_mesurement,
        std::complex< double > * minimal_dec_state,
        int * sub_spaces,
        int nr_spaces,
        int size_tot,
        int nr_of_starts )
    {

    double S_ancilla = 0.0;
    double S_tot = 0.0 ;
    double S_conditional = std::numeric_limits<double>::infinity() ;
    double temp_result;
    int size_no_ancilla = size_tot/2;

    std::complex< double > * rho_ancilla = new std::complex< double > [4];
    std::complex< double > * rho_dec_minimum = new std::complex< double > [size_tot*size_tot];
    std::complex< double > * rho_resudual = new std::complex< double > [ size_no_ancilla * size_no_ancilla];
    
    S_tot= state_v_n_entropy ( rho_tot,size_tot );
    
    vector_scalar_mult ( sub_spaces,-1,nr_spaces );

    matrix_partial_trace ( rho_tot,rho_ancilla,sub_spaces, nr_spaces, size_tot );                               // trace out ancilla
    //matrix_show ( rho_ancilla,2 );
    vector_scalar_mult ( sub_spaces,-1,nr_spaces );

    S_ancilla = state_v_n_entropy ( rho_ancilla,2 );

    //std::cout << "S_ancilla= " << S_ancilla << std::endl;

    

    //std::cout << "s_tot=" << S_tot << std::endl;

    col_vec max ( 4 );
    col_vec min ( 4 );
    col_vec ini ( 4 );

    max = 2*M_PI , 2*M_PI, 2*M_PI, 2*M_PI;
    min = 0,0,0,0;
    try
        {
        for ( int i = 0; i < nr_of_starts; ++i )
            {
            ini = gsl_rng_uniform ( _R_G ) *2*M_PI,gsl_rng_uniform ( _R_G ) *M_PI*2,gsl_rng_uniform ( _R_G ) *2*M_PI,gsl_rng_uniform ( _R_G ) *2*M_PI;

           temp_result = dlib::find_min_bobyqa ( opt_conditional_entropy_qubit_ancilla ( rho_tot,rho_ancilla,minimal_mesurement,sub_spaces,nr_spaces,size_tot ),
                                    ini,
                                    9,
                                    min,
                                    max,
                                    0.1,
                                    1e-8,
                                    500 );

	    //std::cout << "s_condition " <<  temp_result << std::endl;
	 //   std::cout << "s_ancilla " << S_ancilla << std::endl;
	  //  std::cout << "s_tot " <<  S_tot << std::endl;

            if ( temp_result < S_conditional )
                {

               // std::cout << "hit_QD-"<< i << '\t' << "sconditional=" << S_conditional-temp_result << "Result=" << S_ancilla - S_tot + S_conditional  << std::endl;
                S_conditional = temp_result;

                }
            }
        }
    catch ( dlib::bobyqa_failure e )
        {
        S_conditional = std::numeric_limits<double>::infinity();
        }
        
      //  std::cout << "s_ancilla " << S_ancilla << std::endl;
     //   std::cout << "s_tot " <<  S_tot << std::endl;
      //  std::cout << "s_condition " <<  S_conditional << std::endl;
   
    double result = S_ancilla - S_tot + S_conditional;

    delete[] rho_ancilla;
    delete[] rho_dec_minimum;
    delete[] rho_resudual;

    return result;

    }


double state_geometric_discord_qubit_ancilla ( const std::complex< double > * const rho_tot,
        unitary_jarlskog &minimal_mesurement,
        std::complex< double > * minimal_dec_state,
        int * sub_spaces,
        int nr_spaces,
        int size_tot,
        int nr_of_starts )
    {

    double result ( std::numeric_limits<double>::infinity() );
    double temp_result ( 0 );

    col_vec max ( 4 );
    col_vec min ( 4 );
    col_vec ini ( 4 );

    max = 2*M_PI , 2*M_PI, 2*M_PI , 2*M_PI;
    min = 0,0,0,0;
    try
        {

        for ( int i = 0; i < nr_of_starts; ++i )

            {
            ini = gsl_rng_uniform ( _R_G ) *2*M_PI,gsl_rng_uniform ( _R_G ) *M_PI/2,gsl_rng_uniform ( _R_G ) *2*M_PI,gsl_rng_uniform ( _R_G ) *2*M_PI;

            dlib::find_min_bobyqa ( opt_hilbert_schmidt_distance ( rho_tot,minimal_mesurement,sub_spaces,nr_spaces,size_tot ),
                                    ini,
                                    9,
                                    min,
                                    max,
                                    0.2,
                                    1e-6,
                                    500 );

            state_decohere_on_subspace ( rho_tot,minimal_dec_state,minimal_mesurement.representation,size_tot,sub_spaces,nr_spaces );
            temp_result = state_hilbert_schmidt_distance ( const_cast<std::complex< double >*> ( rho_tot ),minimal_dec_state,size_tot );
            if ( temp_result < result )
                {
                 std::cout << "hit_GD-"<< i <<'\t'<< result << '\t' <<  temp_result <<std::endl;
                result = temp_result;

                }
            }
        }
    catch ( dlib::bobyqa_failure e )
        {
        result = temp_result;
        }

    return result;

    }

double state_relative_entropy_of_discord_qubit_ancilla ( const std::complex< double > * const rho_tot,
        unitary_jarlskog &minimal_mesurement,
        std::complex< double > * minimal_dec_state,
        int * sub_spaces,
        int nr_spaces,
        int size_tot,
        int nr_of_starts )
    {

    double result ( std::numeric_limits<double>::infinity() );
    double temp_result ( 0 );

    col_vec max ( 4 );
    col_vec min ( 4 );
    col_vec ini ( 4 );

    max = 2*M_PI , 2*M_PI, 2*M_PI , 2*M_PI;
    min = 0,0,0,0;
    try
        {

        for ( int i = 0; i < nr_of_starts; ++i )

            {
            ini = gsl_rng_uniform ( _R_G ) *2*M_PI,gsl_rng_uniform ( _R_G ) *M_PI/2,gsl_rng_uniform ( _R_G ) *2*M_PI,gsl_rng_uniform ( _R_G ) *2*M_PI;

            dlib::find_min_bobyqa ( opt_relative_entropy_qubit_ancilla ( rho_tot,minimal_mesurement,sub_spaces,nr_spaces,size_tot ),
                                    ini,
                                    9,
                                    min,
                                    max,
                                    0.2,
                                    1e-6,
                                    500 );

            state_decohere_on_subspace ( rho_tot,minimal_dec_state,minimal_mesurement.representation,size_tot,sub_spaces,nr_spaces );
            temp_result = state_rel_entropy ( rho_tot,minimal_dec_state,size_tot );
            if ( temp_result < result )
                {
                // std::cout << "hit_RE-"<< i <<'\t'<< result << '\t' <<  temp_result <<std::endl;
                result = temp_result;

                }
            }
        }
    catch ( dlib::bobyqa_failure e )
        {
        result = std::numeric_limits<double>::infinity();
        }

    return result;

    }

double relative_entropy_of_entanglement ( std::complex< double >  * rho, param_state_multipartite_mixed& sigma, int nr_trials )
    {
    col_vec max ( sigma.nr_parameters );
    col_vec min ( sigma.nr_parameters );
    col_vec ini ( sigma.nr_parameters );

    double result ( std::numeric_limits<double>::infinity() );
    double temp_result ( 0 );


    for ( int trials = 0; trials < nr_trials; ++trials )
        {

        try
            {
            for ( int i = 0; i < sigma.nr_parameters; ++i )
                {
                min ( 0,i ) = sigma.par_min[i];
                ini ( 0,i ) = sigma.par_max[i] * gsl_rng_uniform ( _R_G );
                max ( 0,i ) = sigma.par_max[i];
                //std::cout << min (0,i) << '\t' << ini(0,i) << '\t' << max (0,i) << std::endl;
                }

            temp_result = dlib::find_min_bobyqa ( opt_relative_entropy_of_entanglement ( rho,sigma ),
                                                  ini,
                                                  2*sigma.nr_parameters+1,
                                                  min,
                                                  max,
                                                  0.3,
                                                  1e-6,
                                                  5000 );
            }
        catch ( dlib::bobyqa_failure e )
            {
            //cout << "error ";
            }

        if ( temp_result < result )
            {
             std::cout << result << std::endl;
            result = temp_result;
            }
        }
    return result;
    }

double state_opt_trace_distance_qubit_ancilla ( const std::complex< double > * const rho_tot,
        unitary_jarlskog &minimal_mesurement,
        std::complex< double > * minimal_dec_state,
        int * sub_spaces,
        int nr_spaces,
        int size_tot )
    {
    double result ( 0 );

    col_vec max ( 4 );
    col_vec min ( 4 );
    col_vec ini ( 4 );

    max = 2*M_PI , 2*M_PI, 2*M_PI , 2*M_PI;
    min = 0,0,0,0;
    ini = gsl_rng_uniform ( _R_G ) *2*M_PI,gsl_rng_uniform ( _R_G ) *M_PI/2,gsl_rng_uniform ( _R_G ) *2*M_PI,gsl_rng_uniform ( _R_G ) *2*M_PI;

    dlib::find_min_bobyqa ( opt_trace_distance_qubit_ancilla ( rho_tot,minimal_mesurement,sub_spaces,nr_spaces,size_tot ),
                            ini,
                            9,
                            min,
                            max,
                            0.2,
                            1e-8,
                            500 );

    state_decohere_on_subspace ( rho_tot,minimal_dec_state,minimal_mesurement.representation,size_tot,sub_spaces,nr_spaces );

    result = trace_distance ( rho_tot,minimal_dec_state,size_tot );

    }
