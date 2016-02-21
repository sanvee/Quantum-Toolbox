#include "opt-func.h"


void convert_col_arr ( const col_vec col, double * mat )
    {

    for ( int i=0; i< col.size() ; ++i )

        {
        mat[i] = col ( 0,i );
        }

    }

void convert_col_arr ( const col_vec col, std::complex < double > * mat )
    {

    for ( int i=0; i< col.size() ; ++i )

        {

        mat[i].real ( col ( 0,i ) );
        }

    }

double opt_relative_entropy_qubits::operator() ( const col_vec params ) const
    {

    std::complex< double > * minimum;
    double * param;

    minimum = new std::complex < double > [4];
    param = new  double [3];
    convert_col_arr ( params,param );
    qubit::calc_rep_angle ( minimum, param );
    //std::cout << state_rel_entropy ( target_state,minimum, 2 ) << std::endl;
    double result = state_rel_entropy ( target_state,minimum, 2 );

    delete[] minimum;
    delete[] param;

    return result;

    }

double opt_relative_entropy_qubits_single_var::operator() ( const double p ) const
    {

    std::complex< double > * minimum;
    double * params;

    minimum = new std::complex < double > [4];
    params = new  double [3];

    vector_copy ( parameters,params,3 );

    params[param_nr]=p;

    qubit::calc_rep_angle ( minimum, params );
    //vector_write_to_file( params, state_rel_entropy ( target_state,minimum, 2 ) , "../converg.txt",3);
    //std::cout << state_rel_entropy ( target_state,minimum, 2 ) << std::endl;

    double result = state_rel_entropy ( target_state,minimum, 2 );

    delete[] minimum;
    delete[] params;

    return result;

    }


//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


opt_conditional_entropy_qubit_ancilla::opt_conditional_entropy_qubit_ancilla ( const std::complex < double > * state,
        int * dec_vector,
        int nr_subspaces,
        int dim_state )
    {

    rho = ( new std::complex< double >[dim_state*dim_state]() );
    dec_vec  = ( new int[nr_subspaces]() );

    size = dim_state;
    matrix_copy ( state, rho, size );

    nr_sub = nr_subspaces;
    vector_copy ( dec_vector, dec_vec, nr_sub );


    };




double opt_conditional_entropy_qubit_ancilla::operator() ( const col_vec p ) const
    {

    su2 von_neumann_set;

    std::complex< double > * decohered_state = new std::complex< double > [size*size];
    std::complex< double > * rho_no_ancilla = new std::complex< double >[ ( size-2 ) * ( size-2 )];

    convert_col_arr ( p,von_neumann_set.para_euler );
    vector_show<double> ( von_neumann_set.para_euler,3 );

    von_neumann_set.calculate_reprsentation_para_euler();

    state_decohere_on_subspace ( rho,decohered_state,von_neumann_set.representation,size,dec_vec,nr_sub );

    matrix_partial_trace ( decohered_state,rho_no_ancilla,dec_vec,nr_sub, size );

    double result = state_v_n_entropy ( rho_no_ancilla,size-2 );

    delete[] decohered_state;
    delete[] rho_no_ancilla;

    return result;

    }

//-----------------opt relative entropy-------------------------------------------------------------------------


opt_relative_entropy_qubit_ancilla::opt_relative_entropy_qubit_ancilla ( const std::complex < double > * state,
        int * dec_vector,
        int nr_subspaces,
        int dim_state )
    {

    rho = ( new std::complex< double >[dim_state*dim_state]() );
    dec_vec  = ( new int[nr_subspaces]() );

    size = dim_state;
    matrix_copy ( state, rho, size );

    nr_sub = nr_subspaces;
    vector_copy ( dec_vector, dec_vec, nr_sub );


    };




double opt_relative_entropy_qubit_ancilla::operator() ( const col_vec p ) const
    {

    double result;
    su2 von_neumann_set;

    std::complex< double > * decohered_state = new std::complex< double > [size*size];

    convert_col_arr ( p,von_neumann_set.para_euler );
    vector_show<double> ( von_neumann_set.para_euler,3 );

    von_neumann_set.calculate_reprsentation_para_euler();

    state_decohere_on_subspace ( rho,decohered_state,von_neumann_set.representation,size,dec_vec,nr_sub );


    result = state_rel_entropy ( rho,decohered_state,size );
    
    std::cout << result << std::endl;
    
    delete[] decohered_state;

    return result;

    }







//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


double state_quantum_discord_qubit_ancilla ( const std::complex< double > * const rho_tot,
        su2 &minimal_mesurement,
        std::complex< double > * minimal_dec_state,
        int * sub_spaces,
        int nr_spaces,
        int size_tot )
    {

    double S_ancilla ( 0 ),S_tot ( 0 ),S_conditional ( 0 );

    int size_no_ancilla = size_tot/2;

    std::complex< double > * rho_ancilla = new std::complex< double > [4];
    std::complex< double > * rho_dec_minimum = new std::complex< double > [size_tot*size_tot];
    std::complex< double > * rho_resudual = new std::complex< double > [ size_no_ancilla * size_no_ancilla];

    vector_scalar_mult ( sub_spaces,-1.0,3 );
    matrix_partial_trace ( rho_tot,rho_ancilla,sub_spaces, nr_spaces, size_tot );
    matrix_show ( rho_ancilla,2 );
    vector_scalar_mult ( sub_spaces,-1.0,3 );

    S_ancilla = state_v_n_entropy ( rho_ancilla,2 );

    std::cout << "S_ancilla= " << S_ancilla << std::endl;

    S_tot= state_v_n_entropy ( rho_tot,size_tot );

    std::cout << "s_tot=" << S_tot << std::endl;

    col_vec max ( 3 );
    col_vec min ( 3 );
    col_vec ini ( 3 );

    max = M_PI/2 , 2*M_PI, 2*M_PI;
    min = 0,0,0;
    ini = 0.1,1,2;

    dlib::find_min_bobyqa ( opt_conditional_entropy_qubit_ancilla ( rho_tot,sub_spaces,nr_spaces,size_tot ),
                            ini,
                            9,
                            min,
                            max,
                            0.2,
                            1e-8,
                            500 );

    minimal_mesurement.para_euler[0]=ini ( 0,0 );
    minimal_mesurement.para_euler[1]=ini ( 0,1 );
    minimal_mesurement.para_euler[2]=ini ( 0,2 );

    minimal_mesurement.calculate_reprsentation_para_euler();

    std::cout << "1=" << ini ( 0,0 ) << std::endl;
    std::cout << "2=" << ini ( 0,1 ) << std::endl;
    std::cout << "3=" << ini ( 0,2 ) << std::endl;

    state_decohere_on_subspace ( rho_tot,minimal_dec_state,minimal_mesurement.representation,size_tot,sub_spaces,nr_spaces );

    matrix_show ( const_cast <std::complex< double >*> ( rho_tot ),size_tot );

    matrix_partial_trace ( minimal_dec_state,rho_resudual,sub_spaces,nr_spaces,size_tot );

    matrix_show ( minimal_dec_state,size_tot );
    matrix_show ( rho_resudual,size_no_ancilla );

    S_conditional = state_v_n_entropy ( rho_resudual,size_no_ancilla );

    std::cout << " S_cond=" << S_conditional << std::endl;

    matrix_show ( minimal_mesurement.representation,2 );

    double result = S_ancilla - S_tot + S_conditional;

    delete[] rho_ancilla;
    delete[] rho_dec_minimum;
    delete[] rho_resudual;

    return result;

    }


double state_relative_entropy_of_discord_qubit_ancilla ( const std::complex< double > * const rho_tot,
        su2 &minimal_mesurement,
        std::complex< double > * minimal_dec_state,
        int * sub_spaces,
        int nr_spaces,
        int size_tot )
    {

    double result ( 0 );

    col_vec max ( 3 );
    col_vec min ( 3 );
    col_vec ini ( 3 );

    max = M_PI/2 , 2*M_PI, 2*M_PI;
    min = 0,0,0;
    ini = 0.1,1,2;

    dlib::find_min_bobyqa ( opt_relative_entropy_qubit_ancilla ( rho_tot,sub_spaces,nr_spaces,size_tot ),
                            ini,
                            9,
                            min,
                            max,
                            0.2,
                            1e-8,
                            500 );

    minimal_mesurement.para_euler[0]=ini ( 0,0 );
    minimal_mesurement.para_euler[1]=ini ( 0,1 );
    minimal_mesurement.para_euler[2]=ini ( 0,2 );

    minimal_mesurement.calculate_reprsentation_para_euler();

    std::cout << "1=" << ini ( 0,0 ) << std::endl;
    std::cout << "2=" << ini ( 0,1 ) << std::endl;
    std::cout << "3=" << ini ( 0,2 ) << std::endl;

    state_decohere_on_subspace ( rho_tot,minimal_dec_state,minimal_mesurement.representation,size_tot,sub_spaces,nr_spaces );
    
    result = state_rel_entropy(rho_tot,minimal_dec_state,size_tot);


    return result;

    }









