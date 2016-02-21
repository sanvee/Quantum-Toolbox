#ifndef __OPTFUNC__
#define __OPTFUNC__

#include "quit_toolbox.h"
#include "qubit.h"
#include "su2.h"
#include "/home/stylx/projects/add_libs/dlib-18.18/dlib/optimization.h"
#include "/home/stylx/projects/add_libs/dlib-18.18/dlib/matrix/matrix.h"

typedef dlib::matrix<double,0,1> col_vec;

void convert_col_arr ( const col_vec col, std::complex< double > * mat );
void convert_col_arr ( const col_vec col, double * mat );

double state_quantum_discord_qubit_ancilla ( const std::complex< double >*const rho_tot,
					     su2& minimal_mesurement,
					     std::complex< double >* minimal_dec_state,
					     int* sub_spaces,
					     int nr_spaces,
					     int size_tot );

double state_relative_entropy_of_discord_qubit_ancilla ( const std::complex< double > * const rho_tot,
        su2 &minimal_mesurement,
        std::complex< double > * minimal_dec_state,
        int * sub_spaces,
        int nr_spaces,
        int size_tot );



class opt_v_n_entropy
{

public:

    opt_v_n_entropy ( ) { // test müll

    }
    double operator() ( const double arg ) const {
      
      double a = 0.5 * ( 1+arg );
      double b = 0.5 * ( 1-arg );

        //std::cout << arg << std::endl;

        return ( -a*log2 ( a )-b*log2 ( b ) );

    }

};

class opt_relative_entropy_qubits
{

public:

    opt_relative_entropy_qubits ( const std::complex< double > * const in ) {
      
        matrix_copy (in, target_state, 2); 
      
    }

    double operator() ( const col_vec params ) const;



private:

    std::complex< double >  target_state[4];

};

class opt_relative_entropy_qubits_single_var
{

public:

    opt_relative_entropy_qubits_single_var ( std::complex< double > * in, double * p, int nr ) {
        matrix_copy ( in, target_state, 2 );
        param_nr = nr;
        vector_copy ( p,parameters,3 );
    }

    double operator() ( const double p ) const;



private:

    std::complex< double >  target_state[4];
    int param_nr;
    double parameters[3] = {0.0,0.0,0.0};

};

class opt_conditional_entropy_qubit_ancilla
{

public: 
  
  std::complex< double > * rho;
  int size;
  int * dec_vec;
  int nr_sub;
   
  opt_conditional_entropy_qubit_ancilla (const std::complex < double > * state,
				 int * dec_vector,
				 int nr_subspaces,
				 int dim_state);
  
  ~opt_conditional_entropy_qubit_ancilla()
  {
    delete[] rho;
    delete[] dec_vec;
  }
  
  double operator() ( const col_vec p ) const;
 
  
  
  
};

class opt_relative_entropy_qubit_ancilla
{

public: 
  
  std::complex< double > * rho;
  int size;
  int * dec_vec;
  int nr_sub;
   
  opt_relative_entropy_qubit_ancilla (const std::complex < double > * state,
				 int * dec_vector,
				 int nr_subspaces,
				 int dim_state);
  
  ~opt_relative_entropy_qubit_ancilla()
  {
    delete[] rho;
    delete[] dec_vec;
  }
  
  double operator() ( const col_vec p ) const;
 
  
  
  
};




#endif // __OPTFUNC__
