#ifndef __OPTFUNC__
#define __OPTFUNC__

#include "quit_toolbox.h"
#include "quit_unitary.h"
#include "parametrized_state.h"
#include "/home/stylx/projects/add_libs/dlib-18.18/dlib/optimization.h"
#include "/home/stylx/projects/add_libs/dlib-18.18/dlib/matrix/matrix.h"

typedef dlib::matrix<double,0,1> col_vec;

void convert_col_arr ( const col_vec col, double * mat );

double state_quantum_discord_qubit_ancilla ( const std::complex< double >*const rho_tot, unitary_jarlskog& minimal_mesurement, std::complex< double >* minimal_dec_state, int* sub_spaces, int nr_spaces, int size_tot, int nr_of_starts );

double state_relative_entropy_of_discord_qubit_ancilla ( const std::complex< double >*const rho_tot, unitary_jarlskog& minimal_mesurement, std::complex< double >* minimal_dec_state, int* sub_spaces, int nr_spaces, int size_tot, int nr_of_starts );

double state_opt_trace_distance_qubit_ancilla( const std::complex< double > * const rho_tot, unitary_jarlskog &minimal_mesurement,std::complex< double > * minimal_dec_state, int * sub_spaces, int nr_spaces, int size_tot );

double state_geometric_discord_qubit_ancilla ( const std::complex< double > * const rho_tot, unitary_jarlskog &minimal_mesurement, std::complex< double > * minimal_dec_state, int * sub_spaces, int nr_spaces, int size_tot, int nr_of_starts );

double relative_entropy_of_entanglement( std::complex< double >  * rho, param_state_multipartite_mixed& sigma, int nr_trials);



class opt_conditional_entropy_qubit_ancilla
{

public: 
  
  std::complex< double > * rho;
  std::complex< double > * decohered_state;
  std::complex< double > * ancilla;
  std::complex< double > * decohered_ancilla;
  std::complex< double > * rho_no_ancilla;
  int size;
  int size_no_ancilla;
  
  int * dec_vec;
  int * dec_anc ;
  int nr_sub;
  unitary_jarlskog &von_neumann_set;
   
  opt_conditional_entropy_qubit_ancilla (const std::complex< double >* state, std::complex< double >* _ancilla, unitary_jarlskog& mesurement, int* deco_vec, int nr_subspaces, int dim_state );
  
  ~opt_conditional_entropy_qubit_ancilla()
  {

    delete[] decohered_state;
    delete[] decohered_ancilla;
    delete[] rho_no_ancilla;
  }
  
  double operator() ( const col_vec p ) const;
  
};


class opt_relative_entropy_qubit_ancilla
{

public: 
  
  std::complex< double > * rho;
  std::complex< double > * decohered_state;
  int size;
  int * dec_vec;
  int nr_sub;
  unitary_jarlskog& von_neumann_set;
   
  opt_relative_entropy_qubit_ancilla (const std::complex < double > * state,
				 unitary_jarlskog &mesurement,     
				 int * dec_vector,
				 int nr_subspaces,
				 int dim_state);
  
  ~opt_relative_entropy_qubit_ancilla()
  {
    delete[] rho;
    delete[] decohered_state;
    delete[] dec_vec;
  }
  
  double operator() ( const col_vec p ) const;
  
};
//-----------------------------------------------------------------------------------------------------------

class opt_hilbert_schmidt_distance
{
  
public:
  
  std::complex< double > * rho;
  std::complex< double > * decohered_state;
  int size;
  int * dec_vec;
  int nr_sub;
  unitary_jarlskog& von_neumann_set;
   
  opt_hilbert_schmidt_distance (const std::complex < double > * state,
				 unitary_jarlskog &mesurement,     
				 int * dec_vector,
				 int nr_subspaces,
				 int dim_state);
  
  ~opt_hilbert_schmidt_distance()
  {
    delete[] rho;
    delete[] decohered_state;
    delete[] dec_vec;
  }
  
  double operator() ( const col_vec p ) const;
 
};


class opt_relative_entropy_of_entanglement
{

public: 
  
  std::complex< double > * rho;
  param_state_multipartite_mixed& separable_state;
   
  opt_relative_entropy_of_entanglement (std::complex < double > * state,
					param_state_multipartite_mixed& sigma);
  
  ~opt_relative_entropy_of_entanglement()
  {
    
  }
  
  double operator() ( const col_vec p ) const;
};
  

class opt_trace_distance_qubit_ancilla
{

public: 
  
  std::complex< double > * rho;
  std::complex< double > * decohered_state;
  int size;
  int * subspaces;
  int nr_sub;
  unitary_jarlskog& von_neumann_set;
   
  opt_trace_distance_qubit_ancilla (const std::complex < double > * state,
				 unitary_jarlskog &mesurement,     
				 int * _subspaces,
				 int nr_subspaces,
				 int dim_state);
  
  ~opt_trace_distance_qubit_ancilla()
  {
    delete[] rho;
    delete[] subspaces;
    delete[] decohered_state;
  }
  
  double operator() ( const col_vec p ) const;
};
  

#endif // __OPTFUNC__
