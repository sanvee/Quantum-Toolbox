#ifndef __QUIT_JARLSKOG__
#define __QUIT_JARLSKOG__

#include "quit_toolbox.h"


class jarlskog_matrix
{

public:

    int size;
    std::complex< double > * vectors_normalized;
    std::complex< double > * vectors_not_normalized;
    double * norms;
    std::complex< double > * eigenvalues;
    std::complex< double > * representation;

private:

    std::complex< double > * temp_cvec;
    std::complex< double > * temp_V;
    std::complex< double > * temp_V_dag;
    std::complex< double > * temp_representation;
    double c,s;

public:

    jarlskog_matrix() {

        size = 0;
        vectors_normalized=NULL;
        norms=NULL;
        eigenvalues;
        representation=NULL;
        temp_cvec=NULL;
        temp_V=NULL;
        temp_V_dag=NULL;
        temp_representation=NULL;

    };

    jarlskog_matrix ( int s ) :
        vectors_normalized ( new std::complex< double >[ ( s-1 ) * ( s-1 )]() ),
        norms ( new double [ ( s-1 )]() ),
        eigenvalues ( new std::complex< double >[ ( s )]() ),

        representation ( new std::complex< double >[s*s]() ),
        temp_cvec ( new std::complex< double >[s-1]() ),
        temp_V ( new std::complex< double >[s*s]() ),
        temp_V_dag ( new std::complex< double >[s*s]() ),
        temp_representation ( new std::complex< double >[s*s]() ),
        size ( s ) {};

    ~jarlskog_matrix() {

        delete[] vectors_normalized;
        delete[] norms;
        delete[] eigenvalues;
        delete[] representation;

        delete[] temp_cvec;
        delete[] temp_V;
        delete[] temp_V_dag;
        delete[] temp_representation;
    };
    //zuweisungs operator definieren !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   
    void calculate_normalized ( void );
    void calculate_repesentation ( void );
    void set_eigen ( std::complex< double > * eigen_cvec );
    void set_eigen_random ( void );
    void set_vectors_random ( void );
    void set_vectors ( int i, int j );
    void set_norms ( double * norms );
    void set_norms_random ();
    void set_random();
    void minimize_rel_entropy ( std::complex< double >* ref, double tolerance );

private:
 void j_expand ( int klein, int gross );
    void rec_expand ( int n );
    void create_V_n_n ( int n );
    double partial_eigen_sum ( int a, int b );
    double intersection_search_rel_ent ( int a, int b );

};

class parametrized_biseparable_state
{

public:

    int size;
    int size_A;
    int size_B;
    int nr_factors;
    double * factors;
    std::vector < jarlskog_matrix > substates_A;
    std::vector < jarlskog_matrix > substates_B;
    
    parametrized_biseparable_state()
    {
      
      size = 0;
      size_A = 0;
      size_B = 0;
      nr_factors = 0;
      factors = NULL;  
      
      
    }

    parametrized_biseparable_state (int sizea ,int sizeb) :
    size_A (sizea),
    size_B (sizeb),
    size (sizea*sizeb),
    nr_factors(size*size+1),
    factors( new double [nr_factors]()),
    substates_A (nr_factors , jarlskog_matrix (size_A)),
    substates_B (nr_factors , jarlskog_matrix (size_B))
    {};
  
    
    void minimize_relative_entropy_of_entanglement ( std::complex< double > * ref, double tolerance);
    
private:
  
  void fill_random ();
 
};


#endif // __QUIT_JARLSKOG__ 
