#ifndef __PARAMETRIZED_STATE__
#define __PARAMETRIZED_STATE__

#include "quit_toolbox.h"
#include "quit_unitary.h"

class param_state
{

public:

    int dim;
    int nr_parameters;
    double * parameters;
    double * par_max;
    double * par_min;
    double * simplex_angles;

    std::complex< double > * eigen;
    std::complex< double > * representation;
    std::complex< double > * temp1;

    unitary_jarlskog UN;

    param_state() :
        dim ( 0 ),
        nr_parameters ( 0 ),
        UN() {
        parameters = NULL;
        eigen = NULL;
	par_max = NULL;
	par_min = NULL;
        simplex_angles = NULL;
        representation = NULL;
        temp1 = NULL;
    }

    param_state ( int n ) {
        initialize ( n );
    };

    ~param_state() {
        delete[] eigen;
        delete[] parameters;
        delete[] representation;
	delete[] par_max;
	delete[] par_min;
        delete[] temp1;
        delete[] simplex_angles;

    }

    void initialize ( int size );
    void set_limits ();
    void set_paramters ( void );
    void calculate_simplex ( void );
    void calculate_representation ( void );

};

class param_state_multipartite
{
public:

    int nr_subspaces;
    int nr_parameters;
    double * par_max;
    double * par_min;
    int dim;

    int * subspaces;
    double * parameters;
    param_state * states;
    std::complex< double > * representation;
    std::complex< double > * buffer;

private:

    std::complex< double > * temp1;


public:

    param_state_multipartite () :
        nr_subspaces ( 0 ),
        nr_parameters ( 0 ),
        dim ( 0 ) {
        subspaces = NULL;
        states = NULL;
        parameters = NULL;
	par_max = NULL;
	par_min = NULL;
        representation = NULL;
        temp1 = NULL;
        buffer = NULL;
    };

    ~param_state_multipartite () {

        delete[] subspaces;
        delete[] parameters;
	delete[] par_max;
	delete[] par_min;
        delete[] representation;
        delete[] states;
        delete[] temp1;
        delete[] buffer;

    };

    param_state_multipartite ( int n, int * subs )

    {
        initialize ( n,subs );
    };

    void initialize ( int n, int * subs );
    void set_limits ();
    void set_paramters ( void );
    void calculate_representation ( void );
    
};

class param_state_multipartite_mixed
{

public:

    int nr_subspaces;
    int nr_parameters;
    int dim;
    double * par_max;
    double * par_min;
    int nr_mixtures;

    int * subspaces;
    double * parameters;
    double * simplectic_prob;
    double * theta_simplex;
    
    param_state_multipartite * product_states;

    std::complex< double > * representation;
    std::complex< double > * buffer;

private:

    std::complex< double > * temp1;

public:

    param_state_multipartite_mixed () :
        nr_subspaces ( 0 ),
        nr_parameters ( 0 ),
        nr_mixtures ( 0 ),
        dim ( 0 ) {
        subspaces = NULL;
        product_states = NULL;
	theta_simplex = NULL;
        parameters = NULL;
	par_max = NULL;
	par_min = NULL;
        representation = NULL;
        temp1 = NULL;
        buffer = NULL;
	simplectic_prob = NULL;
    };

    param_state_multipartite_mixed ( int n , int * subs , int d ) // n = dimension d = nr of nr_mixtures
    {
      
        initialize ( n, subs ,d );
    }


    ~param_state_multipartite_mixed () {

        delete[] subspaces;
        delete[] parameters;
        delete[] representation;
	delete[] par_max;
	delete[] par_min;
        delete[] temp1;
        delete[] buffer;
	delete[] simplectic_prob;
	delete[] product_states;
	delete[] theta_simplex;
    };
    
    void initialize ( int n , int * subs , int d );
    void set_paramters ( void );
    void set_limits ();
    void calculate_simplex (void);
    void calculate_representation ( void );

};


#endif  // __PARAMETRIZED_STATE__
