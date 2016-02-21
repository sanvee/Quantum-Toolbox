#ifndef __QUIT_UNITARY__
#define __QUIT_UNITARY__

#include "quit_toolbox.h"


class unitary_jarlskog
{

public:

    int dim;
    
    double * parameters; 
    
    /*parameter order is :
     * 
     [ 0 ; n-1 ] phases ( n phases ) 
     [ n ; 2n-2] norms  ( n-1 norms for the Jarlskog vectors )
     [ 2n-1;n^2+n-2] spherical coordinates of the complex unit Jarskog vectors z_n: n-2 vectors of dimension 2 to n-1 + 1 "vector" of dim 1 = 1 -> n^2-n-2 parameters
     we have in total n^2 + n - 3 parmeters for n < 2 ! n = 2 -> 4 parameters!
    *
    */
    
    std::complex< double > * z_n; // <---
    double * theta_z_n;
    double * norms;
    double * phases;

    std::complex< double > * representation;

private:

    std::complex< double > * temp_cvec;
    std::complex< double > * temp_V;
    std::complex< double > * temp_V_dag;
    std::complex< double > * temp_representation;
    std::complex< double > temp;
    double c,s;

public:

    unitary_jarlskog() {

        dim = 0;
        z_n = NULL;
        norms=NULL;
	phases = NULL;
        representation=NULL;

        temp_cvec=NULL;
        temp_V=NULL;
        temp_V_dag=NULL;
        temp_representation=NULL;

    };

    unitary_jarlskog ( int s ) :
    
    
    
        z_n ( new std::complex< double >[ ( s*(s-1))]() ),    
        theta_z_n ( new double[ s*( s-1 )-1]() ),
        norms ( new double [ ( s-1 )]() ),
        phases ( new double [ s ]() ),
        
        representation ( new std::complex< double >[s*s]() ),
        
        parameters ( new double [s*(s+1)-2]() ),
        
        temp_cvec ( new std::complex< double >[s-1]() ),
        temp_V ( new std::complex< double >[s*s]() ),
        temp_V_dag ( new std::complex< double >[s*s]() ),
        temp_representation ( new std::complex< double >[s*s]() ),
        dim ( s )
	{};

    ~unitary_jarlskog() {

        delete[] z_n;
        delete[] norms;
        delete[] representation;

        delete[] temp_cvec;
        delete[] temp_V;
        delete[] temp_V_dag;
        delete[] temp_representation;
    };

    void calculate_repesentation ( void );
    void set_z_n ( void );
    void set_paramters ( double * p );
   // void set_vectors ( int i, int j );
    void set_norms ( double * norms );
    void set_norms_random ();
    void set_random();
    void calculate_n_sphere_vectors();
   

private:
  
    void j_expand ( int klein, int gross );
    void create_V_n_n ( int n );

};







#endif // __QUIT_UNITARY__