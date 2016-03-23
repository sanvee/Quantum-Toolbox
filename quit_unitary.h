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
     [ 0 ; n-1 ] phases ( n phases ) 0-2PI
     [ n ; 2n-2] norms  ( n-1 norms for the Jarlskog vectors ) in [0;pi/2]
     [ 2n-1;n^2+n-2] spherical coordinates of the complex unit Jarskog vectors z_n: n-2 vectors of dimension 2 to n-1 + 1 "vector" of dim 1 = 1 -> n^2-n-2 parameters
     we have in total n^2 + n - 3 parmeters for n < 2 ! n = 2 -> 4 parameters!
    *
    */

    std::complex< double > * z_n;	// <--- array with the complex vectors on the complex unit sphere.
    double * z_n_coefs;			// <---- temp array for z_n needed for calculations
    double * theta_z_n;
    double * norms;
    double * phases;

    std::complex< double > * representation;
    std::complex< double > * temp_representation;
    int nr_parameters,nr_theta,nr_phases,nr_norms,nr_zn;


private:

    std::complex< double > * temp_cvec;
    std::complex< double > * temp_V;
    std::complex< double > * temp_V_dag;
    std::complex< double > temp;
    double c,s;



public:

    unitary_jarlskog() {

        dim = 0;
	nr_parameters = 0;
	nr_phases = 0;
	nr_norms = 0;
	nr_theta = 0;
	nr_zn = 0;
	temp = 0;
	c= 0.0;
	s= 0.0;
        z_n = NULL;
        norms=NULL;
        phases = NULL;
        representation=NULL;
	parameters=NULL;
        temp_cvec=NULL;
        temp_V=NULL;
        temp_V_dag=NULL;
        temp_representation=NULL;
	z_n_coefs=NULL;
	theta_z_n=NULL;

    };

    unitary_jarlskog ( int s ) :
        dim ( s ),
        nr_parameters ( s * s ),
        nr_phases ( s ),
        nr_norms ( s-1 ),
        nr_theta ( ( s-1 ) * ( s-1 ) ),
        nr_zn ( s * ( s-1 ) /2 ),

        phases ( new double [s]{} ),
        z_n ( new std::complex< double >[s * ( s-1 ) /2]{} ),
        z_n_coefs ( new double [  s * ( s-1 )-2]{} ),
        theta_z_n ( new double[ ( s-1 ) * ( s-1 )]{} ),
        norms ( new double [s-1]{} ),

        representation ( new std::complex< double >[s*s]{} ),

        parameters ( new double [s * s]{} ),

        temp_cvec ( new std::complex< double >[s-1]{} ),
        temp_V ( new std::complex< double >[s*s]{} ),
        temp_V_dag ( new std::complex< double >[s*s]{} ),
        temp_representation ( new std::complex< double >[s*s]{} ) {
	  
    };

    ~unitary_jarlskog() {

        delete[] z_n;
        delete[] z_n_coefs;
        delete[] norms;
        delete[] theta_z_n;
        delete[] phases;
        delete[] representation;
        delete[] parameters;

        delete[] temp_cvec;
        delete[] temp_V;
        delete[] temp_V_dag;
        delete[] temp_representation;
    };
    
    
    
    
    void initialize (int s);
    void calculate_representation ( void );
    void set_z_n ( void );
    void set_paramters ();
    // void set_vectors ( int i, int j );
    void calculate_n_sphere_vectors();


private:

    void j_expand ( int klein, int gross );
    void create_V_n_n ( int n );

};







#endif // __QUIT_UNITARY__
