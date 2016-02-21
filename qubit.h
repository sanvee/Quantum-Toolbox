#ifndef __QUBIT__
#define __QUBIT__

#include "quit_toolbox.h"




class qubit
{


public:

    double bloch_v[3];
    double bloch_para[3];
    
    std::complex< double > eigenvalues [2];
    std::complex< double > representation [4];

private:

    std::complex< double > eigbuff [2];
    std::complex< double > repbuff [4];

public:

    qubit() {
      
        bloch_v[0]=0;
        bloch_v[1]=0;
        bloch_v[2]=0;

        bloch_para[0]=0;
        bloch_para[1]=0;
        bloch_para[2]=0;

        eigenvalues[0]=std::complex< double > ( 0,0 );
        eigenvalues[1]=std::complex< double > ( 0,0 );

        representation[0]=std::complex< double > ( 0,0 );
        representation[1]=std::complex< double > ( 0,0 );
        representation[2]=std::complex< double > ( 0,0 );
        representation[3]=std::complex< double > ( 0,0 );
	
	
	repbuff[0]=std::complex< double > ( 0,0 );
        repbuff[1]=std::complex< double > ( 0,0 );
        repbuff[2]=std::complex< double > ( 0,0 );
        repbuff[3]=std::complex< double > ( 0,0 );
	
	eigbuff[0]=std::complex< double > ( 0,0 );
        eigbuff[1]=std::complex< double > ( 0,0 );
	

    }
    
    
    
    void calc_rep_bloch();
    
    static void calc_rep_bloch(std::complex< double >* state, double * para );
    
    void calc_rep_angle();
    
    static void calc_rep_angle( std::complex< double >* state, const double*const p );
    
    double v_n_entropy_para();
    
    static double v_n_entropy_para( double* para );
    
};

#endif  // __QUIT_TOOLBOX__

