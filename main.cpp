//-lstdc++ -lm -lopenblas -lgsl -llapacke -std=c++11//
// use gcc.4.9 // /usr/bin/gcc-4.9

#include "quit_toolbox.h"
#include "quit_import.h"
#include "opt-func.h"
#include "quit_unitary.h"
#include "parametrized_state.h"

using namespace std;

std::complex < double > * _BB_A ( double a,double p);
std::complex< double >* _A_BB ( double a, double p );

std::complex< double > * _BB_A_M ( double p );
std::complex< double > * _A_BB_M ( double p );

complex< double > _a ( 1.0/6.0,0.0 );


int main ( int argc, char **argv )
    {
         std::complex< double >  eta [8*8]   =
        {
        _a,_0,_0,_0,_0,_0,_0,_a,
        _0,_a,_0,_0,_0,_0,_0,_0,
        _0,_0,_a,_0,_0,_0,_0,_0,
        _0,_0,_0,_0,_0,_0,_0,_0,
        _0,_0,_0,_0,_0,_0,_0,_0,
        _0,_0,_0,_0,_0,_a,_0,_0,
        _0,_0,_0,_0,_0,_0,_a,_0,
        _a,_0,_0,_0,_0,_0,_0,_a
        };
	
	  std::complex< double >  singlet [8*8]   =
        {
        _a,_0,_0,_0,_0,_0,_0,_a,
        _0,_a,_0,_0,_0,_0,_0,_0,
        _0,_0,_a,_0,_0,_0,_0,_0,
        _0,_0,_0,_0,_0,_0,_0,_0,
        _0,_0,_0,_0,_0,_0,_0,_0,
        _0,_0,_0,_0,_0,_a,_0,_0,
       _0,_0,_0,_0,_0,_0,_a,_0,
        _a,_0,_0,_0,_0,_0,_0,_a
        };


    std::ofstream TEXTFILE;
    TEXTFILE.open ( "../depol_disc_neg.txt" ); //,std::ios_base::app );
    std::complex< double > * rho;

    int __1[2] = {8};
    int _4_2[2] = {4,-2};
    int _2_4[2] = {-2,4};
    int _2_2_2[3] = {2,2,2};
    
    
    std::complex< double > dec_state [8*8]; 
    std::complex< double > *test = new std::complex< double > [8*8](); 
    std::complex< double > * abb; 
    std::complex< double > * bba;  
    unitary_jarlskog u2(2);
    
   // param_state_multipartite_mixed _2_2_2_(3,_2_2_2,3);
   // param_state_multipartite_mixed _2_4_(2,_2_4,3);
  //  param_state_multipartite_mixed _4_2_(2,_4_2,3);
   
    
    int n = 100;
    
    for ( int i = 0 ; i < n; ++i )
    {
      double p = 1-1.0*(double)i/(double)n;
      
      TEXTFILE << p << '\t'; 
      
      cout << i << endl;
      
      //abb = _A_BB_M(p*M_PI/2);
      
    //  bba = _BB_A_M(p*M_PI/2);
      
      state_depolarizing(eta,test,p,8);
	  //matrix_show(test,8);
	//TEXTFILE << (double)i/n << '\t';
	TEXTFILE << state_quantum_discord_qubit_ancilla(test,u2,dec_state,_4_2,2,8,100) << '\t';
      
      
      
    //  cout << matrix_trace_nxn(bba,8) << endl;
    //  cout << matrix_trace_nxn(abb,8) << endl;
      
  //    matrix_show(abb,8);
   //   matrix_show(bba,8);

  //   TEXTFILE << state_quantum_discord_qubit_ancilla(abb,u2,dec_state,_2_4,2,8,100) << '\t';
     TEXTFILE << state_relative_entropy_of_discord_qubit_ancilla(test,u2,dec_state,_2_4,2,8,50) << '\t';
    // TEXTFILE << relative_entropy_of_entanglement(bba,_2_4_,50)-relative_entropy_of_entanglement(abb,_4_2_,50)<< '\t';
      
    // TEXTFILE << state_geometric_discord_qubit_ancilla(abb,u2,dec_state,_2_4,2,8,5000)<<'\t';
    // TEXTFILE << state_quantum_commutance(abb,2,4,8) <<'\t';
      
     TEXTFILE << state_logarithmic_negativity(test,_2_4,2,8) - state_logarithmic_negativity(test,_4_2,2,8) << endl;
     
   //  delete[] abb;
    //  delete[] bba;
        }

    TEXTFILE.close();
    gsl_rng_free ( _R_G );

    return 0;
    }

std::complex< double > * _A_BB ( double a, double p )
    {
    double c = cos ( a );
    double s = sin ( a );


    std::complex < double > * A_BB = new std::complex < double > [8*8]
        {
        0.25+c*c*_1,         _0,         _0,-0.25*_i+c*s,       -0.25*_i,          _0,         _0,    -0.25*_1,

                 _0,     0.5*_1,         _0,          _0,             _0,      0.5*_1,         _0,          _0,

                 _0,         _0,         _0,          _0,             _0,          _0,         _0,          _0,       

        c*s+0.25*_i,         _0,         _0, s*s+0.25*_1,       0.25*_1,          _0,         _0,    -0.25*_i,


            0.25*_i,         _0,         _0,     0.25*_1,        0.25*_1,          _0,         _0,    -0.25*_i,

                 _0,     0.5*_1,         _0,          _0,             _0,      0.5*_1,         _0,          _0,

                 _0,         _0,         _0,          _0,             _0,          _0,         _0,          _0,

           -0.25*_1,         _0,         _0,     0.25*_i,        0.25*_i,          _0,         _0,     0.25*_1
        };

    matrix_scalar_mult ( A_BB,1.0/(2.0+p),8 );

    return A_BB;
    }

std::complex< double > * _BB_A ( double a ,double p)
    {
    double c = cos ( a );
    double s = sin ( a );


    std::complex < double > * BB_A = new std::complex < double > [8*8]
        {
        0.25+p*c*c*_1,   -0.25*_i,         _0,          _0,             _0,          _0,-0.25*_i+p*c*s*_1,    -0.25*_1,

            0.25*_i,    0.25*_1,         _0,          _0,             _0,          _0,        0.25*_1,    -0.25*_i,

                 _0,         _0,     0.5*_1,      0.5*_1,             _0,          _0,             _0,          _0,

                 _0,         _0,     0.5*_1,      0.5*_1,             _0,          _0,             _0,          _0,


                 _0,         _0,         _0,          _0,             _0,          _0,             _0,          _0,

                 _0,         _0,         _0,          _0,             _0,          _0,             _0,          _0,

     0.25*_i+p*c*s*_1,    0.25*_1,         _0,          _0,             _0,          _0,    0.25+p*s*s*_1,    -0.25*_i,

           -0.25*_1,    0.25*_i,         _0,          _0,             _0,          _0,        0.25*_i,     0.25*_1
        };

    matrix_scalar_mult ( BB_A,1.0/(2.0+p),8 );

    return BB_A;

    }
    
    std::complex< double > * _A_BB_M ( double p )
    {
    double c = cos ( p );
    double s = sin ( p );
      
    std::complex < double > * A_BB = new std::complex < double > [8*8]
        {
     0.5*c*c+0.5*_1,         _0,          _0,-0.5*c*c*_i+0.5,    -0.5*c*s*_i,          _0,         _0,    -0.5*c*s*_1,

                _0,      c*c*_1,          _0,             _0,             _0,      c*s*_1,         _0,          _0, 

                _0,          _0,          _0,             _0,             _0,          _0,         _0,          _0,       

    0.5+0.5*c*c*_i,          _0,          _0, 0.5+0.5*c*c*_1,     0.5*s*c*_1,          _0,         _0,    -0.5*c*s*_i,


        0.5*s*c*_i,          _0,          _0,     0.5*s*c*_1,     0.5*s*s*_1,          _0,         _0,    -0.5*s*s*_i,

                _0,      c*s*_1,          _0,             _0,             _0,      s*s*_1,         _0,          _0,

                _0,          _0,          _0,             _0,             _0,          _0,         _0,          _0,

       -0.5*s*c*_1,          _0,          _0,      0.5*c*s*_i,    0.5*s*s*_i,          _0,         _0,     0.5*s*s*_1
	     
	};

    matrix_scalar_mult ( A_BB,1.0/3.0,8 );

    return A_BB;
    }

std::complex< double > * _BB_A_M ( double p )
    {
      
    double c = cos ( p );
    double s = sin ( p );

    std::complex < double > * BB_A = new std::complex < double > [8*8]
        { 
    0.5+0.5*c*c*_1, -0.5*s*c*_i,         _0,          _0,             _0,          _0,-0.5*c*c*_i+0.5*_1,    -0.5*s*c*_1,

        0.5*s*c*_i,  0.5*s*s*_1,         _0,          _0,             _0,          _0,        0.5*s*c*_1,    -0.5*s*s*_i,

                 _0,         _0,     c*c*_1,      c*s*_1,             _0,          _0,             _0,          _0,

                 _0,         _0,     c*s*_1,      s*s*_1,             _0,          _0,             _0,          _0,


                 _0,         _0,         _0,          _0,             _0,          _0,           _0,          _0,

                 _0,         _0,         _0,          _0,             _0,          _0,           _0,          _0,

  0.5*c*c*_i+0.5*_1,    0.5*s*c*_1,         _0,          _0,             _0,          _0,   0.5+0.5*c*c*_1,   -0.5*s*c*_i,

        -0.5*s*c*_1,    0.5*s*s*_i,         _0,          _0,             _0,          _0,       0.5*s*c*_i,    0.5*s*s*_1
        };

    matrix_scalar_mult ( BB_A,1.0/3.0,8 );

    return BB_A;

    }



