//-lstdc++ -lm -lopenblas -lgsl -llapacke -std=c++11//

#include "quit_toolbox.h"
#include "/home/stylx/projects/add_libs/scs-master/include/scs.h"
#include "quit_import.h"
#include <curses.h>
#include <regex>
#include <string>
#include <algorithm>
#include "opt-func.h"
#include "qubit.h"
#include "quit_jarlskog.h"
#include "quit_unitary.h"


std::complex< double > _a ( 1.0/6,0.0 );

int main ( int argc, char **argv )
    {

    std::ofstream TEXTFILE;
    TEXTFILE.open ( "../test_tab.txt" );

    std::ofstream FILE;
    FILE.open ( "../test_tab.txt" );
    std::complex< double > * IMPORTED_MATRIX;
    int SIZE_IMPORTED;
    IMPORTED_MATRIX = Mathematica_import_matrix_from_file ( "../import.txt",SIZE_IMPORTED );
    
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

	su2 eta_mimimal_mesurement;
	std::complex< double >  eta_minimal[8*8];
	int eta_spaces[3]={2,2,-2};
	
	std::cout << state_relative_entropy_of_discord_qubit_ancilla(eta,eta_mimimal_mesurement,eta_minimal,eta_spaces,3,8) << std::endl;
	
	matrix_show(eta_mimimal_mesurement.representation,2);
	
	matrix_show(eta_minimal,8);

	

    TEXTFILE.close();
    gsl_rng_free ( _R_G );

    return 0;
    }
