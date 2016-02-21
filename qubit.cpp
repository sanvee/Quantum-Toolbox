#include "quit_toolbox.h"
#include "qubit.h"



const std::complex< double > imag_u = std::complex< double > (0,1);

std::complex< double > sigma_x [4]={  std::complex< double > (0,0) , std::complex< double > (1,0),
				      std::complex< double > (1,0) , std::complex< double > (0,0)};
					
std::complex< double > sigma_y [4]={  std::complex< double > (0,0) , std::complex< double > (0,-1),
				      std::complex< double > (0,1) , std::complex< double > (0,0)};

std::complex< double > sigma_z [4]={  std::complex< double > (1,0) , std::complex< double > (0,0),
				      std::complex< double > (0,0) , std::complex< double > (-1,0)};	    

void qubit::calc_rep_bloch()
    {
      representation[0].real(0.5 * (1 + bloch_v[2]));
      
      representation[1].real(0.5 * bloch_v[0]);
      representation[1].imag(-0.5 * bloch_v[1]);
      
      representation[2]=std::conj(representation[1]);
      
      representation[3].real(0.5 * (1 - bloch_v[2]));

    }
    
void qubit::calc_rep_bloch( std::complex< double >* state, double * para)
    {
      state[0].real(0.5 * (1 +  para[2]));
      
      state[1].real(0.5 *  para[0]);
      state[1].imag(-0.5 *  para[1]);
      
      state[2]= std::conj(state[1]);
      state[3].real(0.5 * (1 -  para[2]));

    }
    
void qubit::calc_rep_angle()
    {
      
      bloch_v[0] = bloch_para[0]*sin(bloch_para[1])*cos(bloch_para[2]);
      bloch_v[1] = bloch_para[0]*sin(bloch_para[1])*sin(bloch_para[2]);
      bloch_v[2] = bloch_para[0]*cos(bloch_para[1]); 
      calc_rep_bloch();     

    }
void qubit::calc_rep_angle ( std::complex< double > * state,const double * const p )
    {
      double * param = new double[3];
      
      
      param[0]=p[0]*sin( p[1])*cos( p[2]);
      param[1]=p[0]*sin( p[1])*sin( p[2]);
      param[2]=p[0]*cos( p[1]);
      calc_rep_bloch(state,param);  
      
      delete[] param;

    }

    
    
double qubit::v_n_entropy_para ()
    {
      calc_rep_angle();
      return state_v_n_entropy(representation,repbuff,eigbuff,2);

    }
    
double qubit::v_n_entropy_para ( double * para )
    {
      

    }





