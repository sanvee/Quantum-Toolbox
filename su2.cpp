#include "su2.h"



void su2::calculate_reprsentation_para_euler()
{
  double c0 = cos(para_euler[0]); 
  double s0 = sin(para_euler[0]);
  
  double c1 = cos(para_euler[1]); 
  double s1 = sin(para_euler[1]);
  
  double c2 = cos(para_euler[2]); 
  double s2 = sin(para_euler[2]);
  
  representation[0] = std::complex< double > ( c0 * c1 , c0 * s1 );
  representation[1] = std::complex< double > ( s0 * c2 , s0 * s2 );
  representation[2] =-std::conj(representation[1]);
  representation[3] = std::conj(representation[0]);
  
}