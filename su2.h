#ifndef __SU2__
#define __SU2__

#include "quit_toolbox.h"

class su2
{

public:

    std::complex< double > representation [4];
    std::complex< double > para_caley_klein[2];
    double para_euler[3];   // Max PI/2, 2PI, 2PI.

    su2() : representation (), para_caley_klein (), para_euler()
    {
    };
    
    void calculate_reprsentation_para_euler(); // (theta,alpha,gamma) -->  (cos (theta)* exp ( i alpha ) --- e.c.t.


};



#endif // __SU2__
