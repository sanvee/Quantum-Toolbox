#include "quit_unitary.h"


void unitary_jarlskog::calculate_repesentation ( void ) //<----------------vervollständigen.
    {
    
    c = cos ( norms[0] );
    s = sin ( norms[0] );
    
    representation[0] = c;
    representation[1] = s * z_n[0];
    representation[2] =-s * z_n[0];
    representation[3] = c;

    for ( int n = 2; n < dim; ++n )
        {
	  
        j_expand ( n,n+1 );
    
        representation [n * ( n+2 )] = 1; // (n+1)*(n+1)-1 letzes element der matrix 

        create_V_n_n ( n );  // n = current size - 1
        matrix_mlt ( representation,temp_V,temp_representation,n+1 );
        matrix_copy ( temp_representation,representation,n+1 );
        }

    for ( int n = 0; n < dim; ++n )
        {
        temp = exp ( std::complex< double > ( 0 , phases[n] ) );

        for ( int j = 0 ; j < dim ; ++j )
            {

            representation[ j + dim * n] *= temp;

            }
        }
    }

void unitary_jarlskog::create_V_n_n ( int n )

    {

    matrix_initialize_zero ( temp_V,n+1 );
    
    
    //nix gut  -> matrix_get_partial_column ( n-1,n,z_n,temp_cvec,dim-1 );
    //besser -> get_sub_array (n-1,n,z_n,temp_---------Vervolstaändigen----------------
    
    
    
    temp_cvec = &z_n[n*n-n-1];     
    
    vector_dyadic_self2 ( temp_cvec,temp_V,n,n+1 );

    c = cos ( norms[n-1] );
    s = sin ( norms[n-1] );

    matrix_scalar_mult ( temp_V, - ( 1-c ),n+1 );

    for ( int i = 0; i < n; ++i )
        {

        temp_V[i* ( n+2 )]= _1+temp_V[i* ( n+2 )];

        }

    vector_scalar_mult ( temp_cvec,s,n );

    for ( int i = 0; i < n; ++i )

        {

        temp_V[n+i* ( n+1 )] = temp_cvec[i];
        temp_V[i+n* ( n+1 )] = -std::conj ( temp_cvec[i] );


        }

    temp_V[n * ( n+2 )]= c;

    }


void unitary_jarlskog::j_expand ( int klein, int gross )
    {

    matrix_initialize_zero ( temp_representation,gross );

    for ( int i = 0; i < klein; ++i )
        {
        for ( int j = 0; j < klein; ++j )
            {

            temp_representation[i+j*gross]=representation[i+j*klein];

            }

        }

    for ( int i = 0; i < gross*gross; ++i )
        {
        representation[i]=temp_representation[i];
        }

    }

void unitary_jarlskog::set_paramters ( double * p )
    {
      
      for (int i = 0 ; i < s ; ++i)
      {
	phases[i]=parameters[i];	
      }
      
      int a = 0;
      for (int i = s ; i < 2 * s - 2; ++i)
      {
	norms[a]=parameters[i];
	++a;
      }
      
      a = 0;
      for (int i = 2 * s - 2 ; i < s * s + s - 2; ++i)
      {
	theta_z_n[a] = parameters[i];
	++a;
      }
      

    }
    
void unitary_jarlskog::calculate_n_sphere_vectors()
    {
      
      z_n[0]= exp (std::complex< double > (0,theta_z_n[0]));
      
      
      for (int i = 2; i < s; ++i)
      {
	for (int k = 0; k <= 2 * i ; ++k)
	{
	  
	  c = cos(theta_z_n[ i * i - i - 1 + k ]);
	  s = sin(theta_z_n[ i * i - i - 1 + k ]);  // i*i-i-1 = index in theta_z_n wo die winkel zum i-ten vector stehen. 2*i ist die länge des Abschnitts.  
	  
	  for (int l = 0; l < i; ++ l)
	  {
	    
	  z_n[i*(i-1)+l];// <------------------------------------------------------
	  
	    
	  }
	  
	  
	}
	
	
	
      }

    }

