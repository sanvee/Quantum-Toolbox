#include "quit_unitary.h"

using namespace std;

void unitary_jarlskog::calculate_representation ( void ) //<----------------vervollständigen.
    {
     
    matrix_initialize_zero(representation,dim);
    matrix_initialize_zero(temp_V,dim);
    matrix_initialize_zero(temp_representation,dim);
    matrix_initialize_zero(temp_V_dag,dim);
    set_paramters();
    calculate_n_sphere_vectors();
   // std::cout << "z_n...";
   // vector_show(z_n,nr_zn);

    c = cos ( norms[0] );
    s = sin ( norms[0] );

    representation[0] = c;
    representation[1] = s * exp ( std::complex< double > ( 0,theta_z_n[0] ) );;
    representation[2] =-s * exp ( std::complex< double > ( 0,-theta_z_n[0] ) );;
    representation[3] = c;
    
   // matrix_show(representation,2);

    for ( int n = 2; n < dim; ++n )
        {

        j_expand ( n,n+1 ); //expands repr. adding one dim .. left and bottom zeros.
	//matrix_show(representation,n+1);

        representation [n * ( n+2 )] = _1; // (n+1)*(n+1)-1 letzes element der matrix
	//matrix_show(representation,n+1);
        create_V_n_n ( n );  // n = current size - 1
	//std::cout << "V-n_n" << std::endl;
	//matrix_show(temp_V,n+1);
        matrix_mlt ( representation,temp_V,temp_representation,n+1 );
        matrix_copy ( temp_representation,representation,n+1 );
	//matrix_show(representation,n+1);
	
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


    // std::cout << n* ( n-1 )/2 << " " << n << std::endl;
    std::complex <double> * p = &z_n[n* ( n-1 ) /2];
    
  //  std::cout << "vector";
 //   vector_show(p,n);

    vector_dyadic_self2 ( p,temp_V,n,n+1 );
//    cout << "dyadic" << endl;
//    matrix_show(temp_V,n+1);

    c = cos ( norms[n-1] );
    s = sin ( norms[n-1] );

    matrix_scalar_mult ( temp_V, - ( 1-c ),n+1 );
    
//    cout << "mal -(1-c)" << endl;
//    matrix_show(temp_V,n+1);
    
    for ( int i = 0; i < n; ++i )
        {

        temp_V[i* ( n+2 )]= _1+temp_V[i* ( n+2 )]; //<----diagonale immer i*n+1 du eierkopf :-D

        }
  //      cout << "minus _unitiy _" << endl;
//matrix_show(temp_V,n+1);
    vector_scalar_mult ( p,s,n );

    for ( int i = 0; i < n; ++i )

        {

        temp_V[n+i* ( n+1 )] = p[i];
        temp_V[i+n* ( n+1 )] = -std::conj ( p[i] );

        }
       // cout << "die vectoren kommen hinzu" << endl;
       // matrix_show(temp_V,n+1);
    temp_V[n * ( n+2 )]= c;
//matrix_show(temp_V,n+1);
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

void unitary_jarlskog::set_paramters ()
    {

    for ( int i = 0 ; i < dim ; ++i )
        {
        phases[i]=parameters[i];
        }

    int a = 0;
    for ( int i = dim ; i < dim + dim - 1; ++i )
        {
        norms[a]=parameters[i];
        ++a;
        }

    a = 0;
    for ( int i = dim + dim -1 ; i < dim * dim ; ++i )
        {
        theta_z_n[a] = parameters[i];
        ++a;
        }

    }

void unitary_jarlskog::calculate_n_sphere_vectors()
    {
    int i_theta_l,i_theta,i_coefs,i_coefs_l;

    for ( int i = 0 ; i <  dim * ( dim-1 )-2; ++i ) 
	{
            z_n_coefs[i]=1.0;
        }
  

    //std::cout << theta_z_n[0] << std::endl;
    //std::cout << z_n[0] << std::endl;

    for ( int i = 1; i < dim-1; ++i ) // dim = 3 -> erster index der winkel ist 1 bis 3 nächste 4 -> i*i.
        {
        i_theta = i*i;
        i_theta_l = 2*i+1;
        i_coefs = ( i+1 ) *i-2;
        i_coefs_l = 2* ( i+1 );

        //std::cout << i_theta << " " << i_theta_l << "   " << i_coefs << " " << i_coefs_l << std::endl;


        for ( int k = 0; k < i_theta_l; ++k )
            {

            c = cos ( theta_z_n [ i_theta + k ] );
            s = sin ( theta_z_n [ i_theta + k ] );

            z_n_coefs[i_coefs+k]*=c;

            //std::cout << i_theta + k <<std::endl;


            for ( int j = k+1; j < i_coefs_l; ++j )
                {

                z_n_coefs[i_coefs+j] *= s;
                // std::cout <<" "<< i_coefs+j;


                }
            //  std::cout << std::endl;
            }
        }
    if ( dim > 2 )
        {
        for ( int i = 0; i < ( ( dim - 1 ) * dim - 2 ) /2; ++i )
            {

           //  std::cout << " " << i+1 << "->" << i*2 << "," << i*2+1 << std::endl;
            z_n[i+1] = std::complex< double > ( z_n_coefs[2*i] , z_n_coefs[2*i+1] );

            }
          //  vector_show(z_n,2);
        }

    }
void unitary_jarlskog::initialize ( int s )
    {
       dim = s;
        nr_parameters = s * s ;
        nr_phases = s ;
        nr_norms = s-1 ;
        nr_theta = ( s-1 ) * ( s-1 ) ;
        nr_zn = s * ( s-1 ) /2 ;

        phases = new double [s]{};
        z_n = new std::complex< double >[s * ( s-1 ) /2]{};
        z_n_coefs = new double [  s * ( s-1 )-2]{};
        theta_z_n = new double[ ( s-1 ) * ( s-1 )]{};
        norms = new double [s-1]{};

        representation = new std::complex< double >[s*s]{};

        parameters = new double [s * s]{};

        temp_cvec = new std::complex< double >[s-1]{};
        temp_V = new std::complex< double >[s*s]{};
        temp_V_dag = new std::complex< double >[s*s]{};
        temp_representation = new std::complex< double >[s*s]{};

    }




