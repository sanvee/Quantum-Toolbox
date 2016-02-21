  int main (){
   gsl_rng * R_N_G = initialise_random_number_gernerators();
  
  std::ofstream TEXTFILE,
  NUM_DEBUG;
  TEXTFILE.open("../test_tab.txt");
  NUM_DEBUG.open("../suspekt.txt",std::ios_base::app);
  
  

  
  
  std::complex< double > _a = (1.0,0,0);
  std::complex< double > _n = (0.0,0.0); 
  
  std::complex< double >  eta [8*8]   = 
  {10,_n,_n,_n,_n,_n,_n,_n,
    _n,22,_n,_n,_n,_n,_n,_n,
    _n,_n,33,_n,_n,_n,_n,_n,
    _n,_n,_n,44,_n,_n,_n,_n,
    _n,_n,_n,_n,55,_n,_n,_n,
    _n,_n,_n,_n,_n,66,_n,_n,
    _n,_n,_n,_n,_n,_n,77,_n,
    _n,_n,_n,_n,_n,_n,_n,88};
    
    std::complex< double >  eta_prime [16*16];
    
    //show_matrix(eta,8);
    
    int dim_qubit = 2;
    
    std::complex< double > ONE_8x8[8*8];
    matrix_initialize_unity(ONE_8x8,8);
    
    
    
    
    std::complex< double > sigma_A [4*4],
    sigma_B [6*6],
    sigma_C [6*6]={1,1,1,1,1,1,
      2,2,2,2,2,2,
      3,3,3,3,3,3,
      4,4,4,4,4,4,
      5,5,5,5,5,5,
      6,6,6,6,6,6,},
      
      temp4x4 [4*4],
      sigma   [9*9],
      eigen1  [9],
      eigen2  [8],
      temp1   [8*8],
      temp2   [8*8],
      total_op[8*8],
      rn_untr [8*8],
      basis_v [8*8];
      
      
      
      
      
      
      // generate separable state sigma_A x sigma_B x sigma_C
      
      
      
      
      int subspaces[2] = {3,-2};
      int sep_struct [2] = {3,2};
      
      int tr1 [3] = {-2,2,2};
      int tr2 [3] = {2,-2,2};
      int tr3 [3] = {2,2,-2};
      
      double rel_ent = 0;
      double rel_ent_min = 10;
      int nr_errors = 0;
      std::cout.precision(10);
      /*  
       *  for (int i = 0 ; i < 100000; ++i){
       * 
       *    //state_decohere_on_subspace(eta,eta_prime,8,subspace,3,temp1,temp2,total_op,rn_untr,basis_v); 
       *    //rel_ent = state_rel_entropy(eta,eta_prime,eigen1,eigen2,8);
       *    
       *    //sample_separable_density_matrix_hs(sigma,8,sep_struct,3);
       *    //rel_ent = state_rel_entropy(eta,sigma,eigen1,eigen2,8);
       *    
       *    //sample_separable_density_matrix_up(sigma,6,sep_struct,2);
       *    sample_density_matrix_up(sigma,2);
       *    //sample_pure_density_matrix(sigma,6);
       *    //sample_separable_pure_density_matrix(sigma,6,sep_struct,2);
       *    
       *    //state_partial_transpose_nxn(sigma,sigma_B,6,subspaces,2);
       *    //state_partial_trace(sigma,6,sigma_B,3,subspaces,2);
       * 
       *    matrix_eigenvalues(sigma,eigen1,2);
       *    //vector_cutoff_neg(eigen1,3);
       *    //sort_vec_asc(eigen1,3);
       *    
       *    //TEXTFILE << state_negativity(sigma_B,6) << std::endl;
       *    
       *    
       *    for (int j = 0; j < 2; ++j){
       *    TEXTFILE << eigen1[j].real() << "\t";
}
TEXTFILE << std::endl;


//show_vector(eigen1,2);

/*if (rel_ent <  rel_ent_min && rel_ent  > 0.0){
 *      
 *    rel_ent_min = rel_ent;
 *    std::cout << "i "<< i<< " " << rel_ent << std::endl;
 *      
}

if (rel_ent < 1.0/3.0){
  ++nr_errors;
  std::cout << "i=" << i << " errno=" << nr_errors << std::endl; 
  std::cout << rel_ent << std::endl; 
  show_vector(eigen2,8);
  std::cout << "ERR" << std::endl; 
  
  for (int j = 0; j < 8; ++j){
    
    std::cout.precision(20);
    
    NUM_DEBUG << eigen2[j].real() << '\t';
}
NUM_DEBUG << std::endl;
}

//if (nr_errors > 100000)
//return 0;



if (rel_ent > 0){
  for (int j = 0; j < 1; ++j){
    std::cout.precision(20);  
    TEXTFILE << rel_ent << '\t';
}
TEXTFILE << std::endl;

TEXTFILE.close();
NUM_DEBUG.close();
}
}
*/