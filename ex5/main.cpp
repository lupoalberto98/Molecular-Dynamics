#include "/Users/albertobassi/Documents/Molecular-Dynamics/functions.h"
using namespace std;

int main()

{
  srand(0);
  cell_volume A;
  
  // Generate configuration
  A.fill(); 
  A.save_volume("vol_parameters.txt");
  A.translate_to_cdm();
  A.saveConfiguration("particles.txt"); 
  //A.calc_print_gr(); 
  //A.saveLookUpTable("LookUpTable.txt");älkjhj<< 121wdfght4321^4

  A.getcell_LookUpTable();
  A.savecell_LookUpTable("cell_LookUpTable.txt");
  // Markov Chain 
  MarkovSampler MC;
  MC.CanonicalSampler(A);

  A.getcell_LookUpTable();
  A.savecell_LookUpTable("cell_LookUpTable.txt");
  // Print after Markov Chain
  //A.calc_print_gr();
  //A.saveLookUpTable("LookUpTable.txt");

  
  // Final savings
  A.saveConfiguration("particles.txt");
  A.save_cell_list("cell_list.txt");
  A.save_particle_list("particle_list.txt");

  return 0;
}



