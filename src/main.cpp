#include "header.h"
using namespace std;

int main()

{
  srand(0);
  cell_volume A;
  A.__init__(1000, 10, 0.3, 1.0, 1.0);
  
  // Generate the configuration 
  A.fill(); 
  A.save_volume("vol_parameters.txt");

  // Initial calculations
  A.fill_lists();
  A.getcell_LookUpTable();
  A.get_kinetic_en();
  A.calculate_potential(1, 1);
  A.calculate_forces(1,1);
  cout<<"Initial potential "<<A.potential<<endl;
  cout<<"Initial kinetic energy "<<A.kinetic_en<<endl;

  // Initial savings
  A.save_cell_list("cell_list.txt");
  A.save_particle_list("particle_list.txt");
  A.savecell_LookUpTable("cell_LookUpTable.txt");
  A.saveConfiguration("particles.txt"); 
  
  
  // Uniform the system
  MarkovSampler MC;
  MC.__init__(5000, 0.1);
  MC.CanonicalSampler(A);

  
  
  // Molecular dynamics
  A.md_dynamics(10, 0.001, 1, 1);
  
  
  // Final calculations
  A.fill_lists();
  A.getcell_LookUpTable();
  A.get_kinetic_en();
  A.calculate_potential(1, 1);
  A.calculate_forces(1,1);
  cout<<"Final potential "<<A.potential<<endl;
  cout<<"Final kinetic energy "<<A.kinetic_en<<endl;

  // Final savings
  A.save_cell_list("cell_list.txt");
  A.save_particle_list("particle_list.txt");
  A.savecell_LookUpTable("cell_LookUpTable.txt");
  A.saveConfiguration("particles.txt"); 
  

  return 0;
}



