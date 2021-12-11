#include "/Users/albertobassi/Documents/Computational_Physics/functions.h"
using namespace std;

int main()

{
  srand(0);
  lattice ising;
  
  ising.visualize();
  ising.setPbc();
  cout<<"Energy of configuration: "<<ising.energy()<<endl;

  return 0;
}



