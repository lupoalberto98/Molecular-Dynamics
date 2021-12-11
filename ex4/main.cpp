#include "/Users/albertobassi/Documents/Computational_Physics/functions.h"
using namespace std;

int main()

{
  srand(0);
  volume A;
  MarkovSampler MC;
  A.fill();
  A.saveConfiguration("particles.txt"); // starting configuration
  A.getLookUpTable();
  A.saveLookUpTable("LookUpTable.txt");
  
  unsigned moves;
  double delta;
  cout<<"Insert number of moves"<<endl;
  cin>>moves;
  cout<<"Insert delta"<<endl;
  cin>>delta;
  
  MC.CanonicalSampler(moves, delta, A);
  A.getLookUpTable();
  A.saveConfiguration("particles.txt"); // Final configuration
  A.saveLookUpTable("LookUpTable.txt");
  
  return 0;
}



