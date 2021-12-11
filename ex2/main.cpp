#include "/Users/albertobassi/Documents/Molecular-Dynamics/functions.h"
using namespace std;

int main()

{

  // x is the even sequence, y the odd one
  int N_seq, k;
  double  u1, u2, z1, z2;
  double sigma, mean_val;
  vector <double> z1_seq, z2_seq;

  default_random_engine generator;
  uniform_real_distribution<double> uniform(0.0,1.0);
  
  cout<<"Insert number of points"<<endl;
  cin>>N_seq;
  cout<<"Insert std"<<endl;
  cin>>sigma;
  cout<<"Insert mean"<<endl;
  cin>>mean_val;
  cout<<"Insert k for two point correlation function"<<endl;
  cin>> k;

  ofstream gaus("gaus.txt");
  
  for(int i=0; i<N_seq; ++i){
    //Generate the uniform distribution
    u1 = uniform(generator);
    u2 = uniform(generator);

    // Compute gaussian
    z1 = sqrt(-2*log(u1))*cos(2*M_PI*u2)*pow(sigma,2)+mean_val;
    z2 = sqrt(-2*log(u1))*sin(2*M_PI*u2)*pow(sigma,2)+mean_val;

    //cout<<z1<<z2<<endl;
    // Push back
    z1_seq.push_back(z1);
    z2_seq.push_back(z2);
    
    // Print on a file
    gaus<<z1<<"   "<<z2<<endl;
  
  }

  // Compute average and variance
  cout<<"Average(z1): "<<mean(z1_seq)<<"Variance: "<<var(z1_seq)<<endl;
  cout<<"Average(z2): "<<mean(z2_seq)<<"Variance: "<<var(z2_seq)<<endl;

  // Compute two point correlation function
  cout<<"Two point correlation function(z1):  "<<two_point_corr(z1_seq,k)<<endl;
  cout<<"Two point correlation function(z2):  "<<two_point_corr(z2_seq,k)<<endl;

}



