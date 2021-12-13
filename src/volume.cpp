#include "header.h"

using namespace std;

double mean(vector <double> seq){
  int N = seq.size();
  double sum=0;
    for (int i=0;i<N;++i)
    {
        sum += seq[i];
    }
    return sum/N;
}

double var(vector <double> seq){
   int N = seq.size();
   double m = mean (seq);
   double sum = 0;
    for (int i=0;i<N;++i)
    {
        sum += pow(seq[i]-m,2);
    }
    return sqrt(sum/(N-1));
}


double two_point_corr(vector <double> seq, int k){

  int N = seq.size();
  double sum = 0.0;
  double mu = mean(seq);
  double variance = var(seq);
  
  for(int i=0; i<N-k; ++i){
    sum += (seq[i]-mu)*(seq[i+k]-mu);}

  return sum/(variance*(N-k));
}



// Class lattice 
lattice::lattice(){
  cout<<"Insert the dimension of the lattice [N_1, N_2]"<<endl;
  cin>>N_1;
  cin>>N_2;
  Pbc = true;
  J = 1.0;
  default_random_engine generator;
  uniform_int_distribution<int> int_uniform(0,1);

  for(unsigned j=0; j< N_1; ++j){
    for(unsigned i=0; i<N_2;++i){
      (*grid).push_back(2*int_uniform(generator)-1);
    }
  }
}

lattice::~lattice(){
  delete grid;
}

void lattice::visualize(){
  cout<<"Configuration:"<<endl;
  for(unsigned i=0; i<N_2; ++i){
    for(unsigned j=0; j<N_1; ++j){
      cout<<(*grid)[j*N_2+i];
    }
    cout<<endl;
  }
}

void lattice::setPbc(){
  cout<<"Periodic boundary conditions?"<<endl;
  cin>>Pbc;
}
void lattice::setJ(){
  cout<<"Constant J?"<<endl;
  cin>>J;
}

double lattice::energy(){
  double temp = 0.0;
  for(unsigned j=0; j< N_1-1; ++j){
    for(unsigned i=0; i<N_2-1;++i){
      temp += (*grid)[j*N_2+i]*((*grid)[(j+1)*N_2+i]+ (*grid)[j*N_2+i+1]);
    }
  }

  for(unsigned i=0; i<N_2-1; ++i){
    temp += (*grid)[(N_1-1)*N_2+i]*(*grid)[(N_1-1)*N_2+i+1];
    if(Pbc){
      temp += (*grid)[i]*(*grid)[(N_1-1)*N_2+i];
    }
  }
  for(unsigned j=0; j<N_1-1; ++j){
    temp += (*grid)[j*N_2+N_2-1]*(*grid)[(j+1)*N_2+N_2-1];
    if(Pbc){
      temp += (*grid)[j*N_2]*(*grid)[j*N_2+N_2-1];
    }
  }
  if(Pbc){
    temp += (*grid)[N_2-1]*(*grid)[N_1*N_2-1];
    temp += (*grid)[(N_1-1)*N_2]*(*grid)[N_1*N_2-1];
  }

return -J*temp;
}

// Class vec

vec::vec(){
  // Constructor
  x = 0.0;
  y = 0.0;
  z = 0.0;
}

void vec::__init__(){
  // Overloaded constructor
  x = 0.0;
  y = 0.0;
  z = 0.0;
}

double vec::distanceFromOrigin(){
  return sqrt(x*x + y*y + z*z);
}

double vec::distance(vec B){
  // Distace between two vectors in R^3
  double rx = x - B.x;
  double ry = y - B.y;
  double rz = z - B.z;
  return sqrt(rx*rx + ry*ry + rz*rz);
}

vec getVector(const vec&A, const vec& B){
  // Get the vector  pointing from B to A
  vec C;
  C.x = A.x - B.x;
  C.y = A.y - B.y;
  C.z = A.z - B.z;
  return C;
}

vec getVersor(const vec&A, const vec& B){
  // Get the versor pointing from B to A
  vec C = getVector(A,B);
  double norm = C.distanceFromOrigin();
  C.x = C.x/norm;
  C.y = C.y/norm;
  C.z = C.z/norm;
  return C;
}

// Class particle

particle::particle(){
  sigma = 1.0;
  mass = 1.0;
  kinetic_en = 0.0;
  vx = 0.0;
  vy = 0.0;
  vz = 0.0;
  fx = 0.0;
  fy = 0.0;
  fz = 0.0;
}

vec particle::to_vector(){
  vec temp;
  temp.x = x;
  temp.y = y;
  temp.z = z;
  return temp;
}

void particle::get_kinetic_en(){
  kinetic_en = mass*(vx*vx + vy*vy + vz*vz)*0.5;
}


// Class volume

volume::volume(){
  /* Inizialize the volume with zero particles
  */
  N = 0;
  phi = 0.0;
  T = 1.0;
  rho = 0.0;
  L = 10;
  kinetic_en = 0.0;
}

volume::~volume(){
  delete[] configuration;
  delete[] LookUpTable;
  debug.close();
}

void volume::__init__(const unsigned& N_target, const double& phi_target, const double& T_target, const double& sigma_target){
  /* Overload initializator  for volume
  Parameters:
  N_target = (unsigned) number of particles
  phi_target = (double) packing factor
  T_target = (double) temeprature
  sigma_target = (double) diameter of particles
  */

  // Set constant parameters
  N = N_target;
  phi = phi_target;
  rho = 6.0*phi/M_PI;
  L = sigma_target*cbrt((double) N/rho);
  T = T_target;
  kinetic_en = 3.0/2.0*T*N; // Set the kinetic energy according to the input temperature
  
  // Check if the packing factor is sufficiently small
  unsigned N_max = L/sigma_target;
  if(N > N_max*N_max*N_max)
  {debug<<"Error: packing factor too high. Program aborted."<<endl;
  abort();}


  // Open debug file
  debug.open("debug_log.txt");
  if(debug.is_open()==false){
    cout<<"Error: impossible to open debug file."<<endl;
  }

  // Allocate memory
  configuration = new particle[N];
  LookUpTable = new double[N*(N-1)];

}

void volume::save_volume(string filename){
  ofstream out(filename);
  out<<phi<<" "<<rho<<" "<<L<<" "<<kinetic_en<<endl;
  out.close();
}

void volume::applyPbc(vec& A){
  if(A.x > L/2) {A.x -= L;}
  if(A.x < -L/2) {A.x += L;}
  if(A.y > L/2) {A.y -= L;}
  if(A.y < -L/2) {A.y += L;}
  if(A.z > L/2) {A.z -= L;}
  if(A.z < -L/2) {A.z += L;}
}

double volume::distance(const vec &A, const vec &B){
  // Compute distance with periodic boundary conditions
  vec C = getVector(A, B);
  applyPbc(C);
  return C.distanceFromOrigin();
}

void volume::fill(){
  // Fill the volume compactly with N particles of random velocities such that the total momentum of the system is 0 and the temperature is T
  clock_t start, end; // Clocks to compute total time
  particle temp; // Temporary particle
  unsigned N_max = L/temp.sigma;
  // Random generator
  default_random_engine generator;
  uniform_real_distribution<double> uniform(-0.5,0.5);
  
 

  // Generate N non-overlapping particles with centers in the box [-L/2,L/2]^3
  start = clock();
  cout<<"Generating "<<N<<" non-overlapping particles..."<<endl;
  for(unsigned n=0; n<N; ++n ){
    // Compute the positions in a cubic centred lattice
    unsigned i = (n%N_max)%(N_max*N_max);
    unsigned j = (n/N_max)%N_max;
    unsigned k = (n/(N_max*N_max))%N_max;
    temp.x = -L/2 + temp.sigma/2+ i*temp.sigma;
    temp.y = -L/2 + temp.sigma/2 + j*temp.sigma;
    temp.z =  -L/2 + temp.sigma/2 + k*temp.sigma;
    

    // Generate random velocities 
    temp.vx = uniform(generator);
    temp.vy = uniform(generator);
    temp.vz = uniform(generator);
  
    // Save the particle in the configuration
    configuration[n] = temp;
  }

  // Translate the system to zero momentum and zero cdm
  translate_to_cdm();
  translate_to_mom_null();

  // Compute kinetic energy and intial temperature
  get_kinetic_en();
  double T0 = 2.0/3.0*kinetic_en/((double) N);

  // Rescale velocities 
  for(unsigned n=0; n<N; ++n){
    configuration[n].vx *= sqrt(T/T0);
    configuration[n].vy *= sqrt(T/T0);
    configuration[n].vz *= sqrt(T/T0);
  }

  // Compute right kinetic energy
  get_kinetic_en();

  end = clock();
  // Output the number of particles accepted
  cout<<"Total time to generate "<<N<<" particles: "<<1.0*(end-start)/CLOCKS_PER_SEC<<" s"<<endl;
  
}

/*
void volume::add_rdm_particle(double &sigma_target, default_random_engine& generator, unsigned &counter){
  // Add a particle of diameter sigma to the volume that non-overlap with the others
  uniform_real_distribution<double> uniform(-L/2,L/2);

  // Check if it does not overlap with the previously generated points
  unsigned N_start = N; // Number of particles in the volume
  while(getNumberParticles()<N_start+1){
    // Generate a random particle with diameter sigma, uniformly distributed in the box
    particle temp;
    ++counter;
    temp.x = uniform(generator);
    temp.y = uniform(generator);
    temp.z = uniform(generator);
    temp.sigma = sigma_target; 
   
    int nonoverlapping = 1; 
    for(unsigned j=0; j<N_start;++j){
      double actual_distance = distance(configuration[N_start-j-1],temp);
      double min_distance = (temp.sigma + configuration[N_start-j-1].sigma)*0.5;
      if(actual_distance > min_distance){nonoverlapping *= 1;}
      else{nonoverlapping *= 0;}
    }
    if(nonoverlapping==1){configuration.push_back(temp);}
  }
}


void volume:: rdm_fill(){
  // Fill the volume with randomly generating particles
  clock_t start, end;
  unsigned N_target;
  double phi_target;
  // Random generator
  default_random_engine generator;
  
  configuration.clear();// Empty the volume before
  cout<<"Insert the number of particles N"<<endl;
  cin>>N_target;
  cout<<"Insert desired packing factor"<<endl;
  cin>>phi_target;

  if(phi_target>0.38)
  {cout<<"Error: packing factor too high."<<endl;
  abort();}

  phi = phi_target; // Match the packing factor

  // Set the diameter ton match packing factor
  double sigma_target = 2*cbrt(phi_target*L*L*L*3/(N_target*4*M_PI));

  unsigned counter = 0;
  // Generate N non-overlapping particles with centers in the box [-L/2,L/2]^3
  start = clock();
  cout<<"Generating "<<N_target<<" non-overlapping particles..."<<endl;
  while(getNumberParticles()<N_target){
    add_rdm_particle(sigma_target, generator, counter);
  }
  end = clock();
  // Output the number of particles accepted
  cout<<"Total time to generate "<<N_target<<" particles: "<<1.0*(end-start)/CLOCKS_PER_SEC<<endl;
  cout<<"Accepted particles: "<<N_target<<" out of "<<counter<<endl;

}

*/

void volume::saveConfiguration(string filename){
  ofstream out(filename);
  for(unsigned n=0; n<N; ++n){
    out<<configuration[n].sigma<<" "<<configuration[n].mass<<" "<<configuration[n].x<<" "<<configuration[n].y<<" "<<configuration[n].z<<" ";
    out<<configuration[n].vx<<" "<<configuration[n].vy<<" "<<configuration[n].vz<<" ";
    out<<configuration[n].fx<<" "<<configuration[n].fy<<" "<<configuration[n].fz<<endl;
  }
  out.close();
}

void volume::readConfiguration(string filename){
  ifstream in(filename);
  double x;
  unsigned i = 0;
  particle temp;
  
  while(in>>x){
    if(i%4 == 0){temp.sigma = x;}
    else if(i%4 == 1){temp.x = x;}
    else if(i%4 == 2){temp.y = x;}
    else if(i%4 ==3) {
      temp.z = x;
      configuration[i/4] = temp;}
    ++i;
  }
  in.close();
}

void volume::get_phi(){
  double V_particles = 0.0;
  for(unsigned i=0; i<N; ++i){
    V_particles += 1.0*4/3*M_PI*pow(configuration[i].sigma/2, 3);
  }
  phi = V_particles/(L*L*L*1.0);  
}

void volume::get_cdm(){
  // Compute the center of mass of a given configuration
  cdm.__init__(); // Reinitialize cdm to origin
  if(N>0){
    for(unsigned i=0; i<N; ++i){
    cdm.x += configuration[i].x;
    cdm.y += configuration[i].y;
    cdm.z += configuration[i].z;
    }
    cdm.x /= (double) N;
    cdm.y /= (double) N;
    cdm.z /= (double) N;
  }
}

void volume::translate_to_cdm(){
  // Translate the system into the cdm
  get_cdm();
  for(unsigned i=0; i<N; ++i){
    configuration[i].x -= cdm.x;
    configuration[i].y -= cdm.y;
    configuration[i].z -= cdm.z;
  }
}

void volume::get_mom(){
  /* Compute the total momentum of the system
  */
  mom.__init__();
  if(N>0){
    for(unsigned i=0; i<N; ++i){
      mom.x += configuration[i].mass*configuration[i].vx;
      mom.y += configuration[i].mass*configuration[i].vy;
      mom.z += configuration[i].mass*configuration[i].vz;
    }
    mom.x /= (double) N;
    mom.y /= (double) N;
    mom.z /= (double) N;
  }
}

void volume::translate_to_mom_null(){
  // Translate systems velocities such that the total mometum is zero
  get_mom();
  for(unsigned i=0; i<N; ++i){
    configuration[i].vx -= mom.x/configuration[i].mass;
    configuration[i].vy -= mom.y/configuration[i].mass;
    configuration[i].vz -= mom.z/configuration[i].mass;
  }
}

particle volume::getParticle(const unsigned &i){
  return configuration[i];
}

void volume::swapParticles(const unsigned &i, const unsigned &j){
  particle *temp = new particle;
  *temp = getParticle(i);
  configuration[i] = getParticle(j);
  configuration[j] = *temp;
  delete temp;
}

void volume::getLookUpTable(){
  // Compute the look up table of all inter-particle distances as a vector
  // Save the distances between particle i and all the others sequentially, counting each distances twice in the end
 
  for(unsigned i=0; i<N; ++i){
    for(unsigned j=0; j<N; ++j){
      if(j<i){LookUpTable[i*(N-1)+j] = (distance(getParticle(i), getParticle(j)));}
      if(i<j){LookUpTable[i*(N-1)+j-1] = (distance(getParticle(i), getParticle(j)));}
    }
  }
}

void volume::saveLookUpTable(string filename){
  ofstream out(filename);
  for(unsigned i=0; i<N; ++i){
    for(unsigned j=0; j<N-1; ++j){
        out<<LookUpTable[i*(N-1)+j]<<" ";
    }
    out<<endl;
  }
  out.close();
}

double volume::getXminAverage(){
  double average = 0.0;
  for(unsigned i=0; i<N; ++i){
    average += *min_element(&(LookUpTable[i*(N-1)]), &(LookUpTable[i*(N-1)+N-2]));
  }
  return average/((double) N);
}

void volume::calc_print_gr(){
  //double rho = 6.0*phi/M_PI;
  unsigned N_hist = 200;
  double delta_r = L/2.0/((double) N_hist);
  string filename = "RDF.txt";

  // Compute all inter_particle distances
  getLookUpTable();

  // Save on the file
  ofstream out(filename);
  for(unsigned i=0; i<N_hist; ++i){
    double r_i = i*delta_r;
    double V_r = 4.0/3.0*M_PI*(pow(r_i+delta_r,3.0) - pow(r_i,3.0));
    unsigned particles_i = 0;
    for(unsigned j=0; j<N*(N-1); ++j){
      if(r_i <= LookUpTable[j] && LookUpTable[j] < r_i+delta_r){++particles_i;}
    }
    out<<r_i<<" "<<particles_i/((double) V_r*rho*N)<<endl;
  }

  out.close();
}

void volume::get_kinetic_en(){
  kinetic_en = 0.0;
  for(unsigned i=0; i<N; ++i){
    configuration[i].get_kinetic_en();
    kinetic_en += configuration[i].kinetic_en;
  }
}




