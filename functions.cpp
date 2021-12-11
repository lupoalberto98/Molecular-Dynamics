#include "functions.h"

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
  // Random generator
  default_random_engine generator;
  uniform_real_distribution<double> uniform(-0.5,0.5);
  
  
  unsigned N_max = L/temp.sigma;
  
  if(N > N_max*N_max*N_max)
  {debug<<"Error: packing factor too high. Program aborted."<<endl;
  abort();}

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





// Class cell_volume

cell_volume::cell_volume(){
  /* Initialize the cell volume
  Parameters:
  M_max = (unsigned) maximum number of cells per direction
  rc = cell length
  */
  M_max = 1;
  rc = L;
  potential = 0.0;
}

cell_volume::~cell_volume(){
  delete[] cell_list;
  delete[] particle_list;
  delete[] cell_LookUpTable;
}

void cell_volume::__init__(const unsigned& N_target, const unsigned& M_max_target, const double& phi_target, const double& T_target, const double& sigma_target){
  /* Overloaded initializator for cell volume
  Parameters:
   Parameters:
  N_target = (unsigned) number of particles
  M_max_target = (unsigned) maximum number of cells per direction
  phi_target = (double) packing factor
  T_target = (double) temeprature
  sigma_target = (double) diameter of particles
  */


  // Set constant parameters
  N = N_target;
  M_max = M_max_target;
  phi = phi_target;
  rho = 6.0*phi/M_PI;
  L = sigma_target*cbrt((double) N/rho);
  rc = L/((double) M_max);
  T = T_target;
  kinetic_en = 3.0/2.0*T*N; // Set the kinetic energy according to the input temperature
  
  

  // Open debug file
  debug.open("debug_log.txt");
  if(debug.is_open()==false){
    cout<<"Opening error (Backtrace: __init__() -> open debug file): impossible to open debug file."<<endl;
  }
  

  // Allocate memory
  configuration = new particle[N];
  LookUpTable = new double[N*(N-1)];
  cell_list = new vector<unsigned>[M_max*M_max*M_max];
  particle_list = new unsigned[N];
  cell_LookUpTable = new vector<double>[N];
}

vec cell_volume::get_center_cell(const unsigned &m){
  /* Return the coordinates of the celle of index m
  Parameters:
  m = (unsigned) index of the cell
  Returns:
  center = (vec) coordinates of the center
  */
  vec center;
  unsigned i = m%M_max;
  unsigned j = (m/M_max)%M_max;
  unsigned k = (m/(M_max*M_max))%M_max;
  center.x = -L/2 + rc/2+ i*rc;
  center.y = -L/2 + rc/2 + j*rc;
  center.z =  -L/2 + rc/2 + k*rc;
  return center;
}

unsigned cell_volume::get_index_cell(vec &A){
  /* Give the index of the  cell in which a vector is located
  Parameters:
  A = (vec) point in the volume
  Returns:
  i = (unsigned) index of the cell in which A is located
  */
  applyPbc(A); // First trasnslate back the vector into the volume
  unsigned i =  (A.x+L/2)/rc;
  unsigned j = (A.y+L/2)/rc;
  unsigned k = (A.z+L/2)/rc;
  return i + j*M_max + k*M_max*M_max;
}

void cell_volume::empty_cell_list(){
  /* Empty the cell list
  */
  for(unsigned m=0; m<M_max*M_max*M_max; ++m){
    cell_list[m].clear();
  }
}

void cell_volume::fill_lists(){
  /* Save in particle list the cell indexes of the particles
  Save in cell list the indexes of particles in the cell
   */
  empty_cell_list();
  for(unsigned n=0; n<N; ++n){
    unsigned m = get_index_cell(configuration[n]);
    *(particle_list+n) = m;
    cell_list[m].push_back(n);
  }
}

void cell_volume::save_particle_list(string filename){
  /* Save particle_list in a file
  Parameters:
  filename = (string) name of the file
  */
  ofstream out(filename);
  for(unsigned n=0; n<N; ++n){
    out<<particle_list[n]<<endl;
  }
  out.close();
}

void cell_volume::save_cell_list(string filename){
  /* Save cell_list in a file
  Parameters:
  filename = (string) name of the file
  */
  ofstream out(filename);
  for(unsigned m=0; m<M_max*M_max*M_max; ++m){
    for(unsigned j=0; j<cell_list[m].size(); ++j){
      out<<cell_list[m][j]<<" ";
    }
    out<<endl;
  }
  out.close();
}

unsigned* cell_volume::nearest_cells(const unsigned &m){
  /* Return the pointer to an array containing the indexes of 
  all the neighbouring cells, with pbc, of m (included)
  Parameters:
  m = (unsigned) index of the cell
  Returns:
  indexes = (unsigned*) poiinter to an array containing the indexes of neighbouring cells
  */
  unsigned *indexes = new unsigned[27];
  vec *center = new vec;
  for(int i=-1; i<2; ++i){
    for(int j=-1; j<2; ++j){
      for(int k=-1; k<2; ++k){
        *center = get_center_cell(m);
        (*center).x += (double) i*rc;
        (*center).y += (double) j*rc;
        (*center).z += (double) k*rc;
        applyPbc(*center); 
        unsigned ind = get_index_cell(*center);
        *(indexes+i*9 + j*3 +k + 13) = ind;
      }
    }
  }
  delete center;
  return indexes;  
}

void cell_volume::getcell_LookUpTable(){
  /* Determine if the particle of index i has an overlap with the particles in the neighbouring cells
  First fill_lists() must be called
  */
  for(unsigned i=0; i<N; ++i){
    cell_LookUpTable[i].clear();
    unsigned m = get_index_cell(configuration[i]);
    unsigned *indexes = new unsigned[27];
    indexes = nearest_cells(m);
    for(unsigned l=0; l<27; ++l){
      unsigned ind = *(indexes+l);
      for(unsigned n=0; n<cell_list[ind].size(); ++n){
          unsigned j = cell_list[ind][n];
          if(i != j){
            double dr = distance(configuration[i], configuration[j]);
            if(dr<(configuration[i].sigma+configuration[j].sigma)*0.49999999999999){
              debug<<"Overlap error (Backtrace: getcell_LookUpTable() -> distance dr): particle "<<i+1<<" overlaps with particle "<<j+1<<endl;
              debug<<"Distance is "<<dr<<endl;
            }
            cell_LookUpTable[i].push_back(dr);
          }
      }
    }
  }
}

void cell_volume::savecell_LookUpTable(string filename){
  ofstream out(filename);
  for(unsigned n=0; n<N; ++n){
    for(unsigned j=0; j<cell_LookUpTable[n].size(); ++j){
      out<<cell_LookUpTable[n][j]<<" ";
    }
    out<<endl;
  }
  out.close();
}

int cell_volume::determine_overlap(const unsigned&i){
  /* Determine if the particle of index i has an overlap with the particles in the neighbouring cells
  First fill_lists() must be called
  Parameters:
  i = (const unsigned&) index of the particle, passed by reference
  Returns:
  nonoverlapping = (int) if 1 there is an overlap with a particle in theneighbouring cells at lest, otherwise 0.
  */
  unsigned m = get_index_cell(configuration[i]);
  unsigned *indexes = new unsigned[27];
  indexes = nearest_cells(m);
  int nonoverlapping = 1;
  for(unsigned l=0; l<27; ++l){
    unsigned ind = *(indexes+l);
    for(unsigned n=0; n<cell_list[ind].size(); ++n){
        unsigned j = cell_list[ind][n];
        if(i != j){
          double dr = distance(configuration[i], configuration[j]);
          if(dr >=(configuration[i].sigma + configuration[j].sigma)/2){nonoverlapping *= 1;}
          else {nonoverlapping *= 0;}
        }
    }
    
  }
  return nonoverlapping;
}

double cell_volume::getXminAverage(){
  /* Compute the minimum average distance between particles by looking at neighbouring cells in the cell_LookUpTable
  First getcell_LookUpTable() must be called
  Returns:
  xMinAverage = (double) the average minuminum distance
  */
 double xMinAverage = 0.0;
 for(unsigned n=0; n<N; ++n){
    xMinAverage += *min_element((*(cell_LookUpTable+n)).begin(), (*(cell_LookUpTable+n)).end());
 }
 return xMinAverage/((double) N);
}

void cell_volume::calculate_potential(const double& eps, const double& sig){
  /* Computes the WCA potential of the system
  First getcell_LookUpTable() must be called
  Parameters:
  eps = (double) energy scale of the system
  sig = (double) length scale of the system
  Returns:
  void -> update potential
  */
  double r_cut = pow(2.0,1.0/6.0)*sig;
  potential = 0.0;
  for(unsigned n=0; n<N; ++n){
    for(unsigned j=0; j<cell_LookUpTable[n].size(); ++j){
      double r = cell_LookUpTable[n][j];
      if(r<r_cut){
        potential += 4*eps*(pow(sig/r,12.0) - pow(sig/r,6.0)) + eps;
      }
    }
  }
  potential *= 0.5;
}

void cell_volume::calculate_forces(const double& eps, const double& sig){
  /* Compute the forces acting on each particles
  First fill_lists() must be called
  Parameters:
  eps = (double) energy scale of the system
  sig = (double) length scale of the system
  */
  double r_cut = pow(2.0,1.0/6.0)*sig;
  for(unsigned n=0; n<N; ++n){
    // Set forces to zero
    configuration[n].fx = 0;
    configuration[n].fy = 0;
    configuration[n].fz = 0;
    // Get cell indexe and its neighbouring indexes
    unsigned m = get_index_cell(configuration[n]);
    unsigned *indexes = new unsigned[27];
    indexes = nearest_cells(m);
    for(unsigned l=0; l<27; ++l){
      unsigned ind = *(indexes+l);
      for(unsigned i=0; i<cell_list[ind].size(); ++i){
          // Retireve index of the particle
          unsigned j = cell_list[ind][i];
          if(n != j){
            double r = distance(configuration[n], configuration[j]);
            // Compute magnitude of the radial repulsive force
            double fr = 24.0*eps*(2*pow(sig/r,12.0) - pow(sig/r,6.0))/(r*r);
            // Check if fr is too high
            if(isnan(fr) == true){
              debug<<"Overflow error (Backtrace: calculate_forces() -> total radial force): the force that particle "<<j+1<<" excerts on particle "<<n+1<<" is too high."<<endl;
            }
            // Versor connecting j to n
            vec r_vec = getVector(configuration[n], configuration[j]);
            applyPbc(r_vec);
            // Update forces if the distance is less the the cutoff distance
            // Otherwise particles do not interact
            if(r<r_cut){
              configuration[n].fx += r_vec.x*fr;
              configuration[n].fy += r_vec.y*fr;
              configuration[n].fz += r_vec.z*fr;
            } 
          }
      }
    }
  }
}

void cell_volume::md_step(const double& dt, const double& eps, const double& sig){
  /* Update particle positions and velocties according to Verlet algorithm
  Pay attention to call fill_lists() and getcell_LookUpTable() in the right places
  Parameters:
  dt = (double) time step
  eps = (double) energy scale of the system
  sig = (double) length scale of the system
  */
  fill_lists();
  calculate_forces(eps, sig);
  for(unsigned n=0; n<N; ++n){
    // Intermidiate step for velocities
    configuration[n].vx += dt/(2.0*configuration[n].mass)*configuration[n].fx;
    configuration[n].vy += dt/(2.0*configuration[n].mass)*configuration[n].fy;
    configuration[n].vz += dt/(2.0*configuration[n].mass)*configuration[n].fz;
    // Step for positsions
    configuration[n].x += dt*configuration[n].vx;
    configuration[n].y += dt*configuration[n].vy;
    configuration[n].z += dt*configuration[n].vz;
    applyPbc(configuration[n]);
    int nonoverlapping = determine_overlap(n);
    if(nonoverlapping==0){
      debug<<"Overlap error (Backtrace: MD -> overlap after intermidiate step): particle "<<n+1<<" overlaps with some particle"<<endl;
    }
  }

  
  // Update forces
  fill_lists();
  calculate_forces(eps,sig);
  for(unsigned n=0; n<N; ++n){
    // Final step for velocities
    configuration[n].vx += dt/(2.0*configuration[n].mass)*configuration[n].fx;
    configuration[n].vy += dt/(2.0*configuration[n].mass)*configuration[n].fy;
    configuration[n].vz += dt/(2.0*configuration[n].mass)*configuration[n].fz;
  }
}

void cell_volume::md_dynamics(const double& total_time, const double& dt, const double& eps, const double& sig){
  /* Update particle positions and velocties according to Verlet algorithm
  Parameters:
  total_time = (double) total simulation time
  dt = (double) time step
  eps = (double) energy scale of the system
  sig = (double) length scale of the system
  */
  unsigned num_steps = total_time/dt;
  clock_t start, end;
  start = clock();
  cout<<"Performing molecular dynamics..."<<endl;
  ofstream out("energy.txt");
  for(unsigned step=0; step<num_steps; ++step){
    md_step(dt, eps, sig);
    getcell_LookUpTable();
    get_kinetic_en();
    calculate_potential(eps, sig);
    if(step%10 ==0){
      out<<step<<" "<<kinetic_en<<" "<<potential<<" "<<kinetic_en+potential<<endl;
    }
  }
  end = clock();
  cout<<"Total time to perform Molecular Dynamics: "<<1.0*(end-start)/CLOCKS_PER_SEC<<" s"<<endl;
  out.close();

}


// Class MarkovSampler

MarkovSampler::MarkovSampler(){
  moves = 1000;
  delta = 0.1;
}

void MarkovSampler::__init__(const unsigned& moves_target, const double& delta_target){
  moves = moves_target;
  delta = delta_target;
}

void MarkovSampler::CanonicalSampler(volume &A){
  particle *temp = new particle;
  clock_t start, end;
  string filename = "MonteCarloLog.txt";
  
  
  // Random generator
  default_random_engine generator;
  uniform_real_distribution<double> double_uniform(-1.0,1.0);
  uniform_int_distribution<int> int_uniform(0, A.N-1);

  start = clock();
  cout<<"Performing Monte Carlo simulation..."<<endl;
  ofstream out(filename);
  for(unsigned move=0; move<moves; ++move){
    unsigned accepted_attempts = 0;
    for(unsigned attempt=0; attempt<A.N; ++attempt){
      // Generate a direction
      double dx = double_uniform(generator);
      double dy = double_uniform(generator);
      double dz = double_uniform(generator);
      // Rescale to have norm = delta
      double norm = sqrt(dx*dx+dy*dy+dz*dz);
      dx = dx/norm*delta;
      dy = dy/norm*delta;
      dz = dz/norm*delta;
      // Uniformly choose a site i
      unsigned i = int_uniform(generator);
      int nonoverlapping = 1;
      
      // Move the particle in that site for a dstance delta
      *temp = A.configuration[i];
      (*temp).x += dx;
      (*temp).y += dy;
      (*temp).z += dz;

      // Apply periodic boundary conditions
      A.applyPbc(*temp);

      // Run over all the other particles
      for(unsigned j=0;j<A.N; ++j){
        if(j!=i){
          double actual_distance = A.distance(*temp, A.configuration[j]);
          if(actual_distance>((*temp).sigma + (A.configuration[j]).sigma)/2){nonoverlapping *= 1;}
          else{ nonoverlapping *= 0;}
        }
      }
      
      // Checkoverlapping
      if(nonoverlapping==1){
        A.configuration[i] = *temp;
        ++accepted_attempts; 
      }
    }
    
    
    if(move%100 == 0){
      A.saveConfiguration("particles.txt"); 
      A.getLookUpTable();
      double xmin_average = A.getXminAverage();
      double acceptance_rate = accepted_attempts*1.0/A.N;
      out<<move<<" "<<acceptance_rate<<" "<<xmin_average<<endl;
    }
  }
  end = clock();
  
  cout<<"Total time to perform Monte Carlo simulation: "<<1.0*(end-start)/CLOCKS_PER_SEC<<" s"<<endl;
  delete temp;
  out.close();
}

void MarkovSampler::CanonicalSampler(cell_volume &A){
  particle *proposal = new particle;
  particle *temp = new particle;
  clock_t start, end;
  string filename = "MonteCarloLog.txt";
  
  
  // Random generator
  default_random_engine generator;
  uniform_real_distribution<double> double_uniform(-1.0,1.0);
  uniform_int_distribution<int> int_uniform(0, A.N-1);

  start = clock();
  cout<<"Performing Monte Carlo simulation..."<<endl;
  ofstream out(filename);
  for(unsigned move=0; move<moves; ++move){
    unsigned accepted_attempts = 0;
    for(unsigned attempt=0; attempt<A.N; ++attempt){
      // Generate a direction
      double dx = double_uniform(generator);
      double dy = double_uniform(generator);
      double dz = double_uniform(generator);
      // Rescale to have norm = delta
      double norm = sqrt(dx*dx+dy*dy+dz*dz);
      dx = dx/norm*delta;
      dy = dy/norm*delta;
      dz = dz/norm*delta;
      // Uniformly choose a site i
      unsigned i = int_uniform(generator);

      // Move the particle in that site for a dstance delta
      *proposal = A.configuration[i];
      (*proposal).x += dx;
      (*proposal).y += dy;
      (*proposal).z += dz;

      // Apply periodic boundary conditions
      A.applyPbc(*proposal);

      // Swap proposal with configuration[i]
      *temp = A.configuration[i];
      A.configuration[i] = *proposal;
      *proposal = *temp;

      // Determine overlapping with neighbours
      A.fill_lists();
      int nonoverlapping = A.determine_overlap(i);
      
      // Checkoverlapping
      if(nonoverlapping==0){
        A.configuration[i] = *proposal;
      }
      else{++accepted_attempts;}
    }
    
    if(move%100 == 0){
      A.getcell_LookUpTable();
      double xMinAverage = A.getXminAverage();
      double acceptance_rate = accepted_attempts*1.0/A.N;
      out<<move<<" "<<acceptance_rate<<" "<<xMinAverage<<endl;
      //out<<move<<" "<<acceptance_rate<<" "<<xmin_average<<endl;
    }
  }
  end = clock();
  
  cout<<"Total time to perform Monte Carlo simulation: "<<1.0*(end-start)/CLOCKS_PER_SEC<<" s"<<endl;
  delete temp;
  delete proposal;
  out.close();
}

