#include "header.h"

using namespace std;



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

