#include "header.h"

using namespace std;



// Class cell_volume

cell_volume::cell_volume(){
  /**
   * @brief Default inizializator of cell_volume
   * 
   */
  M_max = 1;
  rc = L;
  potential = 0.0;
}

cell_volume::~cell_volume(){
  /**
   * @brief Default destructor of cell_volume
   * 
   */
  delete[] cell_list;
  delete[] particle_list;
  delete[] cell_LookUpTable;
}

void cell_volume::__init__(const unsigned& N_target, const unsigned& M_max_target, const double& phi_target, const double& T_target, const double& sigma_target){
  /**
   * @brief Overloaded initializator of cell_volume class
   * 
   * Before initialization, check if the packing factor is too high, abort the program otherwise.
   * Check if M_max is too small. We want particles not to escape neighbouring cells
   * Finally allocate the memory.
   * 
   * @param N_target is the number of particles that we want to fill the volume with
   * @param M_max_target is the number of particles per coordinate direction that we want to divide the volume with
   * @param phi_target is the desired packing factor
   * @param T_target is the desired Temperature of initial configuration
   * @param sigma_target is the desired diameter of particles
   * 
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

  // Check if the packing factor is sufficiently small
  unsigned N_max = L/sigma_target;
  if(N > N_max*N_max*N_max)
  {cout<<"Error: packing factor too high. Program aborted."<<endl;
  abort();}

  // Check if every particle is sufficiently small to stay inside a cell
  if(rc < sigma_target*0.5){
    cout<<"Warning: rc is too small. Automatically set rc as the minimum value >= sigma/2."<<endl;
    M_max = (int) L*2.0*1/sigma_target;
    rc = L/((double) M_max);
  }
  
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
  /**
   * @brief Copmute the coordinates of the cell labeled with index m
   * 
   * @param m is the index of the cell
   * @return the vector center of the cell
   * 
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
  /**
   * @brief Compute the index of the cell in which a vector is located
   * 
   * @param A is the vector to consider
   * @return the index of its cell
   * 
   */

  applyPbc(A); // First trasnslate back the vector into the volume
  unsigned i =  (A.x+L/2)/rc;
  unsigned j = (A.y+L/2)/rc;
  unsigned k = (A.z+L/2)/rc;
  return i + j*M_max + k*M_max*M_max;
}

void cell_volume::empty_cell_list(){
  /**
   * @brief Empty cell_list
   * 
   */
  for(unsigned m=0; m<M_max*M_max*M_max; ++m){
    cell_list[m].clear();
  }
}

void cell_volume::fill_lists(){
  /**
   * @brief Compute cell_list and particle_list
   * 
   * Save in cell_list the indexes of the particles whose center lies inside.
   * Save in particle_list the index of the cell in which it lies
   */
  empty_cell_list();
  for(unsigned n=0; n<N; ++n){
    unsigned m = get_index_cell(configuration[n]);
    *(particle_list+n) = m;
    cell_list[m].push_back(n);
  }
}

void cell_volume::save_particle_list(string filename){
  /**
   * @brief Save particle_list in a file
   * 
   * @param filename is the name of the file
   * 
   */
  ofstream out(filename);
  for(unsigned n=0; n<N; ++n){
    out<<particle_list[n]<<endl;
  }
  out.close();
}

void cell_volume::save_cell_list(string filename){
  /**
   * @brief Save cell_list in a file
   * 
   * @param filename is the name of the file
   * 
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
  /**
   * @brief Compute neighbouring cells indexes
   * @param m is the index of a cell
   * @return the pointer to an array containing the indexes of all the neighbouring cells, with pbc, of m (included)
   * 
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
  /**
   * @brief Compute look up table
   * 
   * Determine inter particles distances between particles in neighbouring cells
   * Before fill_lists() must be called
   * 
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
  /**
   * @brief Save cell_LookUpTable in a file
   * 
   * @param filename is the name of the file
   * 
   */
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
  /**
   * @brief Determine if a particle overlaps with neighbouring ones
   * 
   * First fill_lists() must be called
   * 
   * @param i is the index of the particle
   * @return 1 if the particle does not overlap, 0 otherwise
   * 
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
  /**
   * @brief Compute minimum average distance
   * 
   * First getcell_LookUpTable() must be called
   * 
   * @return minium average distance between particles 
   */

 double xMinAverage = 0.0;
 for(unsigned n=0; n<N; ++n){
    xMinAverage += *min_element((*(cell_LookUpTable+n)).begin(), (*(cell_LookUpTable+n)).end());
 }
 return xMinAverage/((double) N);
}

void cell_volume::calculate_potential(const double& eps, const double& sig){
  /**
   * @brief Compute WCA potential of the system
   * 
   * getcell_LookUpTable() must be called before
   * 
   * @param eps is the energy scale
   * @param sig is the lenght scale 
   * 
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
  /**
   * @brief Compute total force acting on each particle
   * 
   * fill_lists() must be called before
   * 
   * @param eps is the energy scale of the system
   * @param sig is the length scale of the system 
   * 
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
  /**
   * @brief Update particle positions and velocties according to Verlet algorithm
   * 
   * Warning: till now this dynamics consider particles as points without dimension
   * Pay attention to call fill_lists() and getcell_LookUpTable() in the right order
   * 
   * @param dt is the temporal step
   * @param eps is the energy scale of the system
   * @param sig is the lenght scale of the system
   * 
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
    fill_lists();
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



void cell_volume::md_equilibrate(const unsigned& num_steps, const double& dt, const double& eps, const double& sig, const double& T0){
  /**
   * @brief Perform molecular dynamics for many steps
   * 
   * Equilibrate the system for a certain number of steps, by computing the average kinetic energy after a certain number of steps
   * and rescaling velocities to match the input temperature
   * 
   * @param num_steps is the number of steps during equilibration
   * @param dt is the temporal step
   * @param eps is the energy scale of the system
   * @param sig is the lenght scale of the system
   * @param T0 is the designed temperature for equilibration
   * 
   */

  unsigned step_av = num_steps/10;  
  unsigned save_step = num_steps/1000; 
  double ken_mean = 0;
  unsigned count = 0;
  clock_t start, end;

 
  start = clock();
  cout<<"Equilibrating the system..."<<endl;
  ofstream out("energy.txt");
  for(unsigned step=0; step<num_steps; ++step){
    md_step(dt, eps, sig); // Lists already updated in md_step
    if(step%save_step == 0){
      get_kinetic_en(); // Get the kinetic energy
      getcell_LookUpTable(); // Get look up tables
      calculate_potential(eps, sig); // Compute potential energy

      // Save
      out<<kinetic_en<<" "<<potential<<" "<<kinetic_en+potential<<endl;

      // Update average
      ken_mean += kinetic_en;
      ++ count;

      // Compute the averages every step_av steps
      if(step%step_av == 0 ){
        // Divide for the number of steps
        ken_mean /= (double) count;
      
        // Compute the temperature
        T = 2./(3.*N)*ken_mean;

        // Rescale velocities to match input temperature
        for(unsigned n=0; n<N; ++n){
          configuration[n].vx *= sqrt(T0/T);
          configuration[n].vy *= sqrt(T0/T);
          configuration[n].vy *= sqrt(T0/T);
        }

        // Reset averages to zero
        ken_mean = 0;
        count = 0;
      }
    }
    
  }
  end = clock();
  cout<<"Total time to perform equilibration: "<<1.0*(end-start)/CLOCKS_PER_SEC<<" s"<<endl;
  out.close();

}


complex<double> cell_volume::calculate_ssf(const vec& q, const double& dt, const double& eps, const double& sig){
  /** 
   * @brief Compute the static strcuture factor of WCA potential.
   *
   * Compute the term entering in average brackets that should be averaged during the mlecular dynamics.
   * See Exercises/CS08_2021.pdf for more details.
   * cell_LookUpTable() must be called before.
   * 
   * 
   * 
   * 
  */

  unsigned num_steps = 1000;
  unsigned step_av = num_steps/10;  // typically 1000
  unsigned save_step = num_steps/1000; // typically 10
  complex<double> ssf = 0.0;
  unsigned count = 0;
  const complex<double> I(0.0,1.0);   

  for(unsigned step=0; step<num_steps; ++step){
    md_step(dt, eps, sig);
    if(step%save_step == 0){
      // Loop over all possible pairs of particles in neighbouring cells
      for(unsigned i=0; i<N; ++i){
        for(unsigned j=0; j<cell_LookUpTable[i].size(); ++j){
          ssf += exp(I*dot(q, getVector(configuration[i], configuration[j])));
        }
      }

      ++count;
      if(step%step_av == 0){
        ssf /= count;
        ssf = ssf/((double) N)+1.;
      }
    }
  }
return ssf;
  
}

double cell_volume::calculate_heat_capacity(const unsigned& num_steps, const unsigned& num_save, const double& dt, const double& eps, const double& sig){
  /**
   * @brief Compute the heat capacity by kinetic energy fluctuations
   * Warning: to have meaningful results equilibrate the system before
   * 
   * @param num_steps is the total number of steps
   * @param num_save is how many steps to save the simulation
   * @param dt is the temporal step
   * @param eps is the energy scale of the system
   * @param sig is the lenght scale of the system
   * 
   */

  double ken_mean = 0;
  double ken_square_mean = 0;
  unsigned count = 0;
  clock_t start, end;

 
  start = clock();
  cout<<"Computing heat capacity..."<<endl;
  
  for(unsigned step=0; step<num_steps; ++step){
    md_step(dt, eps, sig); // Lists already updated in md_step
    if(step%num_save == 0){
      get_kinetic_en(); // Get the kinetic energy
      
      // Update average
      ken_mean += kinetic_en;
      ken_square_mean += kinetic_en*kinetic_en;
      ++ count;
    }
  }
  // Compute averages
  ken_mean /= (double) count;
  ken_square_mean /= (double) count;
  double T_mean = 2.*ken_mean/(3.*N);
  double ken_var = ken_square_mean - ken_mean*ken_mean;
  double heat_capacity = 3./2.*N/(1.-2.*ken_var/(3.*N*T_mean*T_mean));
  end = clock();
  cout<<"Total time to compute heat capacity: "<<1.0*(end-start)/CLOCKS_PER_SEC<<" s"<<endl;
  return heat_capacity;
}


double cell_volume::calculate_pressure(const unsigned& num_steps, const unsigned& num_save, const double& dt, const double& eps, const double& sig){
  /**
   * @brief Compute the pressure of the gas by virial expansion
   * 
   * @param num_steps is the total number of steps
   * @param num_save is how many steps to save the simulation
   * @param dt is the temporal step
   * @param eps is the energy scale of the system
   * @param sig is the lenght scale of the system
   * 
   */

  double ken_mean = 0;
  unsigned count = 0;
  double virial = 0;
  clock_t start, end;

 
  start = clock();
  cout<<"Computing pressure..."<<endl;
  for(unsigned step=0; step<num_steps; ++step){
    md_step(dt, eps, sig); // Lists already updated in md_step
    if(step%num_save == 0){
      get_kinetic_en();
      // Update average
      ken_mean += kinetic_en;
      ++ count;
      // Compute virial term
      for(unsigned i=0; i<N-1; ++i){
        for(unsigned j=i+1; j<N; ++j){
          // Copmute distance
          double r = distance(configuration[i], configuration[j]);
          // Compute magnitude of the radial repulsive force
          double fr = 24.0*eps*(2*pow(sig/r,12.0) - pow(sig/r,6.0))/(r*r);
          virial += fr*r;
        }
      }
    }
  }
  ken_mean /= (double) count;
  double T_mean = 2./(3.*N)*ken_mean;
  double pressure = N*T_mean + 1./3.*virial;
  pressure /= L*L*L;
  end = clock();
  cout<<"Total time to compute pressure: "<<1.0*(end-start)/CLOCKS_PER_SEC<<" s"<<endl;

  return pressure;
}

vec cell_volume::calculate_msd_vautocorr(const unsigned& num_steps, const unsigned& num_save, const double& dt, const double& eps, const double& sig){
  /**
   * @brief Compute the mean square displacement, save in a file and compute diffusion copefficient deriving last values
   * and copmute velocity auto-correlation function
   * 
   * @param num_steps is the total number of steps
   * @param num_save is how many steps to save the simulation
   * @param dt is the temporal step
   * @param eps is the energy scale of the system
   * @param sig is the lenght scale of the system
   * 
   */
  vec msd, msd_semilast, diff_coeff, v_autocorr;
  clock_t start, end;
  // Take initial configuration
  vector <particle> initial_configuration(configuration, configuration+N);
  unsigned integer_part = num_steps/num_save;

  start = clock();
  ofstream out("msd_vcorr.txt");
  cout<<"Computing mean square displacement and velocity auto-correlation function..."<<endl;
  for(unsigned step=0; step<num_steps; ++step){
    md_step(dt, eps, sig); // Lists already updated in md_step
    if(step%num_save == 0){
      msd.x = msd.y = msd.z = 0;
      v_autocorr.x = v_autocorr.y = v_autocorr.z = 0;
      for(unsigned n=0; n<N; ++n){
        // Mean square displacement
        msd.x += (configuration[n].x - initial_configuration[n].x)*(configuration[n].x - initial_configuration[n].x);
        msd.y += (configuration[n].y - initial_configuration[n].y)*(configuration[n].y - initial_configuration[n].y);
        msd.z += (configuration[n].z - initial_configuration[n].z)*(configuration[n].z - initial_configuration[n].z);
        // Velocity auto-correlation function
        v_autocorr.x += configuration[n].vx*initial_configuration[n].vx;
        v_autocorr.y += configuration[n].vy*initial_configuration[n].vy;
        v_autocorr.z *= configuration[n].vz*initial_configuration[n].vz;
      }
      // Divide to take the averages
      msd.x /= (double) N;
      msd.y /= (double) N;
      msd.z /= (double) N;
      v_autocorr.x /= (double) N;
      v_autocorr.y /= (double) N;
      v_autocorr.z /= (double) N;

      // Print on a file
      out<<msd.x<<" "<<msd.y<<" "<<msd.z<<" "<<v_autocorr.x<<" "<<v_autocorr.y<<" "<<v_autocorr.z<<endl;

      // Take second and third to last msd to compute diffusion coefficient
      if(step==num_steps-3*integer_part){msd_semilast = msd;}
    }
  }
  // Compute derivative and save in msd (output diffusion coefficient)
  diff_coeff.x = (msd.x-msd_semilast.x)/(4*dt);
  diff_coeff.y = (msd.y-msd_semilast.y)/(4*dt);
  diff_coeff.z = (msd.z-msd_semilast.z)/(4*dt);

  out.close();
  end = clock();
  cout<<"Total time to compute msd and velocity auto-correlation function: "<<1.0*(end-start)/CLOCKS_PER_SEC<<" s"<<endl;
  return diff_coeff;
}

