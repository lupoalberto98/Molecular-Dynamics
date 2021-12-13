#include "header.h"

using namespace std;


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

