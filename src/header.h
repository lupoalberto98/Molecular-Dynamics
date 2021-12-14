#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <random>
//#include <math.h>
#include <string>
#include <stdlib.h>
#include <time.h>
#include <algorithm>
#include <complex>

#ifndef FUNCTIONS_H
#define FUNCTIONS_H

using namespace std;

double LCG(int, int, int, int, int);
double mean(vector <double>);
double var(vector <double>);
double two_point_corr(vector <double>, int);


class lattice
{
    public:
    unsigned N_1, N_2;
    double J;
    bool Pbc;
    vector <int> *grid = new vector<int>;

    lattice();
    ~lattice();
   
    void visualize();
    void setPbc();
    void setJ();
    double energy();
};

class vec
{
    public:
    double x, y, z;

    vec();
    void __init__();
    double distanceFromOrigin();
    double distance(vec);
};

double dot(const vec&, const vec&);
vec getVector(const vec&, const vec&);
vec getVersor(const vec&, const vec&);

class particle : public vec
{
    public:
    double sigma; // Diameter
    double mass; 
    double kinetic_en;
    double vx, vy, vz, fx, fy, fz;
    particle();
    vec to_vector();
    void get_kinetic_en();

};


class volume 
{
    public:
    unsigned N; // Target of particles to generate
    double phi; // Packing factor
    double rho; // Density
    double L; // Box side length
    double T; // Temperature of the system
    double kinetic_en; // Kinetic energy

    ofstream debug;
    vec cdm; // Center of mass
    vec mom; // Total momentum of the system
    particle *configuration; // Sequence of particles
    double *LookUpTable; // Look up table for inter-particle distances
    

    volume();
    ~volume();
    void __init__(const unsigned&, const double&, const double&, const double&);

    void save_volume(string);
    void applyPbc(vec&);
    double distance(const vec&, const vec&);
    void fill();
    //void add_rdm_particle(double&, default_random_engine&, unsigned&);
    //void rdm_fill();
    void destroyParticle();
    void saveConfiguration(string);
    void readConfiguration(string);
    void get_phi();
    void get_cdm();
    void translate_to_cdm();
    void get_mom();
    void translate_to_mom_null();
    particle getParticle(const unsigned&);
    void swapParticles(const unsigned&, const unsigned&);
    void getLookUpTable();
    void saveLookUpTable(string);
    double getXminAverage();
    void calc_print_gr();
    void get_kinetic_en();
};

class cell_volume : public volume
{
    /**
     * @brief 
     * 
     */
    public:
    unsigned M_max;
    double rc;
    double potential;
    vector<unsigned> *cell_list;
    unsigned *particle_list;
    vector<double> *cell_LookUpTable; 

    cell_volume();
    ~cell_volume();
    void __init__(const unsigned&, const unsigned&, const double&, const double&, const double&);

    vec get_center_cell(const unsigned &);
    unsigned get_index_cell(vec&);
    void empty_cell_list();
    void fill_lists();
    void save_particle_list(string);
    void save_cell_list(string);
    unsigned* nearest_cells(const unsigned&);
    void getcell_LookUpTable();
    void savecell_LookUpTable(string);
    int determine_overlap(const unsigned&);
    double getXminAverage();
    void calculate_potential(const double&, const double&);
    void calculate_forces(const double&, const double&);
    void md_step(const double&, const double&, const double&);
    void mix_system(const double&, const double&, const double&, const unsigned&, const double&);
    void md_dynamics(const double&, const double&, const double&, const double&);
    complex<double> calculate_ssf(const vec&);

};

class MarkovSampler 
{
    public:
    unsigned moves;
    double delta;

    MarkovSampler();
    void __init__(const unsigned&, const double&);
    void CanonicalSampler(volume&);
    void CanonicalSampler(cell_volume&);
};

#endif // FUNCTIONS_H
