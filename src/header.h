#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <random>
#include <math.h>
#include <string>
#include <stdlib.h>
#include <time.h>
#include <algorithm>
#include <complex>
#include <cstdlib>

#ifndef FUNCTIONS_H
#define FUNCTIONS_H

using namespace std;

double LCG(int, int, int, int, int);
double mean(vector <double>);
double var(vector <double>);
double two_point_corr(vector <double>, int);


class lattice
{
    /**
     * @brief Z2 lattice
     * 
     * @param N_1 is the horizontal dimension
     * @param N_2 is the vertical dimension
     * @param J is the energy constant of Ising model
     * @param Pbc determines periodic boundary conditions
     * @param grid is the actual lattice
     * 
     */
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
    /**
     * @brief 3d vector
     * 
     * @param x is first coordinate dimension
     * @param y is the second coordinate dimension
     * @param z is the third coordinate dimension
     * 
     */
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
    /**
     * @brief Physical particle
     * 
     * @param sigma is the diameter
     * @param mass is the mass
     * @param kinetic_en is the kinetic energy
     * @param vx is x velocity
     * @param vy is y velocity
     * @param vz is z velocity
     * @param fx is the x force acting on it
     * @param fy is the y force acting on it
     * @param fz is the z force acting on it
     * 
     */
    public:
    double sigma; 
    double mass; 
    double kinetic_en;
    double vx, vy, vz, fx, fy, fz;
    particle();
    vec to_vector();
    void get_kinetic_en();

};


class volume 
{
    /**
     * @brief Physical volume filled with particles
     * 
     * @param N is the number of particles
     * @param phi is the packing factor
     * @param rho is the density
     * @param L is the lenght of cubical box
     * @param T is the temperature of the system
     * @param kinetic_en is the kinetic energy of the configuration
     * @param Cv is the specific heat at constant volume
     * @param debuf is a file where to print debug logs
     * @param cdm is the center of mass of the particles
     * @param mom is the total momentum of the configuration
     * @param configuration is a pointer to the particle list
     * @param LookUpTable is a pointer to the look up table storing all inter-particle distances
     * 
     */
    public:
    unsigned N; 
    double phi; 
    double rho; 
    double L; 
    double T; 
    double kinetic_en; 
    double Cv;

    ofstream debug;
    vec cdm; 
    vec mom;
    particle *configuration; 
    double *LookUpTable; 
    

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
     * @brief Physical volume structured in cells
     * 
     * @param M_max is the maximum number of cells per coordinate dimension
     * @param rc is the lenght of cubical cell
     * @param potential is the potential energy of the system
     * @param cell_list is a pointer to the cell list storing all the indexes of prticles belonging to each cell
     * @param cell_LookUpTable is the inter-particle distance between particles in neighbouring cells
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
    void md_dynamics(const double&, const double&, const double&, const double&);
    complex<double> calculate_ssf(const vec&);

};

class MarkovSampler 
{
    /**
     * @brief Perform Monte Carlo simulation on a volume or cell_volume
     * 
     * @param moves is the number of MC sweeps
     * @param delta is the lenght of which each particle is moved at each MC step
     * 
     */
    public:
    unsigned moves;
    double delta;

    MarkovSampler();
    void __init__(const unsigned&, const double&);
    void CanonicalSampler(volume&);
    void CanonicalSampler(cell_volume&);
};

#endif // FUNCTIONS_H
