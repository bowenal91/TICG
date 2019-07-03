#include "math_vector.h"

using namespace LAMMPS_NS;

class Box {
	public:
       double x; //Dimensions of the periodic domain
       double y;
       double z;

       double chi;          //Flory-Huggins parameter
       double kappa;        //Inverse compressibility
       double A,B,C;        //Phenomenological Landau-deGennes parameters used to describe nematic energy
       double mu, K;        //Potential to describe pairwise nematic interactions
       double Delta;        //Discretization length (in units of persistence length) of the WLC
       double dens;         //Inverse density of the box
       double E_tot, E_nb, E_bond, E_wlc; //Total energy, non-bond energy, harmonic bond energy, twistable worm-like chain energy
       int seed; //Random number seed
       int numBeads; //Number of atoms in the simulation
       int numBonds; //Number of harmonic bonded interactions in the simulation
       int numAngles; //Number of angle interactions in the simulation
       int numMolecules; //Number of molecules in the simulation. Used for large scale MC moves
       int numTypes;     //Number of different types of interacting particles in the simulation  

       double grid_min_x; //These three are all zero unless the grid is moving to eliminate PM0 artifacts
       double grid_min_y;
       double grid_min_z;
       int grid_num_x;    //These three are determined when the simulation is initialized
       int grid_num_y;
       int grid_num_z;
       double grid_dx;    //These three are the discretization of the system and are specified in the input file
       double grid_dy;
       double grid_dz;
       int numGrid;       //Number of total grid sites in the simulation

       int numCycles; //Number of Monte Carlo Cycles
       int numEq; //Number of equilibration cycles at high T
       int numAnneal; //Number of annealing cycles over which T is lowered
       int cycle; //Current cycle
       int numAccepts;
       int numRejects;
       int vizInterval; //Interval after which visualization data is outputted
       int checkpoint; //0 if starting from scratch. 1 if loading a configuration
       int checkInterval; //MC cycle interval after which checkpoints are made

       double fracDisplace, fracTranslate, fracRotate, fracPivot,fracReptate; //Relative occurence of the different Monte Carlo moves
       double displaceMax, rotateMax, translateMax, pivotMax;
       int reptateMax;

       int init_style; //1 for randomly placed but linear, 2 for crystal lattice, 3 for random walk placement
        
       bool hasWall; //There is a wall in the z-direction of the simulation
};
