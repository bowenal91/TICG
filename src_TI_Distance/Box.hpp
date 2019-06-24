#include "math_vector.h"

using namespace LAMMPS_NS;

class Box {
	public:
       double x; //Dimensions of the periodic domain
       double y;
       double z;

       double sigma; //Standard deviation of the gaussian density cloud representing each bead
       double invsigsqr;
       double prefact; //Prefactor for Gaussian that needs to be calculated only once
       double prefact3; //Prefactor for 3-body Gaussian Potential
       double rcut; //Cutoff distance for energy calculations
       int grid_cut; //Distance between grid cells that are cutoff

       double chi;          //Flory-Huggins parameter
       double kappa;        //Inverse compressibility
       double mu, K;        //Potential to describe pairwise nematic interactions
       double dens;         //Inverse density of the box
       double E_tot, E_nb, E_bond, E_wlc; //Total energy, non-bond energy, harmonic bond energy, twistable worm-like chain energy
       int seed; //Random number seed
       int numBeads; //Number of atoms in the simulation
       int numBonds; //Number of harmonic bonded interactions in the simulation
       int numAngles; //Number of angle interactions in the simulation
       int numTwists; //Number of angle interactions in the simulation
       int numAttachments; //Number of angle interactions in the simulation
       int numMolecules; //Number of molecules in the simulation. Used for large scale MC moves
       int numTypes;     //Number of different types of interacting particles in the simulation  

       int grid_width_x; //These things specify how far we have to search with each grid cell when computing the neighbor list
       int grid_width_y;
       int grid_width_z;
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

       double fracDisplace, fracTranslate, fracRotate, fracRotateMoly,fracPivot,fracReptate; //Relative occurence of the different Monte Carlo moves
       double displaceMax, rotateMax, translateMax, pivotMax;
       int reptateMax;

       int init_style; //1 for randomly placed but linear, 2 for crystal lattice, 3 for random walk placement
        
       bool hasWall; //There is a wall in the z-direction of the simulation
       
       int numCNC;
       double cylX, cylY, cylR; //x and y coordinate of cylinder center, and the cylinder radius

       double TI_shift; //Shift for calculating numerical derivative in thermodynamic integration
       double FPlus,FMinus; //Potential energies for increased and decreased distance, respectively. Used in calculating derivatives
       double TI_dF; //Differential of free energy
       int TI_sample_interval; //Interval over which free energy gradient is sampled
};
