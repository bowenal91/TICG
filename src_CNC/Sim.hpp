#include "Objects.hpp"
#include "Box.hpp"
#include "mersenne.h"
#include "math_vector.h"
#include <vector>

using namespace LAMMPS_NS;

class Sim {
	public:
        Box box;
        Bead *beads;
        Grid *grid;
        Bond *bonds;
        Angle *angles;
        Twist *twists;
        CNC *cnc;
        Attachment *attachments;
        CRandomMersenne *RanGen;
   
        std::vector<std::vector<int> > bondList;
        std::vector<std::vector<int> > angleList;
        std::vector<std::vector<int> > twistList;
        std::vector<std::vector<int> > attachmentList;
        std::vector<std::vector<int> > moleculeList;

        double **nus;
        double ***ws;

        //double **chis;              //Matrix containing all binary chi interactions
        double **mus_f;             //Matrix containing nematic Maier-Saupe Parameters (Pi-stacking direction)
        double **mus_u;             //Matrix containing nematic Maier-Saupe Parameters (polymer alignment direction)

        Sim();
        void MCSim();
        void updateT();
        void openFiles();
        void read_input();
        void initialize_data();
        void initialize_system();
        void allocate_memory();
        void calc_random_vector(vector &b);
        void calc_random_normal_vector(vector &b, vector a);
        void calc_cross_vector(vector &c, vector a, vector b);
        double calc_GB_total();
        void update_disk_energies();
        double calc_disk_energy(int);
        double calc_disk_interval_energy(int,int);
        double calc_disk_GB_energy(int);
        double calc_disk_GB_energy_without(int,int);
        double calc_bond_total();
        double calc_wlc_total();
        double eps_func(double,double,double,double,double,double,double,double);
        double sig_func(double,double,double,double,double,double,double,double);
        double eps_func_asymm(double,double,double,double);
        double sig_func_asymm(double,double,double,double);
        bool withinExclude(int,int);
        double calc_dist(vector r1, vector r2);
        void nearest_image_dist(vector &r, vector r1, vector r2);
        void PBC_shift(vector &r_shift, vector r);
        void proj_vector(vector &a, vector n);
        void moveGrid();
        void initialize_grid();
        void build_bead_lists();

        int mapToGrid(LAMMPS_NS::vector r);
        int calc_grid_placement(int,int,int);
        void addBeadToGrid(int);
        void removeBeadFromGrid(int);

        void read_configuration();
        void read_CNC_data();
        void read_topology();
        void print_configuration();
        void generate_lists();
        void generate_matrices();
        double calc_total_energy();
        double calc_total_nonbond_energy();
        double calc_total_triplet_energy();
        double calc_bead_pair_energy(int);
        double calc_bead_triplet_energy(int);
        double calc_pair_energy(int,int);
        double calc_triplet_energy(int,int,int);
        double calc_bond_energy(int);
        double calc_angle_energy(int);
        double calc_twist_energy(int);
        double calc_site_energy(int);
        double calc_bead_energy(int);
        double calc_nonbond_bead_energy(int);
        double calc_grid_energy(int);
        double calc_nematic_pair_energy(int,int);
        double calc_molecule_energy(int);
        void update_all_attachments();
        void update_attachment(int);
        void update_atom_attachments(int);
        double calc_cos_angle(vector r0, vector r1, vector r2);
        void update_u_vectors(int);
        void update_single_u_vector(int);
        double calc_angle(vector r0, vector r1, vector r2);
        bool flagFunc(vector r);  



        void MCMove();
        void MC_displace();
        void MC_rotate();
        void MC_translate_molecule();
        void MC_rotate_molecule();
        void MC_pivot();
        void MC_reptate();


        void shiftCOM();
        void printXYZ();
        void printPSF();
        void printPOV();
        void writeEnergy();
        void dumpData();
        void dumpGrid();

        void printCheckpoint();
        void readCheckpoint();
        void empty_containers();
};
