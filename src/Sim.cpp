#include "time.h"
#include <iostream>
#include <math.h>
#include <cstdlib>
#include <stdio.h>
#include <algorithm>
#include "Sim.hpp"
#include "mersenne.h"
#define PI 3.14159265358979323846


using namespace LAMMPS_NS;
using namespace std;

Sim::Sim() {
    read_input();
    read_configuration();
    read_topology();
    generate_lists();
    generate_matrices();
    printXYZ();
    shiftCOM();
    initialize_grid();
    initialize_data();
    openFiles();
}

void Sim::MCSim() {
    box.numAccepts = 0;
    box.numRejects = 0;
    //printXYZ();
    double numNeighbors = (2*box.grid_width_x+1)*box.grid_dx*(2*box.grid_width_y+1)*box.grid_dy*(2*box.grid_width_z+1)*box.grid_dz*box.dens;
    printf("Average Number of Neighbors: %f\n",numNeighbors);
    for (box.cycle=0;box.cycle<box.numCycles;box.cycle++) {
        updateT();
        if (box.cycle%box.vizInterval == 0) {
            shiftCOM();
            printXYZ();
            dumpData();
            writeEnergy();
            printf("Turn: %d\n",box.cycle);

        }
        MCMove();
        if (box.cycle%box.checkInterval == 0) {
            print_configuration();
        }
    }
    writeEnergy();
    shiftCOM();
    dumpData();
    printXYZ();
    //printPOV();
    printf("Percent Acceptance: %f\n",100*double(box.numAccepts)/double(box.numAccepts+box.numRejects));
    empty_containers();
}

void Sim::openFiles() {
    FILE *dump;
    dump = fopen("config.xyz","w");
    fclose(dump);
    dump=fopen("Energy.txt","w");
    fclose(dump);
    
    dump = fopen("data.txt","w");
    fprintf(dump,"%d\n",box.numBeads);
    fclose(dump);

    dump = fopen("grid.txt","w");
    fprintf(dump,"%d\t%d\t%d\t%f\t%f\t%f\n",box.grid_num_x,box.grid_num_y,box.grid_num_z,box.grid_dx,box.grid_dy,box.grid_dz);
    fclose(dump);
    
    printPSF();

}

void Sim::read_input() {

  FILE *input;
  char tt[2001];
  int ret;
  input=fopen("input", "r");

  if (input!=NULL)
    { 	/******		Reading Input File		********/
      ret=fscanf(input,"%d", &box.seed);		     fgets(tt,2000,input);
      ret=fscanf(input,"%d", &box.checkpoint);		 fgets(tt,2000,input);
      ret=fscanf(input,"%d", &box.numCycles);		 fgets(tt,2000,input);
      ret=fscanf(input,"%d", &box.numEq);		 fgets(tt,2000,input);
      ret=fscanf(input,"%d", &box.numAnneal);		 fgets(tt,2000,input);
      ret=fscanf(input,"%d", &box.vizInterval);		 fgets(tt,2000,input);
      ret=fscanf(input,"%d", &box.checkInterval);    fgets(tt,2000,input);
      ret=fscanf(input,"%lf", &box.fracDisplace);    fgets(tt,2000,input);
      ret=fscanf(input,"%lf", &box.fracRotate);		 fgets(tt,2000,input);
      ret=fscanf(input,"%lf", &box.fracTranslate);   fgets(tt,2000,input);
      ret=fscanf(input,"%lf", &box.fracRotateMoly);   fgets(tt,2000,input);
      ret=fscanf(input,"%lf", &box.x);               fgets(tt,2000,input);
      ret=fscanf(input,"%lf", &box.y);               fgets(tt,2000,input);
      ret=fscanf(input,"%lf", &box.z);               fgets(tt,2000,input);
      ret=fscanf(input,"%lf", &box.sigma);               fgets(tt,2000,input);
      ret=fscanf(input,"%lf", &box.rcut);               fgets(tt,2000,input);
      ret=fscanf(input,"%lf", &box.grid_dx);               fgets(tt,2000,input);
      ret=fscanf(input,"%lf", &box.grid_dy);               fgets(tt,2000,input);
      ret=fscanf(input,"%lf", &box.grid_dz);               fgets(tt,2000,input);
      ret=fscanf(input,"%lf", &box.chi);               fgets(tt,2000,input);
      ret=fscanf(input,"%lf", &box.kappa);               fgets(tt,2000,input);
      ret=fscanf(input,"%lf", &box.mu);               fgets(tt,2000,input);
      ret=fscanf(input,"%lf", &box.K);               fgets(tt,2000,input);
      ret=fscanf(input,"%d",  &box.hasWall);               fgets(tt,2000,input);
      ret=fscanf(input,"%lf",  &box.displaceMax);               fgets(tt,2000,input);
      ret=fscanf(input,"%lf",  &box.translateMax);               fgets(tt,2000,input);
      ret=fscanf(input,"%lf",  &box.rotateMax);               fgets(tt,2000,input);
      ret=fscanf(input,"%lf",  &box.pivotMax);               fgets(tt,2000,input);
      fgets(tt,2000,input);fgets(tt,2000,input);

      fclose(input);
      /***************************************/
    }

  else
    { 
      fprintf(stdout,"\n Input File Not Found\n");
      exit(EXIT_FAILURE);
    }
    


}

void Sim::read_configuration() {
    int i,start,end,anchor,rigid,issite,hassites;
    char tt[2001];
    double total_volume = 0.0;
    

    FILE *atom_file = fopen("beads.dat","r"); //Location of all atom files 
    //fgets(tt,2000,atom_file);             //Section title
    fscanf(atom_file, "%d",&box.numBeads);
    fscanf(atom_file, "%d",&box.numTypes);
    beads = new Bead[box.numBeads];
    for (i=0;i<box.numBeads;i++) {
        fscanf(atom_file,"%d", &beads[i].id);
        fscanf(atom_file,"%lf%lf%lf", &beads[i].r[0], &beads[i].r[1], &beads[i].r[2]);
        fscanf(atom_file,"%lf%lf%lf", &beads[i].rn[0], &beads[i].rn[1], &beads[i].rn[2]);
        fscanf(atom_file,"%lf%lf%lf", &beads[i].u[0], &beads[i].u[1], &beads[i].u[2]);
        fscanf(atom_file,"%lf%lf%lf", &beads[i].f[0], &beads[i].f[1], &beads[i].f[2]);
        fscanf(atom_file,"%lf%lf%lf", &beads[i].v[0], &beads[i].v[1], &beads[i].v[2]);
        fscanf(atom_file,"%lf", &beads[i].vol);
        fscanf(atom_file,"%d", &beads[i].type);
        fscanf(atom_file,"%d", &anchor);
        fscanf(atom_file,"%d", &rigid);
        fscanf(atom_file,"%d", &issite);
        fscanf(atom_file,"%d", &hassites);
        if (anchor==1) {
            beads[i].isAnchor = true;
        } else {
            beads[i].isAnchor = false;
        }
        if (rigid==1) {
            beads[i].isRod = true;
        } else {
            beads[i].isRod = false;
        }
        if (issite==1) {
            beads[i].isSite = true;
        } else {
            beads[i].isSite = false;
        }
        if (hassites==1) {
            beads[i].hasSites = true;
        } else {
            beads[i].hasSites = false;
        }

        PBC_shift(beads[i].r,beads[i].rn); //Makes it easier when writing initial configurations

        total_volume += beads[i].vol;
    }
   
    fclose(atom_file);
     
    box.dens = total_volume/(box.x*box.y*box.z);
}

void Sim::read_topology() {
    int i;
    int start,end;
    char tt[2001];
    FILE *topo_file = fopen("topology.dat","r");

    fgets(tt,2000,topo_file);             //BONDS
    fscanf(topo_file,"%d",&box.numBonds);
    bonds = new Bond[box.numBonds];
    for (i=0;i<box.numBonds;i++) {
        fscanf(topo_file,"%d%d%lf%lf",&bonds[i].id1, &bonds[i].id2, &bonds[i].k_G, &bonds[i].r0);
    }
    
    fgets(tt,2000,topo_file);             //ANGLES
    fgets(tt,2000,topo_file);             //ANGLES
    fscanf(topo_file,"%d",&box.numAngles);
    angles = new Angle[box.numAngles];
    for (i=0;i<box.numAngles;i++) {
        fscanf(topo_file,"%d%d%d%lf%lf",&angles[i].id0,&angles[i].id1,&angles[i].id2,&angles[i].k_b,&angles[i].theta0);
    }

    fgets(tt,2000,topo_file);             //TWISTS
    fgets(tt,2000,topo_file);             //TWISTS
    fscanf(topo_file,"%d",&box.numTwists);
    twists = new Twist[box.numTwists];
    for (i=0;i<box.numTwists;i++) {
        fscanf(topo_file,"%d%d%d%d%lf%lf%lf%lf",&twists[i].id1,&twists[i].id2, &twists[i].id3, &twists[i].id4, &twists[i].a1, &twists[i].a2, &twists[i].a3, &twists[i].a4);
    }

    fgets(tt,2000,topo_file);             //ATTACHMENTS
    fgets(tt,2000,topo_file);             //ATTACHMENTS
    fscanf(topo_file,"%d",&box.numAttachments);
    attachments = new Attachment[box.numAttachments];
    for (i=0;i<box.numAttachments;i++) {
        fscanf(topo_file,"%d%d%lf%lf%lf",&attachments[i].id_core,&attachments[i].id_site,&attachments[i].u_comp,&attachments[i].f_comp,&attachments[i].v_comp);
    }



    fgets(tt,2000,topo_file);             //MOLECULES
    fgets(tt,2000,topo_file);             //MOLECULES
    fscanf(topo_file,"%d",&box.numMolecules);
    moleculeList.clear();
    moleculeList.resize(box.numMolecules);
    for (i=0;i<box.numMolecules;i++) {
        fscanf(topo_file,"%d%d",&start,&end); //Start and end for atoms
        moleculeList[i].push_back(start);
        moleculeList[i].push_back(end);
        fscanf(topo_file,"%d%d",&start,&end); //Start and end for bonds
        moleculeList[i].push_back(start);
        moleculeList[i].push_back(end);
        fscanf(topo_file,"%d%d",&start,&end); //Start and end for angles
        moleculeList[i].push_back(start);
        moleculeList[i].push_back(end);     
        fscanf(topo_file,"%d%d",&start,&end); //Start and end for twists
        moleculeList[i].push_back(start);
        moleculeList[i].push_back(end);     

    }

    fclose(topo_file);
}

void Sim::print_configuration() {
    int i;
    char tt[2001];
    int start,end,anchor,rigid,issite,hassites;
    FILE *atom_file = fopen("beads.dat","w"); //Location of all atom files 
    fprintf(atom_file, "%d\n",box.numBeads);
    fprintf(atom_file, "%d\n",box.numTypes);
    for (i=0;i<box.numBeads;i++) {
        fprintf(atom_file,"%d\t", beads[i].id);
        fprintf(atom_file,"%lf\t%lf\t%lf\t", beads[i].r[0], beads[i].r[1], beads[i].r[2]);
        fprintf(atom_file,"%lf\t%lf\t%lf\t", beads[i].rn[0], beads[i].rn[1], beads[i].rn[2]);
        fprintf(atom_file,"%lf\t%lf\t%lf\t", beads[i].u[0], beads[i].u[1], beads[i].u[2]);
        fprintf(atom_file,"%lf\t%lf\t%lf\t", beads[i].f[0], beads[i].f[1], beads[i].f[2]);
        fprintf(atom_file,"%lf\t%lf\t%lf\t", beads[i].v[0], beads[i].v[1], beads[i].v[2]);
        if (beads[i].isAnchor) {
            anchor = 1;
        } else {
            anchor = 0;
        }
        if (beads[i].isRod) {
            rigid = 1;
        } else {
            rigid = 0;
        }
        if (beads[i].isSite) {
            issite = 1;
        } else {
            issite = 0;
        }
        if (beads[i].hasSites) {
            hassites = 1;
        } else {
            hassites = 0;
        }
        fprintf(atom_file,"%lf\t", beads[i].vol);
        fprintf(atom_file,"%d\t", beads[i].type);
        fprintf(atom_file,"%d\t", anchor);
        fprintf(atom_file,"%d\t", rigid);
        fprintf(atom_file,"%d\t", issite);
        fprintf(atom_file,"%d\n", hassites);
        vec_norm(beads[i].u);
    }

    fclose(atom_file);
}

int Sim::mapToGrid(LAMMPS_NS::vector r) {
    double x,y,z,x1,y1,z1;
    int a,b,c,m;
    x = r[0];
    y = r[1];
    z = r[2];


    a = floor(x/box.grid_dx);
    b = floor(y/box.grid_dy);
    c = floor(z/box.grid_dz);
    
    m = c*box.grid_num_x*box.grid_num_y + b*box.grid_num_x + a;

    return m;
}

void Sim::initialize_grid() {
    int i,j,k,m;
    //Must always choose the larger grid spacing to make sure that nothing goes horribly wrong
    box.grid_num_x = int(floor(box.x / box.grid_dx));
    box.grid_num_y = int(floor(box.y / box.grid_dy));
    box.grid_num_z = int(floor(box.z / box.grid_dz));

    box.grid_dx = box.x / double(box.grid_num_x);
    box.grid_dy = box.y / double(box.grid_num_y);
    box.grid_dz = box.z / double(box.grid_num_z);

    box.grid_width_x = int(ceil(box.rcut/box.grid_dx));
    box.grid_width_y = int(ceil(box.rcut/box.grid_dy));
    box.grid_width_z = int(ceil(box.rcut/box.grid_dz));

    //printf("%f\t%f\t%f\n",box.grid_dx,box.grid_dy,box.grid_dz);

    box.numGrid = box.grid_num_x*box.grid_num_y*box.grid_num_z;

    grid = new Grid[box.numGrid];
    for (k=0;k<box.grid_num_z;k++) {
        for (j=0;j<box.grid_num_y;j++) {
            for (i=0;i<box.grid_num_x;i++) {
                m = calc_grid_placement(i,j,k); 
                grid[m].x = i;
                grid[m].y = j;
                grid[m].z = k;
                grid[m].id = m;
            }
        }
    }
    
    build_bead_lists();

}

void Sim::build_bead_lists() {
    int i,m;
    for (i=0;i<box.numGrid;i++) {
        grid[i].bead_list.head = NULL;
    }
    
    for (i=0;i<box.numBeads;i++) {
        beads[i].next = NULL;
        beads[i].prev = NULL;
    }

    for (i=0;i<box.numBeads;i++) {
        m = mapToGrid(beads[i].r);
        beads[i].grid = m;
        grid[m].bead_list.pushNode(&beads[i]);
    }
    
}

void Sim::generate_lists() { //Generate the lists containing all the bonds/angles/twists associated with each atom
    int i;
    
    bondList.clear();
    angleList.clear();
    twistList.clear();
    attachmentList.clear();
    
    bondList.resize(box.numBeads);
    angleList.resize(box.numBeads);
    twistList.resize(box.numBeads);
    attachmentList.resize(box.numBeads);

    for (i=0;i<box.numBonds;i++) { //Element N of bondList contains a list of id's for all bonds that atom N is a part of
        bondList[bonds[i].id1].push_back(i);
        bondList[bonds[i].id2].push_back(i);
    }

    for (i=0;i<box.numAngles;i++) {
        angleList[angles[i].id0].push_back(i);
        angleList[angles[i].id1].push_back(i);
        angleList[angles[i].id2].push_back(i);
    }

    for (i=0;i<box.numTwists;i++) {
        twistList[twists[i].id1].push_back(i);
        if (twists[i].id1 != twists[i].id2)
            twistList[twists[i].id2].push_back(i);
        twistList[twists[i].id3].push_back(i);
        if (twists[i].id3 != twists[i].id4)
            twistList[twists[i].id4].push_back(i);
    }

    for (i=0;i<box.numAttachments;i++) {
        attachmentList[attachments[i].id_core].push_back(i);
        attachmentList[attachments[i].id_site].push_back(i);
    }



}


void Sim::generate_matrices() {
    int i,j;
    
    chis = new double *[box.numTypes];
    mus_f = new double *[box.numTypes];
    mus_u = new double *[box.numTypes];

    for (i=0;i<box.numTypes;i++) {
        
        chis[i] = new double [box.numTypes];
        mus_f[i] = new double [box.numTypes];
        mus_u[i] = new double [box.numTypes];

    }
   
    FILE *data = fopen("interactions.dat","r");
    for (i=0;i<box.numTypes;i++) {
        for (j=0;j<box.numTypes;j++) {
            fscanf(data,"%lf",&chis[i][j]);
        }
    }
    for (i=0;i<box.numTypes;i++) {
        for (j=0;j<box.numTypes;j++) {
            fscanf(data,"%lf",&mus_f[i][j]);
        }
    }
    for (i=0;i<box.numTypes;i++) {
        for (j=0;j<box.numTypes;j++) {
            fscanf(data,"%lf",&mus_u[i][j]);
        }
    }

    fclose(data);

}


void Sim::initialize_data() {
    double dummy = sqrt(4.0*PI*box.sigma*box.sigma);
    dummy = dummy*dummy*dummy;
    box.prefact = 1.0/(box.dens*dummy);
    box.invsigsqr = 1.0/(box.sigma*box.sigma);
    
    if (box.seed < 0) {
        box.seed = time(NULL);
    }
    RanGen = new CRandomMersenne(box.seed);
    box.E_tot = calc_total_energy(); 
}

void Sim::calc_random_vector(LAMMPS_NS::vector &b) {
    double R1,R2,R3;
    do { //Generate random unit vector and make it u.
        R1 = (2*RanGen->Random()-1);
        R2 = (2*RanGen->Random()-1);
        R3 = R1*R1+R2*R2;
    } while (R3>=1);
               
    b[0] = 2*sqrtl(1.0-R3)*R1;
    b[1] = 2*sqrtl(1.0-R3)*R2;
    b[2] = 1-2.0*R3;


}

void Sim::calc_random_normal_vector(LAMMPS_NS::vector &n, LAMMPS_NS::vector a) { //Uses a rotation matrix to calculate a unit vector, n, normal to input vector a
    double theta, phi, A,B,C,D,R1,x0,y0;
    double ax,ay,az;
    ax = a[0]; ay = a[1]; az = a[2];
    theta = acos(az / sqrt(ax*ax+ay*ay+az*az));
    phi = atan2(ay,ax);
    A = cos(theta);
    B = sin(theta);
    C = cos(phi);
    D = sin(phi);
    R1 = 2*PI*RanGen->Random();
    x0 = cos(R1);
    y0 = sin(R1);
    n[0] = A*C*x0 - D*y0;
    n[1] = A*D*x0 + C*y0;
    n[2] = -1*B*x0;

    
}

void Sim::calc_cross_vector(LAMMPS_NS::vector &c, LAMMPS_NS::vector a, LAMMPS_NS::vector b) {//Cross product
    
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];

}

void Sim::updateT() {
    /*
    int i = box.cycle;
    double dT = (box.T_anneal-box.T_fin)/double(box.numAnneal);
    int numSteps;
    if (i < box.numEq) {
        box.T = box.T_anneal;
    } else if (i > box.numEq && i < box.numEq+box.numAnneal) {
        numSteps = i - box.numEq;
        box.T = box.T_anneal - numSteps*dT;
    } else {
        box.T = box.T_fin;
    }
    */
}

/*********************************ENERGY CALCULATIONS**************************************/
double Sim::calc_total_energy() {
    int i,j;
    double E = 0.0;
    double grid,bond,angle; 
    
    E += calc_total_nonbond_energy();
    
    //grid = E;
    for (i=0;i<box.numBonds;i++) {
        E += calc_bond_energy(i);
    }
    //bond = E - grid;
    for (i=0;i<box.numAngles;i++) {
        E += calc_angle_energy(i);
    }
    
    for (i=0;i<box.numTwists;i++) {
        E += calc_twist_energy(i);
    }



    //angle = E - grid - bond;
    //printf("%f\t%f\t%f\n",grid/double(box.numBeads),bond/double(box.numBonds),angle/double(box.numAngles));
    return E;
    
}

double Sim::calc_total_nonbond_energy() { //Loops through all grid cells and calculates pair interactions only for cells that are larger than current cell. This avoids double counting
    int i,j,k,ii,jj,kk,m1,m2;
    double E = 0.0;
    for (i=0;i<box.numBeads;i++) {
        E += calc_nonbond_bead_energy(i);
    }   

    E *= 0.5;


    return E;

}

double Sim::calc_bond_energy(int i) {
    double r,k_G,r0,E;
    int id1,id2;
    int j;

    k_G = bonds[i].k_G; 
    r0 = bonds[i].r0;
    id1 = bonds[i].id1;
    id2 = bonds[i].id2;
    r = calc_dist(beads[id1].rn,beads[id2].rn);

    E = k_G*(r-r0)*(r-r0);

    return E;
}

double Sim::calc_angle_energy(int i) {
    double k, theta0, theta, E;
    int id0,id1,id2;

    k = angles[i].k_b;
    theta0 = angles[i].theta0;
    id0 = angles[i].id0;
    id1 = angles[i].id1;
    id2 = angles[i].id2;
    
    theta = calc_angle(beads[id0].rn,beads[id1].rn,beads[id2].rn);

    E = 0.5*k*(theta-theta0)*(theta-theta0);
    return E;
}

double Sim::calc_twist_energy(int i) {
    double a1,a2,a3,a4, theta, E;
    int id1,id2,id3,id4,j;
    double value;
    double c11,c22,c33,c12,c13,c23,d12,d13,d23,hi,lo,costheta;

    a1 = twists[i].a1;
    a2 = twists[i].a2;
    a3 = twists[i].a3;
    a4 = twists[i].a4;
    id1 = twists[i].id1;
    id2 = twists[i].id2;
    id3 = twists[i].id3;
    id4 = twists[i].id4;
    LAMMPS_NS::vector r1, r2, r3; //4-body potential between points a,b,c,d. r1 = b-a, r2 = c-b, r3 = d-c.
   
    //Set up each arm of the dihedral, and calculate dot products and theta
    
    if (id1==id2) {
        for (j=0;j<3;j++) {
            r1[j] = -beads[id1].f[j];
        }
    } else {
        for (j=0;j<3;j++) {
            r1[j] = beads[id2].rn[j] - beads[id1].rn[j];
        }
    }

    for (j=0;j<3;j++) {
        r2[j] = beads[id3].rn[j] - beads[id2].rn[j];
    }

    if (id3==id4) {
        for (j=0;j<3;j++) {
            r3[j] = beads[id4].f[j];
        }
    } else {
        for (j=0;j<3;j++) {
            r3[j] = beads[id4].rn[j] - beads[id3].rn[j];
        }
    }

    //Dot products, etc
    c11 = vec_dot(r1,r1);
    c22 = vec_dot(r2,r2);
    c33 = vec_dot(r3,r3);
    c12 = vec_dot(r1,r2);
    c13 = vec_dot(r1,r3);
    c23 = vec_dot(r2,r3);

    d12 = c11*c22 - c12*c12;
    d13 = c11*c33 - c13*c13;
    d23 = c22*c33 - c23*c23;

    hi = -(c23*c12-c13*c22);
    lo = sqrt(d23*d12);

    costheta = hi/lo;

    if (costheta > 1.0) costheta = 1.0;
    if (costheta < -1.0) costheta = -1.0;

    theta = acos(costheta);

    E = 0.5*(a1*(1+cos(theta)) + a2*(1-cos(2*theta)) + a3*(1+cos(3*theta)) + a4*(1-cos(4*theta)));
/*
    if (i==0) {
        printf("%f\t%f\n",theta,E);
    }
    */
    return E;


}

double Sim::calc_site_energy(int i) { //Calculate the energy associated with a site that is not associated with its core atom
    int j,k,core;
    k = attachmentList[i][0]; //Each site should only ever have one attachment - otherwise it would be overspecified
    core = attachments[k].id_core;
    double E = 0.0;

    for (j=0;j<bondList[i].size();j++) {
        k = bondList[i][j];
        if (bonds[k].id1 != core && bonds[k].id2 != core) {
            E += calc_bond_energy(k);
        }
    }

    for (j=0;j<angleList[i].size();j++) {
        k = angleList[i][j];
        if (angles[k].id0 != core && angles[k].id1 != core && angles[k].id2 != core) {
            E += calc_angle_energy(k);
        }
    }

    for (j=0;j<twistList[i].size();j++) {
        k = twistList[i][j];
        if (twists[k].id1 != core && twists[k].id2 != core) {
            E += calc_twist_energy(k);
        }
    }
    
    return E;

}

double Sim::calc_bead_energy(int i) { //Calculate all non-grid energy associated with a particular bead - this energy is in no way dependent on the u orientation vector
    int j,k;
    double E = 0.0;
    for (j=0;j<bondList[i].size();j++) {
        k = bondList[i][j];
        E += calc_bond_energy(k);
    }

    for (j=0;j<angleList[i].size();j++) {
        k = angleList[i][j];
        E += calc_angle_energy(k);
    }
    
    for (j=0;j<twistList[i].size();j++) {
        k = twistList[i][j];
        E += calc_twist_energy(k);
    }
    
    E += calc_nonbond_bead_energy(i);
    return E;
}

int Sim::calc_grid_placement(int a, int b, int c) {
    a = a%box.grid_num_x;
    if (a<0) a += box.grid_num_x;
    b = b%box.grid_num_y;
    if (b<0) b += box.grid_num_y;
    c = c%box.grid_num_z;
    if (c<0) c += box.grid_num_z;
    
    return c*box.grid_num_y*box.grid_num_x + b*box.grid_num_x + a;
}

double Sim::calc_nonbond_bead_energy(int i) {
    int ii,jj,kk,m,x,y,z,j;
    int zmin,zmax;
    m = beads[i].grid;
    x = grid[m].x;
    y = grid[m].y;
    z = grid[m].z;
    int iter = 0;
    Bead *neighbor;

    if (box.hasWall == 1) {
        zmax = min(z+box.grid_width_z,box.grid_num_z-1);
        zmin = max(z-box.grid_width_z,0);
    } else {
        zmin = z-box.grid_width_z;
        zmax = z+box.grid_width_z;
    }
    
    double E = 0.0;

    for (ii=x-box.grid_width_x;ii<=x+box.grid_width_x;ii++) {
        for (jj=y-box.grid_width_y;jj<=y+box.grid_width_y;jj++) {
            for (kk=zmin;kk<=zmax;kk++) {
                m = calc_grid_placement(ii,jj,kk);
                neighbor = grid[m].bead_list.head;
                while(neighbor != NULL) {
                    j = neighbor->id;
                    if (i!=j) E += calc_pair_energy(i,j);
                    neighbor = neighbor->next;
                    iter++;
                }
            }
        }
    }

    return E;
     
}

double Sim::calc_pair_energy(int i, int j) {
    double E,delta,dist;
    double E1;
    int t1,t2;
    LAMMPS_NS::vector rnear;
    t1 = beads[i].type;
    t2 = beads[j].type;

    if (t1 == t2) {
        delta = 1.0;
    } else {
        delta = 0.0;
    }

    nearest_image_dist(rnear,beads[i].r,beads[j].r);

    dist = vec_dot(rnear,rnear);

    if(box.rcut*box.rcut < dist) {
        return 0.0;
    }
    
    E = box.prefact*exp(-0.25*dist*box.invsigsqr);
    
    E *= box.kappa + chis[t1][t2]*(1.0+delta) + 2.0*calc_nematic_pair_energy(i,j);
    
    return E;

}


double Sim::calc_grid_energy(int i) { //Calculate all energy associated with a particular grid - NEED TO UPDATE THIS WHEN CONSIDERING MORE THAN BINARY INTERACTIONS
/*
    double E,E1,E2;
    double vol  = grid[i].vol;
    double phi_tot = 0.0;
    int t1,t2;
    for (t1=0;t1<box.numTypes;t1++) {
        phi_tot += grid[i].phi[t1];
    }

    E1 = 0.5*box.kappa*phi_tot*phi_tot;
    
    for (t1=0;t1<box.numTypes;t1++) {
        for (t2=t1+1;t2<box.numTypes;t2++) {
            E1 += chis[t1][t2]*grid[i].phi[t1]*grid[i].phi[t2];
        }
    }

    E1 *= box.dens*vol;

    E2 = 0.0; 
    int j,k;
    //printf("%d\n",grid[i].bead_list.size());
    for (j=0;j<grid[i].bead_list.size();j++) {
        for (k=j+1;k<grid[i].bead_list.size();k++) {
            E2 += calc_nematic_pair_energy(grid[i].bead_list[j],grid[i].bead_list[k]);
        }
    }
    if (vol == 0.0) {
        E2 = 0.0;
    } else {
        E2 /= (box.dens*vol);
    }
    E = E1+E2;
    //printf("%f\t%f\t%f\t%f\t%f\n",phiA,phiB,E1,E2,E);

    return E;
    */
}

double Sim::calc_nematic_pair_energy(int i, int j) {
    double dot1, dot2, E;
    int t1,t2;
    if (!beads[i].isRod || !beads[j].isRod) return 0.0; 
    dot1 = vec_dot(beads[i].f,beads[j].f);
    dot2 = vec_dot(beads[i].u,beads[j].u);

    t1 = beads[i].type;
    t2 = beads[j].type;
    
    E = -mus_f[t1][t2] * dot1*dot1 - mus_u[t1][t2] * dot2*dot2;

    return E;
    
    
    //LAMMPS_NS::vector rij, cross;
    //nearest_image_dist(rij, beads[i].r, beads[j].r);
    //vec_norm(rij);
    //calc_cross_vector(cross, beads[i].u, beads[j].u);
    //dot2 = vec_dot(cross, rij);

    //E = -box.mu * dot1*dot1 + box.K*dot2*dot1;

}

double Sim::calc_molecule_energy(int m) { //Need to figure this crap out later on - MIGHT NOT BE PARTICULARLY USEFUL FOR GRID ENERGY MODEL
    double E = 0.0;
    int i,j,k;
    int start_atom,end_atom,start_bond,end_bond,start_angle,end_angle,start_twist,end_twist;

    start_atom = moleculeList[m][0];
    end_atom = moleculeList[m][1];
    start_bond = moleculeList[m][2];
    end_bond = moleculeList[m][3];
    start_angle = moleculeList[m][4];
    end_angle = moleculeList[m][5];

    if (end_bond-start_bond != 0) { 
        for (i=start_bond;i<=end_bond;i++) {
            E += calc_bond_energy(i);
        }
    }

    if (end_angle-start_angle != 0) {
        for (i=start_angle;i<=end_angle;i++) {
            E += calc_angle_energy(i);
        }
    }
    
    return E;

}

void Sim::update_all_attachments() {
    int i;
    for (i=0;i<box.numAttachments;i++) {
        update_attachment(i);
    }
}

void Sim::update_attachment(int i) {//Shift the site indicated by the attachment to be oriented correctly according to the core
    int core,site;
    double u,f,v;
    core = attachments[i].id_core;
    site = attachments[i].id_site;
    u = attachments[i].u_comp;
    f = attachments[i].f_comp;
    v = attachments[i].v_comp;

    beads[site].rn[0] = beads[core].rn[0] + u*beads[core].u[0] + f*beads[core].f[0] + v*beads[core].v[0];
    beads[site].rn[1] = beads[core].rn[1] + u*beads[core].u[1] + f*beads[core].f[1] + v*beads[core].v[1];
    beads[site].rn[2] = beads[core].rn[2] + u*beads[core].u[2] + f*beads[core].f[2] + v*beads[core].v[2];

    PBC_shift(beads[site].r,beads[site].rn);

}


void Sim::update_atom_attachments(int i) { //Update all attachments associated with a given core atom
    int j;
    if (attachmentList.size() != 0) {

        for (j=0;j<attachmentList[i].size();j++) {
            update_attachment( attachmentList[i][j] );
        }

    }

}



/*************************************************************************************************/
/*
void Sim::update_u_vectors(int i) { //Update all u-vectors that are dependent on a bead's position: i, i-1, i+1 
    
    int j;

    if (beads[i].isStart) { //Bead is the start of a chain - forward difference
        
        update_single_u_vector(i);
        update_single_u_vector(i+1);

    } else if (beads[i].isEnd) { //Bead is end of a chain - backward difference

        update_single_u_vector(i);
        update_single_u_vector(i-1);

    } else { //Bead is in the middle of a chain - central difference

        update_single_u_vector(i+1);
        update_single_u_vector(i-1);

    }

}

void Sim::update_single_u_vector(int i) {   //Defines the u vector as the central difference of the two nearest connecting beads
    int j;

    if (beads[i].isStart) { //Bead is the start of a chain - forward difference

        for (j=0;j<3;j++) {
            beads[i].u[j] = beads[i+1].rn[j] - beads[i].rn[j];
        }

        vec_norm(beads[i].u);

    } else if (beads[i].isEnd) { //Bead is end of a chain - backward difference

        for (j=0;j<3;j++) {
            beads[i].u[j] = beads[i].rn[j] - beads[i-1].rn[j];
        }

        vec_norm(beads[i].u);

    } else { //Bead is in the middle of a chain - central difference

        for (j=0;j<3;j++) {
            beads[i].u[j] = beads[i+1].rn[j] - beads[i-1].rn[j];
        }

        vec_norm(beads[i].u);

    }

}
*/
double Sim::calc_angle(LAMMPS_NS::vector r0, LAMMPS_NS::vector r1, LAMMPS_NS::vector r2) { //Calculate and return the angle between vectors r2-r0 and r1-r0
    double theta;
    LAMMPS_NS::vector v1,v2;
    v1[0] = r1[0] - r0[0]; 
    v1[1] = r1[1] - r0[1]; 
    v1[2] = r1[2] - r0[2]; 
    v2[0] = r2[0] - r0[0]; 
    v2[1] = r2[1] - r0[1]; 
    v2[2] = r2[2] - r0[2]; 

    vec_norm(v1);
    vec_norm(v2);

    theta = vec_dot(v1,v2);
    if (theta > 1.0) theta = 1.0;
    if (theta < -1.0) theta = -1.0;
    theta = acos(theta);
    return theta;

}

double Sim::calc_cos_angle(LAMMPS_NS::vector r0, LAMMPS_NS::vector r1, LAMMPS_NS::vector r2) { //Calculate and return the angle between vectors r2-r0 and r1-r0
    double costheta;
    LAMMPS_NS::vector v1,v2;
    v1[0] = r1[0] - r0[0]; 
    v1[1] = r1[1] - r0[1]; 
    v1[2] = r1[2] - r0[2]; 
    v2[0] = r2[0] - r0[0]; 
    v2[1] = r2[1] - r0[1]; 
    v2[2] = r2[2] - r0[2]; 

    vec_norm(v1);
    vec_norm(v2);

    costheta = vec_dot(v1,v2);
    if (costheta > 1.0) costheta = 1.0;
    if (costheta < -1.0) costheta = -1.0;
    return costheta;

}

double Sim::calc_dist(LAMMPS_NS::vector r1, LAMMPS_NS::vector r2) {

    double x,y,z,r;
    x = r1[0] - r2[0];
    y = r1[1] - r2[1];
    z = r1[2] - r2[2];

    r = sqrt(x*x+y*y+z*z);
    return r;

}

void Sim::nearest_image_dist(LAMMPS_NS::vector &r, LAMMPS_NS::vector r1, LAMMPS_NS::vector r2) { //r1 - r2 using the nearest image convention
    
    double dx,dy,dz;
    double xhalf = box.x/2;
    double yhalf = box.y/2;
    double zhalf = box.z/2;

    dx = r1[0] - r2[0];
    dy = r1[1] - r2[1];
    dz = r1[2] - r2[2];

    dx = dx > xhalf ? dx - box.x : dx; 
    dx = dx < -xhalf ? dx + box.x : dx; 
    dy = dy > yhalf ? dy - box.y : dy; 
    dy = dy < -yhalf ? dy + box.y : dy; 
    dz = dz > zhalf ? dz - box.z : dz; 
    dz = dz < -zhalf ? dz + box.z : dz; 

    r[0] = dx; r[1] = dy; r[2] = dz;

}

void Sim::PBC_shift(LAMMPS_NS::vector &r_shift, LAMMPS_NS::vector r) {
    r_shift[0] = r[0] - box.x*floor(r[0]/box.x);
    r_shift[1] = r[1] - box.y*floor(r[1]/box.y);
    r_shift[2] = r[2] - box.z*floor(r[2]/box.z);
}

void Sim::proj_vector(LAMMPS_NS::vector &a, LAMMPS_NS::vector n) {//Project vector a onto the plane described by unit normal vector n, then normalize the new a
    double dot;
    int i;
    LAMMPS_NS::vector n2;
    for (i=0;i<3;i++) {
        n2[i] = n[i];
    }
    dot = vec_dot(a,n2);
    for (i=0;i<3;i++) {
        a[i] = a[i] - dot*n2[i];
    }
    vec_norm(a);
}

/***************************************************MONTE CARLO MOVES******************************************/

void Sim::MCMove() {
    double displace = box.fracDisplace;
    double rotate = box.fracDisplace+box.fracRotate;
    double translate = box.fracDisplace+box.fracRotate+box.fracTranslate;
    double rotate_moly = box.fracDisplace+box.fracRotate+box.fracTranslate+box.fracRotateMoly;
   
    double rand = RanGen->Random();

    if (rand < displace) 
        MC_displace();
    else if (rand < rotate) 
        MC_rotate();
    else if (rand < translate) 
        MC_translate_molecule();
    else if (rand < rotate_moly)
        MC_rotate_molecule();

}

//Philosophy: Get indices of all grids affected by the move before making the move. Do not double count affected grids. Need to watch grids that are affected by orientation changes

void Sim::MC_displace() { //Randomly displace all beads in the simulation
    int i,j,k,first,last,m,m2;
    double rx,ry,rz;
    double delta = box.displaceMax; //Maximum displacement distance in any direction
    double pre_energy, post_energy, dE;
    LAMMPS_NS::vector rnew,rper;
    Bead save;
    double a,b;
    bool flag;
    
    for (int move=0;move<box.numBeads;move++) {

        do {
            i = double(box.numBeads)*RanGen->Random();
        } while(i >= box.numBeads);
       
        if (beads[i].isAnchor) continue; //If a bead is attached to a wall then it will not be moved
        if (beads[i].isSite) continue; //If the bead has rigid bonds, then a displace move will not be accepted

        m = beads[i].grid;

        //Figure out which beads and grids are affected and store that data

        rx = (2.0*RanGen->Random()-1.0)*delta;
        ry = (2.0*RanGen->Random()-1.0)*delta;
        rz = (2.0*RanGen->Random()-1.0)*delta;

        save = beads[i];

        rnew[0] = beads[i].rn[0] + rx;
        rnew[1] = beads[i].rn[1] + ry;
        rnew[2] = beads[i].rn[2] + rz;
     
        if (box.hasWall) {
            if(rnew[2] > box.z || rnew[2] < 0.0) {
                continue; //Automatically reject any moves that push it through a wall
            }
        }

        PBC_shift(rper, rnew);  //Gotta figure out what to do with attached sites

        m2 = mapToGrid(rper);
        pre_energy = calc_bead_energy(i);
        if (beads[i].hasSites) {
            for (j=0;j<attachmentList[i].size();j++) {
                k = attachmentList[i][j];
                pre_energy += calc_site_energy(attachments[k].id_site);
            }
        }
        grid[m].bead_list.deleteNode(&beads[i]);

        //Data storage and energy calculations complete - now perform move

        for (j=0;j<3;j++) {
            beads[i].rn[j] = rnew[j];
            beads[i].r[j] = rper[j];
        }

        if (beads[i].hasSites) {
            update_atom_attachments(i);
        }
       
        beads[i].grid = m2;
        grid[m2].bead_list.pushNode(&beads[i]);

        //Move complete now perform post move energy calculations
        post_energy = calc_bead_energy(i);
        if (beads[i].hasSites) {
            for (j=0;j<attachmentList[i].size();j++) {
                k = attachmentList[i][j];
                post_energy += calc_site_energy(attachments[k].id_site);
            }
        }

        dE = post_energy - pre_energy;
        if (RanGen->Random() < exp(-dE)) { //Accept changes
            box.E_tot += dE;
            box.numAccepts++;
        } else { //Reverse changes
            //Restore Beads
            
            grid[m2].bead_list.deleteNode(&beads[i]);
            beads[i] = save;
            if (beads[i].hasSites) {
                update_atom_attachments(i);
            }
            grid[m].bead_list.pushNode(&beads[i]);
            //Restore Grids
            box.numRejects++;
        }
        
        
    }
    
}

void Sim::MC_rotate() {
    int i,j,k,flag,m;
    double max_theta = box.rotateMax;
    double theta;
    double pre_energy,post_energy,dE;
    Bead save;
    LAMMPS_NS::vector axis;
    quaternion quat;

    for (int move = 0;move<box.numBeads;move++) {
        do {
            i = double(box.numBeads)*RanGen->Random();
        } while(i >= box.numBeads);
       
        if (!beads[i].hasSites) continue; //If the bead does not have any sites, then rotating it will not do much 
        if (beads[i].isAnchor) continue; //If a bead is attached to a wall then it will not be moved
        if (beads[i].isSite) continue; //If the bead has rigid bonds, then a displace move will not be accepted

        m = beads[i].grid;
        
        save = beads[i];
        //No need to do anything with the neighborlists because the bead does not move
        theta = max_theta*(2*RanGen->Random()-1);
        calc_random_vector(axis);
        quat[0] = cos(theta/2);
        for (j=0;j<3;j++) {
            quat[j+1] = axis[j]*sin(theta/2);
        }

        pre_energy = calc_bead_energy(i);
        for (j=0;j<attachmentList[i].size();j++) {
            k = attachmentList[i][j];
            pre_energy += calc_site_energy(attachments[k].id_site);
        }

        quat_vec_rot(beads[i].u,save.u,quat);
        quat_vec_rot(beads[i].f,save.f,quat);
        quat_vec_rot(beads[i].v,save.v,quat);

        update_atom_attachments(i);


        post_energy = calc_bead_energy(i);
        for (j=0;j<attachmentList[i].size();j++) {
            k = attachmentList[i][j];
            post_energy += calc_site_energy(attachments[k].id_site);
        }
        dE = post_energy - pre_energy;
        
        if (RanGen->Random() < exp(-dE)) { //Accept changes
            box.E_tot += dE;
            box.numAccepts++;
        } else { //Reverse changes
            //Restore Beads
            beads[i] = save;
            if (beads[i].hasSites) {
                update_atom_attachments(i);
            }
            //Restore Grids
            
            box.numRejects++;
        }
        
        
    }
    
}


void Sim::MC_translate_molecule() { //Randomly translates an entire molecule -- however that is defined
  
    int moly,i,j,k,dir,core,disk_length,ref,first,last,tot_move,nAtoms,m;
    double x,y,z,xCenter,yCenter,zCenter;   
    double pre_energy, post_energy, dE;
    LAMMPS_NS::vector dist1; // Contains distance vector for each disk from center of mass for rotation
    LAMMPS_NS::vector dist2; // Contains new distance vector for each disk from center of mass for rotation
    double max_dist = box.translateMax;
    double a,b;
    bool flag,noWallBreak;


    for (int move=0;move<box.numMolecules;move++) {
        
        noWallBreak = true;

        do {
            moly = double(box.numMolecules)*RanGen->Random();
        } while(moly >= box.numMolecules);
        
        first = moleculeList[moly][0];
        last = moleculeList[moly][1];
        nAtoms = last-first+1;

        Bead *save = new Bead[nAtoms];
        Bead *update_chain = new Bead[nAtoms];

        for (i=0;i<nAtoms;i++) {
            save[i] = beads[first+i];
            update_chain[i] = beads[first+i];
            if (beads[first+i].isAnchor) noWallBreak = false;
        }
       
        x = max_dist*(2*RanGen->Random()-1.0); //Calculation for new center of mass
        y = max_dist*(2*RanGen->Random()-1.0); //Calculation for new center of mass
        z = max_dist*(2*RanGen->Random()-1.0); //Calculation for new center of mass
        
        
        for (j=0;j<nAtoms;j++) {
            update_chain[j].rn[0] += x;
            update_chain[j].rn[1] += y;
            update_chain[j].rn[2] += z;
            PBC_shift(update_chain[j].r,update_chain[j].rn);
            update_chain[j].grid = mapToGrid(update_chain[j].r);
            if ((update_chain[j].rn[2] < 0.0 || update_chain[j].rn[2] > box.z) && box.hasWall) {
                noWallBreak = false;
            }
    

        }

        if (noWallBreak) {
            //Save any grids that will be affected by the move
            pre_energy = 0.0;

            for (i=0;i<nAtoms;i++) {
                pre_energy += calc_nonbond_bead_energy(first+i);
                m = beads[first+i].grid;
                grid[m].bead_list.deleteNode(&beads[first+i]);
                beads[first+i] = update_chain[i];
            }

            post_energy = 0.0;

            for (i=0;i<nAtoms;i++) {
                m = mapToGrid(beads[first+i].r);
                grid[m].bead_list.pushNode(&beads[first+i]);
                post_energy += calc_nonbond_bead_energy(first+i);
            }

            dE = post_energy - pre_energy;
        }
        if (RanGen->Random() < exp(-dE) && noWallBreak) { //Accept changes
            box.E_tot += dE; 
            box.numAccepts++;
        } else { 
            //Restore Beads
            for (i=0;i<nAtoms;i++) {
                m = beads[first+i].grid;
                grid[m].bead_list.deleteNode(&beads[first+i]);
                beads[first+i] = save[i];
                m = beads[first+i].grid;
                grid[m].bead_list.pushNode(&beads[first+i]);
            }
            //Restore Grids
            box.numRejects++;
        }


        delete [] save;
        delete [] update_chain;

    }
    

}

void Sim::MC_rotate_molecule() { //Randomly rotate an entire molecule around its center of mass
    int moly,m,i,j,k,dir,core,disk_length,ref,first,last,tot_move,nAtoms;
    double x,y,z,xCenter,yCenter,zCenter,theta;   
    double pre_energy, post_energy, dE;
    LAMMPS_NS::vector dist1; // Contains distance vector for each disk from center of mass for rotation
    LAMMPS_NS::vector dist2; // Contains new distance vector for each disk from center of mass for rotation
    LAMMPS_NS::vector axis; // Contains new distance vector for each disk from center of mass for rotation
    quaternion quat;
    double max_theta = box.rotateMax;
    double a,b;
    bool flag,noWallBreak;


    for (int move=0;move<box.numMolecules;move++) {
        
        noWallBreak = true;

        do {
            moly = double(box.numMolecules)*RanGen->Random();
        } while(moly >= box.numMolecules);
        
        first = moleculeList[moly][0];
        last = moleculeList[moly][1];
        nAtoms = last-first+1;

        Bead *save = new Bead[nAtoms];
        Bead *update_chain = new Bead[nAtoms];

        for (i=0;i<nAtoms;i++) {
            save[i] = beads[first+i];
            update_chain[i] = beads[first+i];

            if(beads[first+i].isAnchor) noWallBreak = false;
        }
       
        xCenter = 0.0;
        yCenter = 0.0;
        zCenter = 0.0;
        for (j=first;j<=last;j++) {
            xCenter += beads[j].rn[0];
            yCenter += beads[j].rn[1];
            zCenter += beads[j].rn[2];
        }
        xCenter /= double(nAtoms); //Calculation of center of mass of polymer
        yCenter /= double(nAtoms);
        zCenter /= double(nAtoms);
        
        theta = max_theta*(2*RanGen->Random()-1);
        calc_random_vector(axis); //Random unit vector to rotate around
        quat[0] = cos(theta/2); 
        for (j=0;j<3;j++) {
            quat[j+1] = axis[j]*sin(theta/2);
        }
        for (i=0;i<nAtoms;i++) {
            dist1[0] = beads[first+i].rn[0] - xCenter;
            dist1[1] = beads[first+i].rn[1] - yCenter;
            dist1[2] = beads[first+i].rn[2] - zCenter;

            quat_vec_rot(dist2,dist1,quat);

            update_chain[i].rn[0] = xCenter + dist2[0];
            update_chain[i].rn[1] = yCenter + dist2[1];
            update_chain[i].rn[2] = zCenter + dist2[2];

            PBC_shift(update_chain[i].r,update_chain[i].rn);

            quat_vec_rot(update_chain[i].u,save[i].u,quat);
            quat_vec_rot(update_chain[i].f,save[i].f,quat);
            quat_vec_rot(update_chain[i].v,save[i].v,quat);

            update_chain[i].grid = mapToGrid(update_chain[i].r);
            if ((update_chain[i].rn[2] < 0.0 || update_chain[i].rn[2] > box.z) && box.hasWall) {
                noWallBreak = false;
            }

        }
        
        if(noWallBreak) {
            //Save any grids that will be affected by the move
            pre_energy = 0.0;

            for (i=0;i<nAtoms;i++) {
                pre_energy += calc_nonbond_bead_energy(first+i);
                m = beads[first+i].grid;
                grid[m].bead_list.deleteNode(&beads[first+i]);
                beads[first+i] = update_chain[i];
            }

            post_energy = 0.0;

            for (i=0;i<nAtoms;i++) {
                m = mapToGrid(beads[first+i].r);
                grid[m].bead_list.pushNode(&beads[first+i]);
                post_energy += calc_nonbond_bead_energy(first+i);
            }

            dE = post_energy - pre_energy;
        }
        if (RanGen->Random() < exp(-dE) && noWallBreak) { //Accept changes
            box.E_tot += dE; 
            box.numAccepts++;
        } else { 
            //Restore Beads
            for (i=0;i<nAtoms;i++) {
                m = beads[first+i].grid;
                grid[m].bead_list.deleteNode(&beads[first+i]);
                beads[first+i] = save[i];
                m = beads[first+i].grid;
                grid[m].bead_list.pushNode(&beads[first+i]);
            }
            //Restore Grids
            box.numRejects++;
        }
    
        delete [] save;
        delete [] update_chain;
    }


}

/*
void Sim::MC_pivot() {

    int moly,i,j,k,dir,core,disk_length,ref,first,last,tot_move,nAtoms;
    double x,y,z,xCenter,yCenter,zCenter,theta;   
    double pre_energy, post_energy, dE;
    LAMMPS_NS::vector dist1; // Contains distance vector for each disk from center of mass for rotation
    LAMMPS_NS::vector dist2; // Contains new distance vector for each disk from center of mass for rotation
    LAMMPS_NS::vector axis; // Contains new distance vector for each disk from center of mass for rotation
    quaternion quat;
    double max_theta = box.pivotMax;
    std::vector<int> grid_ids; //ID's of grids affected by the monte carlo move
    double a,b;
    bool flag,noWallBreak;


    for (int move=0;move<box.numMolecules;move++) {
        
        noWallBreak = true;

        do {
            moly = double(box.numMolecules)*RanGen->Random();
        } while(moly >= box.numMolecules);
        
        first = moleculeList[moly][0];
        last = moleculeList[moly][1];
        nAtoms = last-first+1;
        ref = first;

        do {
            core = double(nAtoms)*RanGen->Random();
        } while(core >= nAtoms);
        
        core += ref; //id of atom around which rotation will occur

        if (RanGen->Random() > 0.5) {
            dir = 1;
            first = core;
        } else {
            dir = -1;
            last = core;
        }

        nAtoms = last-first+1;

        Bead *save = new Bead[nAtoms];
        Bead *update_chain = new Bead[nAtoms];
        Grid *save_grid = new Grid[2*nAtoms]; //Maximum number of grids affected are nAtoms for initial grids, and nAtoms for final grids

        for (i=0;i<nAtoms;i++) {
            save[i] = beads[first+i];
            update_chain[i] = beads[first+i];
            flag = true;
            for (j=0;j<grid_ids.size();j++) {
                if (grid_ids[j] == save[i].grid) flag = false;
            }
            if (flag) grid_ids.push_back(save[i].grid);

            if (beads[first+i].isAnchor && first+i != core) noWallBreak = false;
        }
       
        theta = max_theta*(2*RanGen->Random()-1);
        calc_random_vector(axis); //Random unit vector to rotate around
        quat[0] = cos(theta/2); 
        for (j=0;j<3;j++) {
            quat[j+1] = axis[j]*sin(theta/2);
        }
        
        for (i=0;i<nAtoms;i++) {
            dist1[0] = beads[first+i].rn[0] - beads[core].rn[0];
            dist1[1] = beads[first+i].rn[1] - beads[core].rn[1];
            dist1[2] = beads[first+i].rn[2] - beads[core].rn[2];

            quat_vec_rot(dist2,dist1,quat);

            update_chain[i].rn[0] = beads[core].rn[0] + dist2[0];
            update_chain[i].rn[1] = beads[core].rn[1] + dist2[1];
            update_chain[i].rn[2] = beads[core].rn[2] + dist2[2];

            PBC_shift(update_chain[i].r,update_chain[i].rn);

            quat_vec_rot(update_chain[i].u,save[i].u,quat);

            update_chain[i].grid = mapToGrid(update_chain[i].r);
            if ((update_chain[i].rn[2] < 0.0 || update_chain[i].rn[2] > box.z) && box.hasWall) {
                noWallBreak = false;
            }
            flag = true;
            for (j=0;j<grid_ids.size();j++) {
                if (grid_ids[j] == update_chain[i].grid) flag = false;
            }
            if (flag) grid_ids.push_back(update_chain[i].grid);

        }
        

        
        if(noWallBreak) {
            //Save any grids that will be affected by the move
            pre_energy = calc_bead_energy(core); //Only bead that has angles change is the one acting as the center of rotation
            for (i=0;i<grid_ids.size();i++) {
                save_grid[i] = grid[grid_ids[i]];
                pre_energy += calc_grid_energy(grid_ids[i]);
            }

            for (i=0;i<nAtoms;i++) {
                removeBeadFromGrid(first+i);
                beads[first+i] = update_chain[i];
            }

            update_single_u_vector(core);

            for (i=0;i<nAtoms;i++) {
                addBeadToGrid(first+i);
            }

            post_energy = calc_bead_energy(core);
            for (i=0;i<grid_ids.size();i++) {
                post_energy += calc_grid_energy(grid_ids[i]);
            }

            dE = post_energy - pre_energy;
        }
        if (RanGen->Random() < exp(-dE) && noWallBreak) { //Accept changes
            box.E_tot += dE; 
            box.numAccepts++;
        } else { 
            //Restore Beads
            for (i=0;i<nAtoms;i++) {
                beads[first+i] = save[i];
            }
            //Restore Grids
            if(noWallBreak) {
                for (i=0;i<grid_ids.size();i++) {
                    grid[grid_ids[i]] = save_grid[i];
                }
            }
            box.numRejects++;
        }
    
        delete [] save;
        delete [] update_chain;
        delete [] save_grid;
        grid_ids.clear();
    }


}
*/

void Sim::MC_reptate() { //Performs a reptation move on the polymer. Does not change any bond lengths, so it can be used for rigid polymers
   /* 
    int moly,i,j,k,dir,core,disk_length,ref,first,last,tot_move,nAtoms;
    double x,y,z,xCenter,yCenter,zCenter;   
    double pre_energy, post_energy, dE;
    LAMMPS_NS::vector dist1; // Contains distance vector for each disk from center of mass for rotation
    LAMMPS_NS::vector axis; // Contains new distance vector for each disk from center of mass for rotation
    double max_dist = box.translateMax;
    std::vector<int> grid_ids; //ID's of grids affected by the monte carlo move
    double a,b,bond_length;
    bool flag,noWallBreak;


    for (int move=0;move<box.numMolecules;move++) {
        
        noWallBreak = true;

        do {
            moly = double(box.numMolecules)*RanGen->Random();
        } while(moly >= box.numMolecules);
        
        first = moleculeList[moly][0];
        last = moleculeList[moly][1];
        nAtoms = last-first+1;

        Bead *save = new Bead[nAtoms];
        Bead *update_chain = new Bead[nAtoms];
        Grid *save_grid = new Grid[nAtoms+1]; //There may be one additional grid point included due to the reptation. All others affected will aleady be there

        for (i=0;i<nAtoms;i++) {
            save[i] = beads[first+i];
            update_chain[i] = beads[first+i];
            flag = true;
            for (j=0;j<grid_ids.size();j++) {
                if (grid_ids[j] == save[i].grid) flag = false;
            }
            if (flag) grid_ids.push_back(save[i].grid);

            if (beads[first+i].isAnchor) noWallBreak = false;
        }
      
        //Need to calculate the correct length of the bond that is being moved from one end to another
        

        if (RanGen->Random() < 0.5) { //New bead is attached onto last bead. First bead in chain is removed
            dir = 1;
            bond_length = 0.0;
            bond_length += (beads[first].xn - beads[first+1].xn) * (beads[first].xn - beads[first+1].xn);
            bond_length += (beads[first].yn - beads[first+1].yn) * (beads[first].yn - beads[first+1].yn);
            bond_length += (beads[first].zn - beads[first+1].zn) * (beads[first].zn - beads[first+1].zn);
        } else {                      //New bead is attached onto first bead. Last bead in chain is removed
            dir = -1;
            bond_length = 0.0;
            bond_length += (beads[last].xn - beads[last-1].xn) * (beads[last].xn - beads[last-1].xn);
            bond_length += (beads[last].yn - beads[last-1].yn) * (beads[last].yn - beads[last-1].yn);
            bond_length += (beads[last].zn - beads[last-1].zn) * (beads[last].zn - beads[last-1].zn);
        }

        //Generate a random vector that is bond_length long
        

        calc_random_vector(axis);

        for (i=0;i<3;i++) {
            axis[i] = bond_length*axis[i];
        }


        x = max_dist*(2*RanGen->Random()-1.0); //Calculation for new center of mass
        y = max_dist*(2*RanGen->Random()-1.0); //Calculation for new center of mass
        z = max_dist*(2*RanGen->Random()-1.0); //Calculation for new center of mass
       

        
        for (j=0;j<nAtoms;j++) {
            update_chain[j].rn[0] += x;
            update_chain[j].rn[1] += y;
            update_chain[j].rn[2] += z;
            PBC_shift(update_chain[j].r,update_chain[j].rn);
            update_chain[j].grid = mapToGrid(update_chain[j].r);
            if ((update_chain[j].rn[2] < 0.0 || update_chain[j].rn[2] > box.z) && box.hasWall) {
                noWallBreak = false;
            }
            flag = true;
            for (i=0;i<grid_ids.size();i++) {
                if (grid_ids[i] == update_chain[j].grid) flag = false;
            }
            if (flag) grid_ids.push_back(update_chain[j].grid);
    

        }

        if (noWallBreak) {
            //Save any grids that will be affected by the move
            pre_energy = 0.0;
            for (i=0;i<grid_ids.size();i++) {
                save_grid[i] = grid[grid_ids[i]];
                pre_energy += calc_grid_energy(grid_ids[i]);
            }

            for (i=0;i<nAtoms;i++) {
                removeBeadFromGrid(first+i);
                beads[first+i] = update_chain[i];
                addBeadToGrid(first+i);
            }

            post_energy = 0.0;
            for (i=0;i<grid_ids.size();i++) {
                post_energy += calc_grid_energy(grid_ids[i]);
            }

            dE = post_energy - pre_energy;
        }
        if (RanGen->Random() < exp(-dE) && noWallBreak) { //Accept changes
            box.E_tot += dE; 
            box.numAccepts++;
        } else { 
            //Restore Beads
            for (i=0;i<nAtoms;i++) {
                beads[first+i] = save[i];
            }
            //Restore Grids
            for (i=0;i<grid_ids.size();i++) {
                grid[grid_ids[i]] = save_grid[i];
            }
            box.numRejects++;
        }
    
        delete [] save;
        delete [] update_chain;
        delete [] save_grid;

        grid_ids.clear();
    }
    
   
*/
}
/**************************************************************************************************************/


/**********************************DATA OUTPUT AND VISUALIZATION***********************************************/
void Sim::shiftCOM() { //Shifts the center of mass of all the nonperiodic coordinates back into the simulation box 
    int i,j,m,start,end,nAtoms,flag;
    double xCenter,yCenter,zCenter,dx,dy,dz;
    for (i=0;i<box.numMolecules;i++) {
        flag = 0;
        xCenter = 0.0;
        yCenter = 0.0;
        zCenter = 0.0;
        start = moleculeList[i][0];
        end = moleculeList[i][1];
        nAtoms = end-start+1;
        for (j=start;j<=end;j++) {
            if (beads[j].isAnchor) {
                flag = 1;
            }
            xCenter += beads[j].rn[0];
            yCenter += beads[j].rn[1];
            zCenter += beads[j].rn[2];
        }
        if (flag == 1) {
            continue;
        }
        xCenter /= double(nAtoms);
        yCenter /= double(nAtoms);
        zCenter /= double(nAtoms);
        dx = box.x*floor(xCenter/box.x);
        dy = box.y*floor(yCenter/box.y);
        dz = box.z*floor(zCenter/box.z);
        for (j=start;j<=end;j++) {
            beads[j].rn[0] -= dx;
            beads[j].rn[1] -= dy;
            beads[j].rn[2] -= dz;
        }
    }
}

void Sim::printXYZ() {
    //Print all the molecules then print out rings to visualize disk-like "atoms"
    
    int i,j;
    FILE *dump;
    dump = fopen("config.xyz","a");

    fprintf(dump,"%d\nConfiguratoin\n",box.numBeads);
    for (i=0;i<box.numBeads;i++) {
        fprintf(dump,"0\t%f\t%f\t%f\n",beads[i].rn[0],beads[i].rn[1],beads[i].rn[2]);
    }
    fclose(dump);
}

void Sim::printPSF() { //Print out the PSF file. Only visualizing bonds here, not twist or angle potentials
      FILE *sout;
      int counter,type,pair,chain_switch,bead_type;
      char t[2]="O";
      sout=fopen("polymer.psf","w");
      fprintf(sout,"*\n*\n*\n*\n*\n\n");
      fprintf(sout,"%7d !NATOMS\n",box.numBeads);
      for(counter=1;counter<=box.numBeads;counter++){
        if(counter<=box.numBeads){

        if (beads[counter-1].type == 0) t[0] = 'O'; 
        if (beads[counter-1].type == 1) t[0] = 'N'; 
        if (beads[counter-1].type == 2) t[0] = 'C'; 
        if (beads[counter-1].type == 3) t[0] = 'C'; 

          bead_type=beads[counter-1].type;

          fprintf(sout, "%8d ", counter);
          fprintf(sout, "POLY ");
          fprintf(sout, "%-4d ",1);
          fprintf(sout, "%s  ", "POL");
          fprintf(sout, "%-5s ",t);
          fprintf(sout, "%3d  ", bead_type);
          fprintf(sout, "%13.6e   ",0.0);
          fprintf(sout, "%7.3lf           0\n", 1.0);
        }

        else{
          t[0]='S'; bead_type=2;
          fprintf(sout, "%8d ", counter);
          fprintf(sout, "SOLV ");
          fprintf(sout, "%-4d ",1);
          fprintf(sout, "%s  ", "SOL");
          fprintf(sout, "%-5s ",t);
          fprintf(sout, "%3d  ", bead_type);
          fprintf(sout, "%13.6e   ",0.0);
          fprintf(sout, "%7.3lf           0\n", 1.0);

        }        
        //fprintf(sout,"%8d\tU\t%4d\tPOL\t%1s\t%1s\t0.0000000\t1.0000\t0\n  ",counter,counter-1,t,t);
      }


      fprintf(sout,"\n%8d !NBOND: bond\n",box.numBonds+box.numAttachments);




      pair=1;
      counter=1;
      int i;
      for(i=0;i<box.numBonds;i++){

        if(pair==5)
          {pair=1;fprintf(sout,"\n");}
        //printf("**counter:\t%d\n",counter);
        fprintf(sout,"%8d%8d",bonds[i].id1+1,bonds[i].id2+1);
        pair++;
      }
      for (i=0;i<box.numAttachments;i++) {
        
        if(pair==5)
          {pair=1;fprintf(sout,"\n");}
        //printf("**counter:\t%d\n",counter);
        fprintf(sout,"%8d%8d",attachments[i].id_core+1,attachments[i].id_site+1);
        pair++;
      }

      fprintf(sout,"\n%8d !NTHETA: angles\n",0);
      fprintf(sout,"\n%8d !NPHI: dihedrals\n",0);
      fprintf(sout,"\n%8d !NIMPR\n",0);
      fprintf(sout,"\n%8d !HDON\n",0);
      fprintf(sout,"\n%8d !HACC\n",0);
      fprintf(sout,"\n%8d !NNB\n",0);


      fclose(sout);

}


void Sim::empty_containers() {
    bondList.clear();
    angleList.clear();
    moleculeList.clear();
}

/*(
void Sim::printPOV() {
    FILE *pov = fopen("disks.xyz","w");
    fprintf(pov, "ITEM: TURN\n%d\nITEM: NUMBER OF ATOMS\n%d\nITEM: BOX BOUNDS\n-10 10\n-10 10\n-10 10\nITEM: ATOMS\n",box.cycle,box.numTotal);

    double scalex = 10.0/(box.x/2);
    double scaley = 10.0/(box.y/2);
    double scalez = 10.0/(box.z/2);

    for (int i=0;i<box.numTotal;i++) {
        fprintf(pov,"%d 1 %f %f %f %f %f %f\n", i, (disk[i].r[0]-box.x/2)*scalex, (disk[i].r[1]-box.y/2)*scaley, (disk[i].r[2]-box.z/2)*scalez, disk[i].f[0], disk[i].f[1], disk[i].f[2]); 
    }
    fclose(pov);
}
*/

void Sim::writeEnergy() {
    
    int i,j,k;
    int t1,t2;
    double bond, angle,twist, compression, chi, nematic;
    double phi_tot, chi_single;
    bond = 0.0; compression = 0.0; nematic = 0.0; chi = 0.0; angle = 0.0; twist = 0.0;
    for (i=0;i<box.numBonds;i++) {
        bond += calc_bond_energy(i); 
    }
    for (i=0;i<box.numAngles;i++) {
        angle += calc_angle_energy(i);
    }
    for (i=0;i<box.numTwists;i++) {
        twist += calc_twist_energy(i);
    }
    chi = calc_total_nonbond_energy();
    
    FILE *en = fopen("Energy.txt","a");
    fprintf(en,"%d\t%f\t%f\t%f\t%f\t%f\t%f\n",box.cycle,bond,angle,twist,chi,bond+angle+twist+chi,box.E_tot);
    fclose(en);
}

void Sim::dumpData() {
    FILE *out = fopen("data.txt","a");
    int i;
    fprintf(out,"%d\n",box.cycle);
    for (i=0;i<box.numBeads;i++) {
        fprintf(out,"%d\t", beads[i].id);
        fprintf(out,"%lf\t%lf\t%lf\t", beads[i].r[0], beads[i].r[1], beads[i].r[2]);
        fprintf(out,"%lf\t%lf\t%lf\t", beads[i].rn[0], beads[i].rn[1], beads[i].rn[2]);
        fprintf(out,"%lf\t%lf\t%lf\t", beads[i].u[0], beads[i].u[1], beads[i].u[2]);
        fprintf(out,"%lf\t%lf\t%lf\t", beads[i].f[0], beads[i].f[1], beads[i].f[2]);
        fprintf(out,"%lf\t%lf\t%lf\t", beads[i].v[0], beads[i].v[1], beads[i].v[2]);
        fprintf(out,"%lf\t", beads[i].vol);
        fprintf(out,"%d\t", beads[i].type);
        fprintf(out,"%d\t", beads[i].isAnchor);
        fprintf(out,"%d\t", beads[i].isRod);
        fprintf(out,"%d\t", beads[i].isSite);
        fprintf(out,"%d\n", beads[i].hasSites);

    }
    fclose(out);
}

void Sim::dumpGrid() {
    /*
    
    FILE* out = fopen("grid.txt","a");
    int i,t;
    fprintf(out,"%d\n",box.cycle);
    for (i=0;i<box.numGrid;i++) {
        for (t=0;t<box.numTypes;t++) {
            fprintf(out,"%f\t",grid[i].phi[t]);
        }
        fprintf(out,"\n");
    }
    fclose(out);
    */
}

/**************************************************************************************************************/


/**************************************************************************************************************/
