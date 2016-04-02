/* 
 * File:   main.cpp
 * Author: Sheridan Beckwith Green (sbg@unc.edu)
 * Part of the UNCC FAST Library
 *
 * Created on May 21, 2015, 10:01 AM
 */

#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
#include <climits>
#include <algorithm>
#include <iterator>
#include <cmath>
#include <float.h>
#include <assert.h>
#include <time.h>
#include <unistd.h>
#include <bmpg_uncc_edu/fast/FASTInitializer.hpp>
#include <bmpg_uncc_edu/fast/ParameterFile.hpp>
#include <bmpg_uncc_edu/chemistry/PDBProtein.hpp>
#include <bmpg_uncc_edu/chemistry/helper/PDBProteinPreprocessor.hpp>
#include "bmpg_uncc_edu/chemistry/PDBAtom.hpp"
#include "bmpg_uncc_edu/forcefield/OPLSForceField.hpp"
#include "bmpg_uncc_edu/forcefield/ForceField.hpp"
#include "bmpg_uncc_edu/chemistry/Element.hpp"
#include <bmpg_uncc_edu/util/logger/LoggerFactory.hpp>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>

using namespace std;
using namespace bmpg_uncc_edu::util;
using namespace bmpg_uncc_edu::chemistry;
using namespace bmpg_uncc_edu::chemistry::hbond;
using namespace bmpg_uncc_edu::fast;
using namespace bmpg_uncc_edu::forcefield;
using namespace bmpg_uncc_edu::chemistry::helper;
using namespace bmpg_uncc_edu::util::logger;

// Macros
#define SQR(x)      ((x)*(x))                        // x^2 

//default constants used throughout
float maxR_ss = 1.8; //this will change depending on input set of vdw radii, largest one in set
float hashSpacing = 3.0;
float R_probe = 1.4;//1.4;//3.15;//1.2;// //angstroms, used for glueing criterion (water TRY 3.15 A)
float R_bw = R_probe * 4.0;//6.0; //angstroms of boundary water
float a_grid = 0.5; //angstroms (grid spacing)
float solvGap = ceil(a_grid+maxR_ss+R_bw); //desired distance from axis to x,y,z min atom coordinate -> should be greater than a_grid + MAXIMUM Rss + Rbw
int * clusterLabel; int * clusterSize; int n_labels; int maxCluster = 1;
int atoms; float ** atomCoords;
//hashmaps used throughout
int *** hashGrid; int * gridChain; int L; //dimensionality of atom hashgrid
float angles[3];
string pdb_file_name;

//for final data output
ofstream voidHistFile;
ofstream microvoidHistFile;
ofstream summaryFile;
ofstream largestMVFile;
const int global_rotation_max = 3;
int global_rot_iter; //needs to be global so we can store in arrays
float glob_prot_vol_gdpts[global_rotation_max]={0}; float glob_mv_gdpts[global_rotation_max]={0}; float glob_cav_gdpts[global_rotation_max]={0}; float glob_solv_gdpts[global_rotation_max]={0};
float glob_mv_clusters[global_rotation_max]={0}; float glob_cav_clusters[global_rotation_max]={0}; float glob_cav_cluster_errors[global_rotation_max]={0};
float glob_largest_mv_gdpts[global_rotation_max]={0}; float glob_largest_cav_gdpts[global_rotation_max]={0}; float glob_gdpts_total[global_rotation_max]={0}; 
float glob_prot_vol[global_rotation_max]={0}; float glob_mv_vol[global_rotation_max]={0}; float glob_cav_vol[global_rotation_max]={0}; float glob_solv_vol[global_rotation_max]={0};
float glob_total_vol[global_rotation_max]={0}; float glob_packing_density[global_rotation_max]={0}; float glob_cpu_time[global_rotation_max]={0};


//variables for geometric method
int max_free = 500000; //this parameter can be changed
int *** free_hashGrid; int ** free_gridChain; //dimensionality of free volume hashgrid
int free_x_dims;
int free_y_dims;
int free_z_dims;

int count1=0; int count2=0; int count3=0; int count4=0;

float old_vdw_radii[100] = {0}; float original_vdw_radii[100] = {0}; float new_vdw_radii[100] = {0}; int vdwTypes = 0;

float R_check = (R_probe + maxR_ss) + 0.51*a_grid; //for checking in geometrical method
float R_cap = R_probe + 1.7321*a_grid;

float a_gridMax = 2.0; //for data analysis
float hydrogenVDW = 1.200; // for inertia rotate data analysis

//spring model parameters
int maxIter = 1000;//1000; //just to ensure that our spring models eventually converge
float E_zero = 0.0002;//0.0001*(R_probe+maxR_ss)*(R_probe+maxR_ss);//0.0003;//0.00025;//0.00024;//0.0001*(R_probe+maxR_ss)*(R_probe+maxR_ss);
float F_zero = 0.008;//0.01*(R_probe + maxR_ss);
float viscosity = 1.25;//0.75;
float max_dr_mag = 0.4;//0.8 * a_grid; //what should this be? experiment with it
float min_dr_mag = 0.0005;//a_grid / 1000.0;
float targetR_probe = R_probe * 1.10;

int largestVoid = 0; int largestMicrovoid = 0;
//only want to make these guys two arrays if we are actually running cluster statistics
//FIXME
int voidHist[120000000] = {0};
int microvoidHist[120000000] = {0};

int statRunNum;
float grid_factor = 0.51;
string connectivity_type;
int connect_num;

float evector[3][3]={0}; float evalue[3]={0}; //for external matrix diagonalization algorithm

//needed global in order to shift when making visualization PDBs
float xshift, yshift, zshift;
PDBAtom * voidAtom = new PDBAtom; PDBAtom * microvoidAtom = new PDBAtom; int voidAtomNum = 0; int microvoidAtomNum = 0;
ofstream timeFile;
ofstream voidVis;
ofstream microvoidVis;
ofstream variableFile;

//for connectivity comparison
float global_microvoid; float global_mv_clusters;

//HEAP MEMORY SO FAR:
//atomCoords (only gets initialized once), hashGrid and gridChain (must use destroyHash()), 


// ----------------------------------------------------------------------------
int dsyevj3(float A[3][3], float Q[3][3], float w[3])
// ----------------------------------------------------------------------------
// Calculates the eigenvalues and normalized eigenvectors of a symmetric 3x3
// matrix A using the Jacobi algorithm.
// The upper triangular part of A is destroyed during the calculation,
// the diagonal elements are read but not destroyed, and the lower
// triangular elements are not referenced at all.
// ----------------------------------------------------------------------------
// Parameters:
//   A: The symmetric input matrix
//   Q: Storage buffer for eigenvectors
//   w: Storage buffer for eigenvalues
// ----------------------------------------------------------------------------
// Return value:
//   0: Success
//  -1: Error (no convergence)
// ----------------------------------------------------------------------------
{
  const int n = 3;
  float sd, so;                  // Sums of diagonal resp. off-diagonal elements
  float s, c, t;                 // sin(phi), cos(phi), tan(phi) and temporary storage
  float g, h, z, theta;          // More temporary storage
  float thresh;
  
  // Initialize Q to the identitity matrix
#ifndef EVALS_ONLY
  for (int i=0; i < n; i++)
  {
    Q[i][i] = 1.0;
    for (int j=0; j < i; j++)
      Q[i][j] = Q[j][i] = 0.0;
  }
#endif

  // Initialize w to diag(A)
  for (int i=0; i < n; i++)
    w[i] = A[i][i];

  // Calculate SQR(tr(A))  
  sd = 0.0;
  for (int i=0; i < n; i++)
    sd += fabs(w[i]);
  sd = SQR(sd);
 
  // Main iteration loop
  for (int nIter=0; nIter < 50; nIter++)
  {
    // Test for convergence 
    so = 0.0;
    for (int p=0; p < n; p++)
      for (int q=p+1; q < n; q++)
        so += fabs(A[p][q]);
    if (so == 0.0)
      return 0;

    if (nIter < 4)
      thresh = 0.2 * so / SQR(n);
    else
      thresh = 0.0;

    // Do sweep
    for (int p=0; p < n; p++)
      for (int q=p+1; q < n; q++)
      {
        g = 100.0 * fabs(A[p][q]);
        if (nIter > 4  &&  fabs(w[p]) + g == fabs(w[p])
                       &&  fabs(w[q]) + g == fabs(w[q]))
        {
          A[p][q] = 0.0;
        }
        else if (fabs(A[p][q]) > thresh)
        {
          // Calculate Jacobi transformation
          h = w[q] - w[p];
          if (fabs(h) + g == fabs(h))
          {
            t = A[p][q] / h;
          }
          else
          {
            theta = 0.5 * h / A[p][q];
            if (theta < 0.0)
              t = -1.0 / (sqrt(1.0 + SQR(theta)) - theta);
            else
              t = 1.0 / (sqrt(1.0 + SQR(theta)) + theta);
          }
          c = 1.0/sqrt(1.0 + SQR(t));
          s = t * c;
          z = t * A[p][q];

          // Apply Jacobi transformation
          A[p][q] = 0.0;
          w[p] -= z;
          w[q] += z;
          for (int r=0; r < p; r++)
          {
            t = A[r][p];
            A[r][p] = c*t - s*A[r][q];
            A[r][q] = s*t + c*A[r][q];
          }
          for (int r=p+1; r < q; r++)
          {
            t = A[p][r];
            A[p][r] = c*t - s*A[r][q];
            A[r][q] = s*t + c*A[r][q];
          }
          for (int r=q+1; r < n; r++)
          {
            t = A[p][r];
            A[p][r] = c*t - s*A[q][r];
            A[q][r] = s*t + c*A[q][r];
          }

          // Update eigenvectors
#ifndef EVALS_ONLY          
          for (int r=0; r < n; r++)
          {
            t = Q[r][p];
            Q[r][p] = c*t - s*Q[r][q];
            Q[r][q] = s*t + c*Q[r][q];
          }
#endif
        }
      }
  }

  return -1;
}

void userInput(int argc, char** argv){
    //bring in parameters from command line
    int c;
    while ((c = getopt (argc, argv, "n:a:p:z:c:")) != -1)
    switch (c){
        case 'n':
            pdb_file_name = optarg;
            break;
        case 'a':
            a_grid = atof(optarg);
            break;
        case 'p':
            R_probe = atof(optarg);
            break;
        case 'z':
            statRunNum = atof(optarg); //for differentiating our statistical files
            break;
        case 'c':
            connectivity_type = optarg;
            if (optarg == "H"){
                connect_num = 7;
            }
            else if(optarg == "L"){
                connect_num = 3;
            }
            else{
                connect_num = 7; //default to High connectivity
            }
            break;
        default:
            abort();
    }
}


//requires that PDB files be minimized and protonated
int atomCount(PDBProtein *protein) {
    return protein->num_atoms();
}

float rotatex(float phi, float theta, float psi, float x, float y, float z){
    //return cos(psi)*((cos(theta)*x)-(sin(theta)*y))+(sin(psi)*z);
    return cos(psi)*cos(theta)*x + (cos(psi)*sin(theta)*sin(phi)- cos(phi)*sin(psi))*y + (sin(psi)*sin(phi) + cos(psi)*cos(phi)*sin(theta))*z;
}

float rotatey(float phi, float theta, float psi, float x, float y, float z){
    //return cos(phi)*((sin(theta)*x)+(cos(theta)*y))-sin(phi)*((cos(psi)*z)-(sin(psi)*((cos(theta)*x)-(sin(theta)*y))));
    return cos(theta)*sin(psi)*x + (cos(psi)*cos(phi) + sin(psi)*sin(theta)*sin(phi))*y + (cos(phi)*sin(psi)*sin(theta) - cos(psi)*sin(phi))*z;
}

float rotatez(float phi, float theta, float psi, float x, float y, float z){
    //return sin(phi)*((sin(theta)*x)+(cos(theta)*y))+cos(phi)*((cos(psi)*z)-(sin(psi)*((cos(theta)*x)-(sin(theta)*y))));
    return -sin(theta)*x + cos(theta)*sin(phi)*y + cos(theta)*cos(phi)*z;
}

void addVDW(float newRad){
    for(int i = 0; i<vdwTypes; i++){
        if(newRad == original_vdw_radii[i]){
            return;
        }
    }
    original_vdw_radii[vdwTypes] = newRad;
    vdwTypes++;
    
    //later write another thing that takes these and maps them to new radii
}

void firstQuadCoordShift(){
    //initialize min and max values for delta, essentially pos. infinity for min, neg. infinity for max. -> guaranteed a change of value
    float xmin = FLT_MAX; float xmax = -FLT_MAX;  
    float ymin = FLT_MAX; float ymax = -FLT_MAX;
    float zmin = FLT_MAX; float zmax = -FLT_MAX;
    
    for(int i=1;i<=atoms;i++){
        //finding max and min values in each dimension
        if (atomCoords[i][1] > xmax) { xmax = atomCoords[i][1];}
        if (atomCoords[i][1] < xmin) { xmin = atomCoords[i][1];}
        if (atomCoords[i][2] > ymax) { ymax = atomCoords[i][2];}
        if (atomCoords[i][2] < ymin) { ymin = atomCoords[i][2];}
        if (atomCoords[i][3] > zmax) { zmax = atomCoords[i][3];}
        if (atomCoords[i][3] < zmin) { zmin = atomCoords[i][3];}
    }
    //storing delta x,y,z in row 1
    float dx, dy, dz;
    atomCoords[0][1] = xmax - xmin; dx = xmax - xmin;
    atomCoords[0][2] = ymax - ymin; dy = ymax - ymin;
    atomCoords[0][3] = zmax - zmin; dz = zmax - zmin;
    
    //implement coordinate shift here
    
    if (xmin > 0) {
        xshift = ceil(xmin/a_gridMax)*a_gridMax;
    } else if(xmin < 0) {
        xshift = floor(xmin/a_gridMax)*a_gridMax;
    }
    if (ymin > 0) {
        yshift = ceil(ymin/a_gridMax)*a_gridMax;
    } else if(ymin < 0) {
        yshift = floor(ymin/a_gridMax)*a_gridMax;
    }
    if (zmin > 0) {
        zshift = ceil(zmin/a_gridMax)*a_gridMax;
    } else if(zmin < 0) {
        zshift = floor(zmin/a_gridMax)*a_gridMax;
    }
    //since adding (MAX-R_ss+a_grid) worth of space on each side, we'll shift such that min atom is at around 3A from each axis (will change)
    //TEST: see what effects different distances from the wall has on the volumes
    for(int i = 1; i<=atoms; i++){
        atomCoords[i][1] = atomCoords[i][1] - (xshift - solvGap);
        atomCoords[i][2] = atomCoords[i][2] - (yshift - solvGap);
        atomCoords[i][3] = atomCoords[i][3] - (zshift - solvGap);
    }  
}

//rotates coordinates and then calls firstQuadCoordShift() to shift back to first quadrant
void rotateCoords(){
    //coordinate rotation
    float tempx; float tempy; float tempz;
    srand(time(NULL));
    angles[0] = 2*M_PI*((float)rand()/RAND_MAX);
    angles[1] = asin(((float)rand()/RAND_MAX));
    angles[2] = M_PI*(2*((float)rand()/RAND_MAX)-1);
    //angles[0] = ((float)rand()/RAND_MAX)*90;
    //angles[1] = ((float)rand()/RAND_MAX)*90;
    //angles[2] = ((float)rand()/RAND_MAX)*90;

    cout << "Angles: " << setprecision(3) << angles[0] << ", " << angles[1] << ", " << angles[2] << "\n";
    
    for(int i=1; i <= atoms; i++){
        tempx = rotatex(angles[0], angles[1], angles[2], atomCoords[i][1], atomCoords[i][2], atomCoords[i][3]);
        tempy = rotatey(angles[0], angles[1], angles[2], atomCoords[i][1], atomCoords[i][2], atomCoords[i][3]);
        tempz = rotatez(angles[0], angles[1], angles[2], atomCoords[i][1], atomCoords[i][2], atomCoords[i][3]);
        atomCoords[i][1] = tempx; atomCoords[i][2] = tempy; atomCoords[i][3] = tempz;
    }
    
    firstQuadCoordShift();
}


//creates 2d array, row 0 holds delta x,y,z values, rows 1 thru atoms holds xyz coordinates of atoms
void getAtomCoords(string pdb_file_name) {
    
    Logger * logger = LoggerFactory::default_logger();
    logger->set_log_level(Logger::CRITICAL);
    
    //loading in parameters and libraries from config file
    ParameterFile* pfile = ParameterFile::instance();
    pfile->load("config.txt");    

    FASTInitializer fi;
    fi.load_user_defined_libraries("config.txt");

    //loading in protein from pdb
    ifstream fin(pdb_file_name.c_str());
    PDBProtein *protein = 0;
    protein = new PDBProtein;
    protein->read(fin);
    fin.close();
    PDBProtein::AtomIterator iter;
    
    //IMPORTANT: compare our OPLS radii to Bondi radii in paper
    //find out if it needs to be preprocessed
    PDBProteinPreprocessor *processor = new PDBProteinPreprocessor();
    processor->process(*protein); //important to get working eventually

    
    //find number of atoms
    atoms = atomCount(protein);
    
    
    //make pointer to atomCoords 2d array:x,y,z coordinates
    atomCoords = (float**) new float* [atoms+1]; //0 to atoms
    for(int i=0;i<=atoms;i++)
      {
        atomCoords[i] = (float*) new float [4]; //0 to 3-> Rss, x,y,z
      }
     //remember to delete this memory after use when running program over many atoms

    OPLSForceField* ff = dynamic_cast<OPLSForceField*>(ForceFieldFactory::default_force_field());
        
    int count = 1;
    const PDBAtom * a;
    OPLS_coul_vdw_record *vdw;
    for (iter = protein->atom_begin(); iter != protein->atom_end(); ++iter) {
            a = *iter;
            atomCoords[count][1] = a->x;
            atomCoords[count][2] = a->y;
            atomCoords[count][3] = a->z;
            //symbol[count] = a->atomic_symbol; //likely no need for this but keep around in case, requires a string array called symbol
            vdw = ff->lookup_coul_vdw_data(a->res, a->atom_name);
            atomCoords[count][0] = vdw->R; //this is Rss.
            addVDW(atomCoords[count][0]); //for doing our transformation
            
            count++;
    }       
}

void generateFreeHash(int* n){
    //design grid and grid_chain for hash map
    //need to figure out dimensionality
    free_x_dims = ceil(n[1]*a_grid/R_probe);
    free_y_dims = ceil(n[2]*a_grid/R_probe);
    free_z_dims = ceil(n[3]*a_grid/R_probe);
    free_hashGrid = (int ***) new int** [free_x_dims]; //initialized to zero, which we need
    for(int i=0; i<free_x_dims; i++){
        free_hashGrid[i] = (int **) new int* [free_y_dims];
        for(int j=0; j<free_y_dims; j++){
            free_hashGrid[i][j] = (int *) new int [free_z_dims];
            for(int k=0; k<free_z_dims; k++){
                free_hashGrid[i][j][k] = 0;
            }
        }
    }
    
    free_gridChain = (int **) new int* [max_free]; //not using index 0, starting at atom 1, according to atomCoords array
    for(int i=0; i<max_free; i++){
        free_gridChain[i] = (int *) new int [4];
        for(int j=0; j<4; j++){
            free_gridChain[i][j] = 0;
        }
    }
}

void generateHash(){
    //design grid and grid_chain for hash map
    L = ceil(cbrt(2*atoms)); //assuming we want a 0.5 average density of cells in hash grid
    hashGrid = (int ***) new int** [L]; //initialized to zero, which we need
    for(int i=0; i<L; i++){
        hashGrid[i] = (int **) new int* [L];
        for(int j=0; j<L; j++){
            hashGrid[i][j] = (int *) new int [L];
            for(int k=0; k<L; k++){
                hashGrid[i][j][k] = 0;
            }
        }
    }
    gridChain = new int[atoms+1]; //not using index 0, starting at atom 1, according to atomCoords array
    for(int iter = 0; iter <=atoms; iter++){
        gridChain[iter] = 0;
    }
    
    //construct grid and grid_chain
    int i; int j; int k; int chainLoc; //locations to be placed in grid
    for(int iter=1; iter<=atoms; iter++){
        //determine hash grid location
        i = (int)floor(atomCoords[iter][1] / hashSpacing) % L;
        j = (int)floor(atomCoords[iter][2] / hashSpacing) % L;
        k = (int)floor(atomCoords[iter][3] / hashSpacing) % L;
        //next, place it into the grid if that locations value is 0, else send it to grid_chain
        if(hashGrid[i][j][k] == 0) {
            hashGrid[i][j][k] = iter;
        } else {
            chainLoc = hashGrid[i][j][k];
            while(gridChain[chainLoc] != 0) {
                chainLoc = gridChain[chainLoc];
            }
            gridChain[chainLoc] = iter;
        }
    }
}


//destroys hash grid as well as cluster labels
void destroyHash(){
    for(int i=0; i<L; i++){
        for(int j=0; j<L; j++){
            delete [] hashGrid[i][j];
        }
        delete [] hashGrid[i];
    }
    delete [] hashGrid;
    delete [] gridChain;
}

//destroys hash grid as well as cluster labels
void destroyFreeHash(){
    for(int i=0; i<free_x_dims; i++){
        for(int j=0; j<free_y_dims; j++){
            delete [] free_hashGrid[i][j];
        }
        delete [] free_hashGrid[i];
    }
    delete [] free_hashGrid;
    for(int i=0; i<max_free; i++){
        delete [] free_gridChain[i];
    }
    delete [] free_gridChain;
}


//if this returns 0, then the gridpoint is farther than input distance (sqDistCheck) from atom, returns 1 then within range, cell is occupied
//TODO: isOccupied calls this but it needs to actually pass in the correct rss squared from the vdw info
int checkDistance(float x1, float y1, float z1, float x2, float y2, float z2, float distToCheck) {
    float distSq = (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2);
    float sqDistCheck = distToCheck*distToCheck;
    if (distSq <= sqDistCheck){ //if the gridpoint is within VDW radii, we return positive
        return 1;
    } else { 
        return 0;
    }
}

//returns the magnitude of a vector
float mag(float x, float y, float z){
    return sqrt(x*x+y*y+z*z);
}

//if 0 is returned, the physical space is occupied by an atom, else 1 is returned and it is either solvent or void
//TODO: find some way to check that these arrays are correctly brought in, since they're on the heap
int isNotOccupied(float x, float y, float z){
    //location of physical gridpoint within our hash grid
    int i = ((int)floor(x / hashSpacing) % L)+L; //add over the mod value to get out of the range of negatives, which mess up mod function
    int j = ((int)floor(y / hashSpacing) % L)+L;
    int k = ((int)floor(z / hashSpacing) % L)+L;
    int atom;
    int bins  = ceil(maxR_ss / hashSpacing);
    //iterating over all 27 surrounding hash lattice cubes searching for atoms
    for(int i2 = i-bins; i2 <= i+bins; i2++) {
        for(int j2 = j-bins; j2 <= j+bins; j2++) {
            for(int k2 = k-bins; k2 <= k+bins; k2++) {
                atom = hashGrid[i2 % L][j2 % L][k2 % L]; 
                if(atom != 0) { //unless there's something in hashGrid here, we move on
                    while(atom != 0){ //start checking through gridChain
                        if(checkDistance(x,y,z,atomCoords[atom][1],atomCoords[atom][2],
                           atomCoords[atom][3],atomCoords[atom][0])==1){ //atomCoords[atom][0] was 0, but now we have it working!!
                            return 0; break;}
                        atom = gridChain[atom];
                    }  
                }
            }
        }
    }
    return 1; //if no confirmed checks, then the physical space is unoccupied
}

int isBoundaryWater(float x, float y, float z){
    //location of physical gridpoint within our hash grid
    int i = ((int)floor(x / hashSpacing) % L)+L; //add over the mod value to get out of the range of negatives, which mess up mod function
    int j = ((int)floor(y / hashSpacing) % L)+L;
    int k = ((int)floor(z / hashSpacing) % L)+L;
    int atom;
    int bins = ceil(R_bw / hashSpacing); //since our hash grid is spaced hashSpacing, we check out this many in each direction
    //iterating over all 27 surrounding hash lattice cubes searching for atoms
    for(int i2 = i-bins; i2 <= i+bins; i2++) {
        for(int j2 = j-bins; j2 <= j+bins; j2++) {
            for(int k2 = k-bins; k2 <= k+bins; k2++) {
                atom = hashGrid[i2 % L][j2 % L][k2 % L]; 
                if(atom != 0) { //unless there's something in hashGrid here, we move on
                    while(atom != 0){ //start checking through gridChain
                        if(checkDistance(x,y,z,atomCoords[atom][1],atomCoords[atom][2],
                           atomCoords[atom][3],R_bw)==1){ //atomCoords[atom][0] was 0, but now we have it working!!
                            return 1; break;}
                        atom = gridChain[atom];
                    }  
                }
            }
        }
    }
    return 0; //if no confirmed checks, then the physical space is irrelevant, counts as bulk solvent
}

int isFreeVolume(float x, float y, float z, float r){
    //location of physical gridpoint within our hash grid
    int i = ((int)floor(x / hashSpacing) % L)+L; //add over the mod value to get out of the range of negatives, which mess up mod function
    int j = ((int)floor(y / hashSpacing) % L)+L;
    int k = ((int)floor(z / hashSpacing) % L)+L;
    int atom;
    int bins = ceil((r + maxR_ss) / hashSpacing);
    //iterating over all 27 surrounding hash lattice cubes searching for atoms
    for(int i2 = i-bins; i2 <= i+bins; i2++) {
        for(int j2 = j-bins; j2 <= j+bins; j2++) {
            for(int k2 = k-bins; k2 <= k+bins; k2++) {
                atom = hashGrid[i2 % L][j2 % L][k2 % L];
                if(atom != 0) { //unless there's something in hashGrid here, we move on
                    while(atom != 0){ //start checking through gridChain
                        if(checkDistance(x,y,z,atomCoords[atom][1],atomCoords[atom][2],
                           atomCoords[atom][3],atomCoords[atom][0]+r)==1){ //if this gridpoint is within distance of an atom
                            return 0; //there is an atom within, so it is NOT free volume
                        }
                        atom = gridChain[atom];
                    }  
                }
            }
        }
    }
    
    return 1; //if no atoms within R_check, this is indeed free volume!
}

//returns -1 if microvoid and 1 if solvent accessible void
int voidType(float x, float y, float z){
    float test[4]; //location of test particle
    test[1] = x; test[2] = y; test[3] = z;//uses indices 1,2,3... index 0 left blank (for purpose of lining up with atomCoords)
    //location of physical gridpoint within our hash grid
    int i = ((int)floor(x / hashSpacing) % L)+L; //add over the mod value to get out of the range of negatives, which mess up mod function
    int j = ((int)floor(y / hashSpacing) % L)+L;
    int k = ((int)floor(z / hashSpacing) % L)+L;
    int atom;
    int bins = ceil((2*R_probe+maxR_ss) / hashSpacing);
    int atomsWithin[1000] = {0}; //contains list of atoms within distance R_ss+2*R_probe (1000 should be large enough)
    int maxAtomsWithin = 0;//keeps track of how many indices are used in this array
    //iterating over all 27 surrounding hash lattice cubes searching for atoms
    for(int i2 = i-bins; i2 <= i+bins; i2++) {
        for(int j2 = j-bins; j2 <= j+bins; j2++) {
            for(int k2 = k-bins; k2 <= k+bins; k2++) {
                atom = hashGrid[i2 % L][j2 % L][k2 % L];
                if(atom != 0) { //unless there's something in hashGrid here, we move on
                    while(atom != 0){ //start checking through gridChain
                        if(checkDistance(x,y,z,atomCoords[atom][1],atomCoords[atom][2],
                           atomCoords[atom][3],atomCoords[atom][0]+(2*R_probe))==1){ //if this gridpoint is within distance of an atom
                            atomsWithin[maxAtomsWithin] = atom; maxAtomsWithin++; //store that atoms index
                        }
                        atom = gridChain[atom];
                    }  
                }
            }
        }
    }
    if(maxAtomsWithin == 0){//if there were no atoms within distance R_ss+(2*R_probe)
        return 1; //if no confirmed checks, then this is a solvent accessible void
    } else{//if there were atoms within this distance, we set up an elastic spring model in case the gridpoint is off center
            
        
        //need to define variables and then write method
        //define variables: E_spring, F_net, F_net_vec, dir, r, unit_vec, dr_vec, viscosity
        float F_net = 0; float F_netVec[4] = {0};//float F_netVec[4]; //float dir[4]; //F_net = F_netVec/dir 
        float prevF_net = 0; float step = 0; //check each pair of F_net values to prevent exact reflection loops
        float r_mag = 0; float r[4] = {0}; float unit_vec[4] = {0}; //unit vec = r_vec/||r|| and r is distance from test point to atom
        float dr_mag = 0; float dr[4] = {0}; float dr_unit[4] = {0}; float prevdr_mag = max_dr_mag; //distance that test will move at each step
        float E_spring = 0; int distCheck;
        //float decay = 1.0;
        int iter = 0;

        do {
            //cout << iter << " ";
            F_netVec[1] = 0; F_netVec[2] = 0; F_netVec[3] = 0;
            for(int count = 0; count < maxAtomsWithin; count++){
                //find distance from test point to atom, magnitude, vector, and unit vector
                r_mag = mag(atomCoords[atomsWithin[count]][1]-test[1],atomCoords[atomsWithin[count]][2]-test[2],atomCoords[atomsWithin[count]][3]-test[3]);
                if(r_mag < R_probe + atomCoords[atomsWithin[count]][0]){ //half spring, so only exerts force if spring is compressed
                    for(int indx = 1; indx <= 3; indx++){
                        r[indx] = atomCoords[atomsWithin[count]][indx]-test[indx];
                        unit_vec[indx] = r[indx] / (r_mag+FLT_MIN);
                        F_netVec[indx] = F_netVec[indx] + (r_mag -(R_probe+atomCoords[atomsWithin[count]][0]))*unit_vec[indx];
                        //cout << F_netVec[indx] << " ";
                    }
                }
            }
            F_net = mag(F_netVec[1],F_netVec[2],F_netVec[3]); //get the final magnitude of force for loop check
            //get F_net and then use it to move our test particle: test = test+dr

            if(abs(F_net - prevF_net) < F_zero){ //if the two forces are nearly identical in magnitude
                step = 0.5;
            }else{
                step = 1.0;
            }
            dr_mag = viscosity*step*F_net;
            
            if(dr_mag > max_dr_mag){
                dr_mag = max_dr_mag; //if it wants to move extremely far away, there's a problem, so we don't let it do that.
            }
            //if it starts to move further again, then it may be stuck in an oscillatory loop
            if(dr_mag >= prevdr_mag){
                dr_mag = dr_mag *0.1*(1.0+((float)rand()/RAND_MAX));////so we break it out of that loop randomly (magnitude from 0.1 to 0.2 * dr_mag)
            }

            for(int indx = 1; indx <= 3; indx++){
                dr_unit[indx] = F_netVec[indx] / (F_net+FLT_MIN); //dr and F_net point along same line
                //cout << dr_unit[indx] << " ";
                dr[indx] = dr_mag * dr_unit[indx];
                test[indx] = test[indx] + dr[indx];
            }
            if(checkDistance(x,y,z,test[1],test[2],
            test[3],targetR_probe)==0){
                return -1; //must be a microvoid
            } //if 0, then it has moved too far away
            
            prevF_net = F_net;
            prevdr_mag = dr_mag;
            iter++;
            //cout << iter << " " << F_net << " " << dr_mag <<" " << F_netVec[1] << " " << F_netVec[2] << " " <<F_netVec[3] << "\n";
            
            //cout << iter << " " << dr_mag/a_grid << "\n";
            //assert(iter < maxIter);
            //cout << "\n";
        } while(dr_mag > min_dr_mag && iter < maxIter);
        
        if(checkDistance(x,y,z,test[1],test[2],
            test[3],R_probe)==0){ //if the test particle is farther than Rprobe from gridpoint)

            return -1; //this grid point is not part of void space (so microvoid)
        }
        
        //cout << maxAtomsWithin << " " << iter << "\n";
        //find E_spring after we do these iterations
        E_spring = 0;
        for(int count = 0; count < maxAtomsWithin; count++){
            r_mag = mag(atomCoords[atomsWithin[count]][1]-test[1],atomCoords[atomsWithin[count]][2]-test[2],atomCoords[atomsWithin[count]][3]-test[3]);
            if(r_mag < R_probe+atomCoords[atomsWithin[count]][0]){
                E_spring = E_spring + (r_mag - (R_probe+atomCoords[atomsWithin[count]][0]))*(r_mag - (R_probe+atomCoords[atomsWithin[count]][0]));
                //cout<<E_spring << "\n";
            }
        }
        
        if(E_spring >= E_zero){ //stressed equilibrium
            return -1; //microvoid space
        }
        //cout << x << " " << y << " " << z << "\n";
        return 1; //grid point volume is part of void space 
    }
}

//returns -1 if microvoid and 1 if solvent accessible void
int withinFreeVolume(float x, float y, float z){
    int free_vol;
    int i = ((int)floor(x / R_probe) % free_x_dims)+free_x_dims; //add over the mod value to get out of the range of negatives, which mess up mod function
    int j = ((int)floor(y / R_probe) % free_y_dims)+free_y_dims;
    int k = ((int)floor(z / R_probe) % free_z_dims)+free_z_dims;
    //iterating over all 27 surrounding hash lattice cubes searching for atoms
    for(int i2 = i-1; i2 <= i+1; i2++) {
        for(int j2 = j-1; j2 <= j+1; j2++) {
            for(int k2 = k-1; k2 <= k+1; k2++) {
                free_vol = free_hashGrid[i2 % free_x_dims][j2 % free_y_dims][k2 % free_z_dims];
                if(free_vol != 0) { //unless there's something in hashGrid here, we move on
                    while(free_vol != 0){ //start checking through gridChain
                        if(checkDistance(x,y,z,free_gridChain[free_vol][1],free_gridChain[free_vol][2],
                           free_gridChain[free_vol][3],R_probe)==1){ //if this gridpoint is within distance of another free volume
                            return 1; //this is within R_probe distance of free volume, so its void
                        }
                        free_vol = free_gridChain[free_vol][0];
                    }  
                }
            }
        }
    }
    return 0; //not within R_probe distance of any free volume, so microvoid
}

//returns -1 if microvoid and 1 if solvent accessible void
int notWithinFreeVolume(float x, float y, float z){
    int free_vol;
    int i = ((int)floor(x / R_probe) % free_x_dims)+free_x_dims; //add over the mod value to get out of the range of negatives, which mess up mod function
    int j = ((int)floor(y / R_probe) % free_y_dims)+free_y_dims;
    int k = ((int)floor(z / R_probe) % free_z_dims)+free_z_dims;
    //iterating over all 27 surrounding hash lattice cubes searching for atoms
    for(int i2 = i-1; i2 <= i+1; i2++) {
        for(int j2 = j-1; j2 <= j+1; j2++) {
            for(int k2 = k-1; k2 <= k+1; k2++) {
                free_vol = free_hashGrid[i2 % free_x_dims][j2 % free_y_dims][k2 % free_z_dims];
                if(free_vol != 0) { //unless there's something in hashGrid here, we move on
                    while(free_vol != 0){ //start checking through gridChain
                        if(checkDistance(x,y,z,free_gridChain[free_vol][1],free_gridChain[free_vol][2],
                           free_gridChain[free_vol][3],R_cap)==1){ //if this gridpoint is within distance of another free volume
                            return 0; //cannot be sure, requires spring model
                        }
                        free_vol = free_gridChain[free_vol][0];
                    }  
                }
            }
        }
    }
    return 1; //microvoid for certain
}

//returns 1 for void and -1 for microvoid
int hybridVoidType(float x, float y, float z){
    if(isFreeVolume(x,y,z,R_probe)){
        count1++;
        return 1; //definitely void
    }
    else if(withinFreeVolume(x,y,z)){
        count2++;
        return 1; //void
    }
    else if(notWithinFreeVolume(x,y,z)){
        count3++;
        return -1; //microvoid
    }
    else{ //nothing worked, must revert to spring model
        count4++;
        return voidType(x,y,z);
    }
}

// uf_find returns the canonical label for the equivalence class containing x
int uf_find(int x) {
  int y = x;
  while (clusterLabel[y] != y){
    y = clusterLabel[y];
    //cout << y << "\n";
  }
  while (clusterLabel[x] != x) {
    int z = clusterLabel[x];
    clusterLabel[x] = y;
    x = z;
  }
  return y;
}

//uf_union joins two equivalence classes and returns the canonical label of the resulting class.
int uf_union(int x, int y) {
  return clusterLabel[uf_find(x)] = uf_find(y);
}

//uf_make_set creates a new equivalence class and returns its label
int uf_make_set(int type) {
  maxCluster ++;
  assert(maxCluster < n_labels);
  clusterLabel[maxCluster] = maxCluster; //may be unnecessary since initialized to clusterLabel[n]=n]
  clusterSize[maxCluster] = type; //add the first point into the new cluster, depending on void type
  //cout << maxCluster<<"\n";
  return maxCluster;
}

int maximum(int x, int y, int z) {
	int max = x; /* assume x is the largest */

	if (y > max) { /* if y is larger than max, assign y to max */
		max = y;
	} /* end if */

	if (z > max) { /* if z is larger than max, assign z to max */
		max = z;
	} /* end if */

	return max; /* max is the largest value */
} /* end function maximum */

//returns 1 if positive, -1 if negative, 0 if 0
int sign(int x){
    if (x < 0){
        return -1;
    }else{
        return 1;
    }
}


void addAtomToPDB(float x, float y, float z, int type){
    //need two MyAtom objects, one for void and one for microvoid?
    if(type == 1){ //void
        voidAtomNum++;
        if(voidAtomNum % 10000 == 0){
            voidAtom->chain_id = (voidAtom->chain_id)+1;
        }
        voidAtom->number = voidAtomNum;
        voidAtom->res_num = voidAtomNum % 10000;
        voidAtom->x = x;
        voidAtom->y = y;
        voidAtom->z = z;
        voidAtom->write(voidVis, false);
        voidVis << "\n";
    }else if(type == -1){ //microvoid
        microvoidAtomNum++;
        if(microvoidAtomNum % 10000 == 0){
            microvoidAtom->chain_id = (microvoidAtom->chain_id)+1;
        }
        microvoidAtom->number = microvoidAtomNum;
        microvoidAtom->res_num = microvoidAtomNum % 10000;
        microvoidAtom->x = x;
        microvoidAtom->y = y;
        microvoidAtom->z = z;
        microvoidAtom->write(microvoidVis, false);
        microvoidVis << "\n";
    }
}

//runs HK algorithm for a new slice, percolating based off of previous slice in z direction
//fills slice with values 0,1, and various other ints representing protein, solvent, and void respectively
int** propagate(int** prevSlice, int** slice, float z, int xIntervals, int yIntervals){
    int neighbors[connect_num]; //these will hold occupation values of 3 adjacent cells
    int nType[connect_num]; //stores void or microvoid type of neighbor
    int hasJoined; //this will be used to determine what to do for union-find each time through
    int type; //determines solvent accessible or microvoid
    for(int i=1; i<=xIntervals; i++){ //starting at gridpoint one over from edge in x and y, because the outer edges are just solvent
        for(int j=1; j<=yIntervals; j++){
             slice[i][j] = isNotOccupied(i*a_grid, j*a_grid, z); //value 0 if atom, 1 if void/solvent
             //begin HK, union find implementation here
             if (!!slice[i][j]) { //if there is NOT an atom here
                 //check whether it is void or microvoid
                 type = voidType(i*a_grid, j*a_grid, z);
                
                 
                //if sign positive, solvent accessible void, negative if microvoid
                neighbors[0] = prevSlice[i][j]; //neighbor below's cluster label, with sign denoting void or microvoid
                neighbors[1] = slice[i][j-1];//neighbor behind's cluster label; thinking of i as the left to right spacial dimension, (i.e. x) 
                neighbors[2] = slice[i-1][j];//and j as the up down, with the slices containing different depths
                nType[0] = sign(clusterSize[prevSlice[i][j]]); 
                nType[1] = sign(clusterSize[slice[i][j-1]]);
                nType[2] = sign(clusterSize[slice[i-1][j]]);
                if(connect_num == 7){
                    neighbors[3] = slice[i-1][j-1];
                    neighbors[4] = prevSlice[i-1][j];
                    neighbors[5] = prevSlice[i][j-1];
                    neighbors[6] = prevSlice[i-1][j-1];
                    nType[3] = sign(clusterSize[slice[i-1][j-1]]);
                    nType[4] = sign(clusterSize[prevSlice[i-1][j]]);
                    nType[5] = sign(clusterSize[prevSlice[i][j-1]]);
                    nType[6] = sign(clusterSize[prevSlice[i-1][j-1]]);
                 }
                 
                 //next, depending on the case, we add this new cell to its appropriate cluster
                 hasJoined = -1; //this flag becomes 1 once this cell has been added into a cluster
                 for(int k=0; k<connect_num; k++){
                     if(neighbors[k] != 0){ //if there isn't an atom in the neighbor cell, its a void
                         if(hasJoined == -1){ //if we haven't yet put this gridpoint into another neighbor's cluster
                             if(type == nType[k]){ //if our gridpoint is the same voidtype as the neighbor, they join
                                 //cout << neighbors[k] << " ";
                                 
                                 
                                slice[i][j] = neighbors[k]; //this gridpoint joins that cluster
                                if(nType[k]>0){
                                    if(uf_find(neighbors[k]) == uf_find(1)){//if cluster 1, check to see if its boundary water
                                        if(!!isBoundaryWater(i*a_grid, j*a_grid, z)){ //only needs check here because other voids MUST be within Rbw
                                            clusterSize[1]++; //its boundary water
                                        } 
                                    }else{
                                        clusterSize[neighbors[k]]++; //add the solvent accessible void into our counter
                                        addAtomToPDB(i*a_grid+(xshift - solvGap),j*a_grid+(yshift - solvGap),z+(zshift - solvGap),type);
                                    }
                                } else{
                                   clusterSize[neighbors[k]]--; //add the microvoid point into our counter
                                   addAtomToPDB(i*a_grid+(xshift - solvGap),j*a_grid+(yshift - solvGap),z+(zshift - solvGap),type);
                                }
                                //cout << neighbors[k] << " " << clusterSize[neighbors[k]] << "\n";
                                hasJoined = neighbors[k]; //the cell has now joined a cluster of neighbor index i
                             }
                         } else{//if the gridpoint has already joined a cluster
                             if(type == nType[k]){ //and the cluster its in is the same type as this neighbor's cluster type
                                 if(hasJoined != neighbors[k]){ //if they're not the same cluster already
                                   uf_union(hasJoined, neighbors[k]);//join the two clusters
                                 }
                             }
                             
                         }
                     }
                         
                 }
                 if(hasJoined == -1){ //if it still hasn't been added, make a new cluster
                     slice[i][j] = uf_make_set(type);
                     addAtomToPDB(i*a_grid+(xshift - solvGap),j*a_grid+(yshift - solvGap),z+(zshift - solvGap),type);
                     //cout << slice[i][j] << " " << clusterSize[slice[i][j]] << "\n";
                 }
                         
                         

             } else{ //if there is a protein atom at this gridpoint
                 clusterSize[0]++; //clusterSize[0] holds total number of protein atoms
             }
             //cout<< slice[i][j] << " ";
        } //cout << "\n";
    } //cout << "\n\n\n";
    return slice;
}

//this subroutine sweeps through our lattice, searching for points defined to be free volume, stores them in hash
void freeVolSweep(int * n){
    clock_t freeVolSweepTime;
    freeVolSweepTime = clock();
    //we will want to initialize our hash grid
    //spacing should be R_check (could it be R_probe?)
    generateFreeHash(n); //we are assuming that this has made our hash grid of proper length
    //now we need to sweep through each gridpoint and see if its free volume or not
    int chainLoc;
    int free_vol_count = 0; //this will be incremented for each free volume we find
    int free_i; int free_j; int free_k;
    
    for(int k=0; k<=n[3]; k++){ //sweep through z direction (IS THIS THE PROPER LENGTH?)
        for(int j=0; j<=n[2]; j++){ //sweep through z direction
            for(int i=0; i<=n[1]; i++){ //sweep through z direction
                if(!!isFreeVolume(i*a_grid,j*a_grid,k*a_grid,R_probe) && !isFreeVolume(i*a_grid,j*a_grid,k*a_grid,R_cap)){ //if free volume but not bulk solvent
                    free_vol_count++;
                    free_i = floor(i*a_grid / R_probe);
                    free_j = floor(j*a_grid / R_probe);
                    free_k = floor(k*a_grid / R_probe);
                    if(free_hashGrid[free_i][free_j][free_k] == 0) {
                        free_hashGrid[free_i][free_j][free_k] = free_vol_count;
                        //cout << free_hashGrid[free_i][free_j][free_k] << "\n";
                    } else {
                        chainLoc = free_hashGrid[free_i][free_j][free_k];
                        while(free_gridChain[chainLoc][0] != 0) {
                            chainLoc = free_gridChain[chainLoc][0];
                        }
                        free_gridChain[chainLoc][0] = free_vol_count;
                    }
                    free_gridChain[free_vol_count][1] = i*a_grid;
                    free_gridChain[free_vol_count][2] = j*a_grid;
                    free_gridChain[free_vol_count][3] = k*a_grid;
                }
            }
        }
    }
    freeVolSweepTime = clock() - freeVolSweepTime;
    cout << "Free Volume Points: " << free_vol_count << "\n";
    printf ("Free Volume Sweep time: %d clicks (%f seconds).\n",(int)freeVolSweepTime,((float)freeVolSweepTime)/CLOCKS_PER_SEC);
}

//runs HK algorithm
//can choose to return value of volume when running time statistics
void hoshenKopelman(){
    //intervals in each coordinate (for grid construction)
    int n[4]; //won't be using n[0]
    for(int i=1; i<=3; i++){
        n[i] = ceil(abs(atomCoords[0][i]/a_grid)) + 2*(solvGap/a_grid); //4 intervals on each side of min and max for solvent (2 Angstroms on each side)
    }
    cout << "Dimensionality: " << n[1] << " x " << n[2] << " x " << n[3] << "\n";
   
    //freeVolSweep(n);
    
    //begin HK
    //write function to check if a grid point is occupied
    //this will be called in the loop for every point we consider
    //we will run along the z axis, scanning in pairs using a propagate function
    int** slice1 = (int **) new int* [n[1]+1]; //this creates [n[1]+1] * [n[2]+1] dimension slices
    int** slice2 = (int **) new int* [n[1]+1];
    for(int i=0; i<=n[1]; i++){
        slice1[i] = (int *) new int [n[2]+1];
        slice2[i] = (int *) new int [n[2]+1];
        for(int j=0; j<=n[2]; j++){
            slice1[i][j] = 1;
            slice2[i][j] = 1;
        }
    }
    
    //bad_alloc throw is here, if a_grid goes too small, this can become a problem for large proteins
    n_labels = (n[1] * n[2] * n[3]) / 2; //this is if cells were occupied in diag. only. max # of clusters possible
    //this is what we will iterate over (starts at 1 since we already know that it exists)
    clusterLabel = (int *) new int[n_labels]; //link list that holds which clusters are in the same equivalence class
    clusterSize = (int *) new int[n_labels]; //collects amount of points in each cluster
    for(int i=0; i<n_labels; i++){
        clusterSize[i] = 0;
        clusterLabel[i] = i;
    }
    clusterSize[1]=0;//n[1]*n[2] + (n[1]+n[2])*2*(2*ceil(n[3]/2.0)-1);
    //cout << (n[1]*n[2] + (n[1]+n[2])*2*(2*ceil(n[3]/2.0)-1))*pow(a_grid,3.0) << "\n"; //do we need to initialize it as 1?
    
    
    //slices have (# of intervals +1) gridpoints to account for endpoints
    //both slices are initialized to ones, so we will keep the outer rim of each at 1s forever, representing solvent along edges
    
    //we can run along z axis at least as long as n_z and when the point is reached where the slice is all 1, we break;
    
    //we will store protein volume in clusterSize[0]
    //likely will need to do something about solvent volume being
    float z1; float z2;
    for(int k=1; k<=ceil(n[3]/2.0)+1; k++){
        //the slices have constant z, but depend differently on the value of k
        z2 = ((2*k)-1)*a_grid; z1 = 2*k*a_grid;
        //alternate slices until we're done.
        slice2 = propagate(slice1, slice2, z2, n[1], n[2]);
        slice1 = propagate(slice2, slice1, z1, n[1], n[2]);
        
        //could throw in something that checks if slice1 equals all 1, to break out of the loop without wasting extra time!
        //need to count up after each run through, and then also run the final canonicalization at the end
    }
    for(int i=0; i<=n[1]; i++){
        delete [] slice1[i];
        delete [] slice2[i];
    }
    delete [] slice1; delete [] slice2;
    
}


//should eventually return something that can be piped out of the program for the final result per runthrough
float clusterSizeCounter(){
    int minVoidSize  = ceil((4.0/3.0)*M_PI*(R_probe*R_probe*R_probe)/(a_grid*a_grid*a_grid)); //smallest size a void can be, removes problem voids
    for(int i=1; i <= maxCluster; i++){
        if(clusterLabel[i] != i){ //if this is not the canonical cluster label for this equivalence class
            clusterSize[clusterLabel[i]] = clusterSize[clusterLabel[i]] + clusterSize[i]; //pour the storage from this cluster label into the next one up the chain towards canonical
            clusterSize[i] = 0;
        }
        //cout << i << " " << clusterLabel[i] << " " << clusterSize[i] << "\n";
    }
    clusterSize[1] = clusterSize[1] + clusterSize[uf_find(1)];//store all boundary water in cluster 1
    clusterSize[uf_find(1)] = 0;
    
    //start putting stuff in here involving finding largest cluster size, collecting distribution statistics, etc.
    int i = 2; int voidvol = 0; int microvoid = 0; 
    int largestVoidInd = 0;
    int largestMicrovoidInd = 0;
    int numberOfVoids = 0; int numberOfMicrovoids = 0;
    int numOfErrorGdpts = 0;
    //pull boundary water volume from its high valued canonical label back into label 1
    
    
    //sum up our cluster sizes
    while(i <= maxCluster){
        
        if(clusterSize[i] > 0){
            if(clusterSize[i] >= minVoidSize){ //if the void is big enough (no false positive voids that are actually microvoids)
                voidvol = voidvol + clusterSize[i];
                numberOfVoids++;
                if(clusterSize[i] > largestVoid){
                    largestVoidInd = i;
                    largestVoid = clusterSize[i];
                }
            }else{ //its a false positive
                numOfErrorGdpts = numOfErrorGdpts + clusterSize[i]; //add to number of error points
                clusterSize[i] = -1 * clusterSize[i]; //put it into microvoids by making it negative
                microvoid = microvoid - clusterSize[i]; //we add it to microvoid instead
                numberOfMicrovoids++;
                //no need to worry about it being the largest either
            }
        }
        else if(clusterSize[i] < 0){
            microvoid = microvoid - clusterSize[i]; //because negative
            numberOfMicrovoids++;
            if(clusterSize[i] < largestMicrovoid){
                largestMicrovoidInd = i;
                largestMicrovoid = clusterSize[i];
            }
        }
        
        i++;

    }
    largestMicrovoid = -1*largestMicrovoid;
   //final void counter:
    int totalVol = clusterSize[0] + microvoid + voidvol; //used to determine packing density.
    int boundaryWater = clusterSize[1];
    int proteinVol = clusterSize[0];
    cout << "R_probe: " << R_probe << "\n";
    cout << "Number of Voids: " << numberOfVoids << "\n";
    cout << "Number of Microvoids: " << numberOfMicrovoids << "\n";
    cout << "Protein Volume = " << clusterSize[0] << " gridpoints,   " << setprecision(3) << clusterSize[0]*pow(a_grid,3.0) << " angstroms cubed.\n";
    cout << "Boundary Water Volume = " << boundaryWater << " gridpoints,   " << setprecision(3) << clusterSize[1]*pow(a_grid,3.0) << " angstroms cubed.\n";
    cout << "Largest Void Volume = " << largestVoidInd << " " << setprecision(3) << largestVoid*pow(a_grid,3.0) << " angstroms cubed.\n";
    cout << "Largest Microvoid Volume = " << largestMicrovoidInd << " " << setprecision(3) << largestMicrovoid*pow(a_grid,3.0) << " angstroms cubed.\n";
    cout << "Void volume: " << voidvol << " gridpoints,   " << setprecision(3) << voidvol*pow(a_grid,3.0) << " angstroms cubed.\n";
    cout << "Microvoid volume: " << microvoid << " gridpoints,   " << setprecision(3) << microvoid*pow(a_grid,3.0) << " angstroms cubed.\n";
    cout << "Void + Microvoid: " << setprecision(3) << (microvoid+voidvol)*pow(a_grid,3.0) << " angstroms cubed.\n";
    cout << "Total Gridpoints: " << clusterSize[0]+clusterSize[1]+voidvol+microvoid << "\n";
    cout << "Packing Density: " << setprecision(3) << (float)clusterSize[0] / totalVol << "\n";
    cout << count1 << " " << count2 << " " << count3 << " " << count4 << "\n";
    

    //for final data run
    glob_prot_vol_gdpts[global_rot_iter] = clusterSize[0]; 
    glob_mv_gdpts[global_rot_iter] = microvoid; 
    glob_cav_gdpts[global_rot_iter] = voidvol; 
    glob_solv_gdpts[global_rot_iter] = boundaryWater;
    glob_mv_clusters[global_rot_iter] = numberOfMicrovoids; 
    glob_cav_clusters[global_rot_iter] = numberOfVoids; 
    glob_cav_cluster_errors[global_rot_iter] = numOfErrorGdpts;
    glob_largest_mv_gdpts[global_rot_iter] = largestMicrovoid; 
    glob_largest_cav_gdpts[global_rot_iter] = largestVoid; 
    glob_gdpts_total[global_rot_iter] = clusterSize[0]+clusterSize[1]+voidvol+microvoid; //void+microvoid+protein+solvent
    glob_prot_vol[global_rot_iter] = clusterSize[0]*pow(a_grid,3.0); 
    glob_mv_vol[global_rot_iter] = microvoid*pow(a_grid,3.0); 
    glob_cav_vol[global_rot_iter] = voidvol*pow(a_grid,3.0); 
    glob_solv_vol[global_rot_iter] = clusterSize[1]*pow(a_grid,3.0);
    glob_total_vol[global_rot_iter] = (clusterSize[0]+clusterSize[1]+voidvol+microvoid)*pow(a_grid,3.0); 
    glob_packing_density[global_rot_iter] = (float)clusterSize[0] / totalVol; //in this case, total volume is the void+microvoid+protein, EXCLUDING solvent accessibility
    
    //histogram generator
    //make this larger, then start storing in this continually, without making it over and over again
    for(int i = 2; i<=maxCluster; i++){
        if(clusterSize[i] > 0){
            voidHist[clusterSize[i]] ++;
        }
        if(clusterSize[i] < 0){
            microvoidHist[-1*clusterSize[i]] ++;
        }
    }
    
    //for largest cluster file
    largestMVFile << largestMicrovoid << "\n"; //ultimately a list of 100 largest microvoids
    
    maxCluster = 1; //need to restart this cycle

    delete [] clusterLabel; delete [] clusterSize;
    //return (float)clusterSize[0] / totalVol; 
    //float * totals = new float[3]; 
    //totals[0] = microvoid*pow(a_grid,3.0); 
    //totals[1] = voidvol*pow(a_grid,3.0); 
    //totals[2] = proteinVol*pow(a_grid,3.0);
    //return totals;
    //return value changes depending on what data we are collecting
    //return (microvoid+voidvol)*pow(a_grid,3.0);
    return (microvoid+voidvol+proteinVol)*pow(a_grid,3.0);
    
    //return (microvoid)*pow(a_grid,3.0);
    //return boundaryWater*pow(a_grid,3.0);
    //return (voidvol)*pow(a_grid,3.0);
    //return (float)clusterSize[0] / totalVol;
}

void histogramGenerator(){
    //ostringstream gridString; gridString << setprecision(3) << a_grid; 
    //ostringstream probeString; probeString << setprecision(3) << R_probe; 
    //voidHistFile.open((pdb_file_name+"-a_grid"+gridString.str()+"-Rprobe"+probeString.str()+"-voidHist.txt").c_str());
    //microvoidHistFile.open ((pdb_file_name+"-a_grid"+gridString.str()+"-Rprobe"+probeString.str()+"-microvoidHist.txt").c_str());
    for(int j = 1; j <=largestVoid; j++){
        if(voidHist[j] > 0){
            voidHistFile << j << " " << voidHist[j] << "\n";
        }
    }
    for(int j = 1; j <=largestMicrovoid; j++){
        if(microvoidHist[j] > 0){
            microvoidHistFile << j << " " << microvoidHist[j] << "\n";
        }
    }
    microvoidHistFile.close();
    voidHistFile.close();
    //print this stuff out to be used in matlab?
}

void histogramGenerator2(){
    for(int j = 1; j <=largestVoid; j++){
        if(voidHist[j] > 0){
            voidHistFile << j << " " << voidHist[j] << "\n";
        }
    }
    for(int j = 1; j <=largestMicrovoid; j++){
        if(microvoidHist[j] > 0){
            microvoidHistFile << j << " " << microvoidHist[j] << "\n";
        }
    }
    //print this stuff out to be used in matlab?
}

//returns mathematical and physical variables to a file
void reportVariables(){
    
    variableFile.open((pdb_file_name + "-variableFile.txt").c_str());
    
    variableFile << "a_grid: " << a_grid << "\n";
    variableFile << "R_probe: " << R_probe << "\n";
    variableFile << "R_bw: " << R_bw << "\n";
    variableFile << "Angles: " << setprecision(3) << angles[0] << ", " << angles[1] << ", " << angles[2] << "\n";
    variableFile << "max_dr_mag: " << max_dr_mag << "\n";
    variableFile << "min_dr_mag: " << min_dr_mag << "\n";
    variableFile << "viscosity: " << viscosity << "\n";
    variableFile << "maxIter: " << maxIter << "\n";
    if(maxR_ss == 2.105){
        variableFile << "OPLS data used.\n";
    } else{
        variableFile << "Bondi radii data used.\n";
    }
    
    variableFile.close();
}

//gives total volume and parameters of run in order, a_grid, R_probe, totalVolume
void totalVolumeRun(){
    clock_t algorithmTime;
    ofstream totalVolume;
    ostringstream gridString; gridString << setprecision(3) << a_grid; 
    ostringstream probeString; probeString << setprecision(3) << R_probe; 
    float vol;
    totalVolume.open((pdb_file_name+"-a_grid"+gridString.str()+"-Rprobe"+probeString.str()+"-totalVolume_new_con.txt").c_str());
    //basic algorithm
    algorithmTime = clock(); //begin timing our algorithm
    generateHash(); //make our hash coded map of atoms
    hoshenKopelman(); //propagate through our grid - percolation
    vol = clusterSizeCounter(); //process void and microvoid statistics
    algorithmTime = clock() - algorithmTime;
    printf ("Algorithm Processing time: %d clicks (%f seconds).\n",(int)algorithmTime,((float)algorithmTime)/CLOCKS_PER_SEC);
    destroyHash();
    totalVolume << a_grid << "\n";
    totalVolume << R_probe << "\n";
    totalVolume << vol;
    totalVolume.close();
}

//runs 100 random orientations of the protein with default Rprobe and a_grid
//generates cluster statistics based on the 100 runs
void clusterStatisticsRun(){   
    ostringstream gridString; gridString << setprecision(3) << a_grid; 
    ostringstream probeString; probeString << setprecision(3) << R_probe; 
    voidHistFile.open((pdb_file_name+"-a_grid"+gridString.str()+"-Rprobe"+probeString.str()+"-voidHist_new_con.txt").c_str());
    microvoidHistFile.open ((pdb_file_name+"-a_grid"+gridString.str()+"-Rprobe"+probeString.str()+"-microvoidHist_new_con.txt").c_str());
    clock_t algorithmTime;
    generateHash();
    int iters = 100;
    for(int i=0; i<iters; i++){ //we can run it 100 times
        algorithmTime = clock();
        generateHash();
        hoshenKopelman(); //propagate through our grid - percolation
        clusterSizeCounter(); //process void and microvoid statistics
        algorithmTime = clock() - algorithmTime;
        rotateCoords(); //rotate for the next time through to run again with new values
        destroyHash(); //since we will be making it again, we need to destroy the memory first
        histogramGenerator2();
    }   
    totalVolumeRun(); //get our volume file as well
    
    microvoidHistFile.close();
    voidHistFile.close();
}

void finalDataCollectionRun(){
    ostringstream gridString; gridString << setprecision(3) << a_grid; 
    ostringstream ratioString; ratioString << setprecision(3) << R_probe / a_grid; 
    voidHistFile.open((pdb_file_name+"-a"+gridString.str()+"-ratio"+ratioString.str()+"-"+connectivity_type+"-cvcstat.txt").c_str());
    microvoidHistFile.open((pdb_file_name+"-a"+gridString.str()+"-ratio"+ratioString.str()+"-"+connectivity_type+"-mvcstat.txt").c_str());
    summaryFile.open((pdb_file_name+"-a"+gridString.str()+"-ratio"+ratioString.str()+"-"+connectivity_type+"-summary.txt").c_str());
    largestMVFile.open((pdb_file_name+"-a"+gridString.str()+"-ratio"+ratioString.str()+"-"+connectivity_type+"-mvLcsze.txt").c_str());  
    
    clock_t algorithmTime;
    for(global_rot_iter=0; global_rot_iter<global_rotation_max; global_rot_iter++){ //we can run it 100 times
        algorithmTime = clock();
        generateHash();
        hoshenKopelman(); //propagate through our grid - percolation
        clusterSizeCounter(); //process void and microvoid statistics, prints out largest mv to largest mv cluster file
        algorithmTime = clock() - algorithmTime; //total alg. time includes hash generation, HK algorithm, and final cluster counting
        glob_cpu_time[global_rot_iter] = ((float)algorithmTime)/CLOCKS_PER_SEC;
        rotateCoords(); //rotate for the next time through to run again with new values
        destroyHash(); //since we will be making it again, we need to destroy the memory first
    }
    histogramGenerator(); // using aggregate cluster statistics instead of single cluster statistics, print to two hist files
    
    //avgs
    float prot_vol_gdpts_avg = 0; float mv_gdpts_avg = 0; float cav_gdpts_avg = 0;
    float solv_gdpts_avg = 0; float mv_clusters_avg = 0; float cav_clusters_avg = 0;
    float cav_cluster_error_avg = 0; float largest_mv_gdpts_avg = 0;
    float largest_cav_gdpts_avg = 0; float gdpts_total_avg = 0;
    float prot_vol_avg = 0; float mv_vol_avg = 0;
    float cav_vol_avg = 0; float solv_vol_avg = 0;
    float total_vol_avg = 0; float packing_density_avg = 0; float cpu_time_avg = 0;
    
    //stddv
    float prot_vol_gdpts_stderr = 0; float mv_gdpts_stderr = 0; float cav_gdpts_stderr = 0;
    float solv_gdpts_stderr = 0; float mv_clusters_stderr = 0; float cav_clusters_stderr = 0;
    float cav_cluster_error_stderr = 0; float largest_mv_gdpts_stderr = 0;
    float largest_cav_gdpts_stderr = 0; float gdpts_total_stderr = 0;
    float prot_vol_stderr = 0; float mv_vol_stderr = 0;
    float cav_vol_stderr = 0; float solv_vol_stderr = 0;
    float total_vol_stderr = 0; float packing_density_stderr = 0; float cpu_time_stderr = 0;
    
    //get avgs first by summing
    for(int i=0; i<global_rotation_max; i++){
        prot_vol_gdpts_avg += glob_prot_vol_gdpts[i];
        mv_gdpts_avg += glob_mv_gdpts[i];
        cav_gdpts_avg += glob_cav_gdpts[i];
        solv_gdpts_avg += glob_solv_gdpts[i];
        mv_clusters_avg += glob_mv_clusters[i];
        cav_clusters_avg += glob_cav_clusters[i];
        cav_cluster_error_avg += glob_cav_cluster_errors[i];
        largest_mv_gdpts_avg += glob_largest_mv_gdpts[i];
        largest_cav_gdpts_avg += glob_largest_cav_gdpts[i];
        gdpts_total_avg += glob_gdpts_total[i];
        prot_vol_avg += glob_prot_vol[i];
        mv_vol_avg += glob_mv_vol[i];
        cav_vol_avg += glob_cav_vol[i];
        solv_vol_avg += glob_solv_vol[i];
        total_vol_avg += glob_total_vol[i];
        packing_density_avg += glob_packing_density[i];
        cpu_time_avg += glob_cpu_time[i];
    }
    
    //then dividing them all by global_rotation_max to get average
    prot_vol_gdpts_avg = prot_vol_gdpts_avg/((float)global_rotation_max);
    mv_gdpts_avg = mv_gdpts_avg/((float)global_rotation_max);
    cav_gdpts_avg = cav_gdpts_avg/((float)global_rotation_max);
    solv_gdpts_avg = solv_gdpts_avg/((float)global_rotation_max);
    mv_clusters_avg = mv_clusters_avg/((float)global_rotation_max);
    cav_clusters_avg = cav_clusters_avg/((float)global_rotation_max);
    cav_cluster_error_avg = cav_cluster_error_avg/((float)global_rotation_max);
    largest_mv_gdpts_avg = largest_mv_gdpts_avg/((float)global_rotation_max);
    largest_cav_gdpts_avg = largest_cav_gdpts_avg/((float)global_rotation_max);
    gdpts_total_avg = gdpts_total_avg/((float)global_rotation_max);
    prot_vol_avg = prot_vol_avg/((float)global_rotation_max);
    mv_vol_avg = mv_vol_avg/((float)global_rotation_max);
    cav_vol_avg = cav_vol_avg/((float)global_rotation_max);
    solv_vol_avg = solv_vol_avg/((float)global_rotation_max);
    total_vol_avg = total_vol_avg/((float)global_rotation_max);
    packing_density_avg = packing_density_avg/((float)global_rotation_max);
    cpu_time_avg = cpu_time_avg/((float)global_rotation_max);
    
    //lastly, we calculate the standard error = corrected sample std dev / sqrt(n)
    for(int i=0; i<global_rotation_max; i++){
        prot_vol_gdpts_stderr += (glob_prot_vol_gdpts[i] - prot_vol_gdpts_avg)*(glob_prot_vol_gdpts[i] - prot_vol_gdpts_avg);
        mv_gdpts_stderr += (glob_mv_gdpts[i] - prot_vol_gdpts_avg)*(glob_mv_gdpts[i] - prot_vol_gdpts_avg);
        cav_gdpts_stderr += (glob_cav_gdpts[i] - cav_gdpts_avg)*(glob_cav_gdpts[i] - cav_gdpts_avg);
        solv_gdpts_stderr += (glob_solv_gdpts[i] - solv_gdpts_avg)*(glob_solv_gdpts[i] - solv_gdpts_avg);
        mv_clusters_stderr += (glob_mv_clusters[i] - mv_clusters_avg)*(glob_mv_clusters[i] - mv_clusters_avg);
        cav_clusters_stderr += (glob_cav_clusters[i] - cav_clusters_avg)*(glob_cav_clusters[i] - cav_clusters_avg);
        cav_cluster_error_stderr += (glob_cav_cluster_errors[i] - cav_cluster_error_avg)*(glob_cav_cluster_errors[i] - cav_cluster_error_avg);
        largest_mv_gdpts_stderr += (glob_largest_mv_gdpts[i] - largest_mv_gdpts_avg)*(glob_largest_mv_gdpts[i] - largest_mv_gdpts_avg);
        largest_cav_gdpts_stderr += (glob_largest_cav_gdpts[i] - largest_cav_gdpts_avg)*(glob_largest_cav_gdpts[i] - largest_cav_gdpts_avg);
        gdpts_total_stderr += (glob_gdpts_total[i] - gdpts_total_avg)*(glob_gdpts_total[i] - gdpts_total_avg);
        prot_vol_stderr += (glob_prot_vol[i] - prot_vol_avg)*(glob_prot_vol[i] - prot_vol_avg);
        mv_vol_stderr += (glob_mv_vol[i] - mv_vol_avg)*(glob_mv_vol[i] - mv_vol_avg);
        cav_vol_stderr += (glob_cav_vol[i] - cav_vol_avg)*(glob_cav_vol[i] - cav_vol_avg);
        solv_vol_stderr += (glob_solv_vol[i] - solv_vol_avg)*(glob_solv_vol[i] - solv_vol_avg);
        total_vol_stderr += (glob_total_vol[i] - total_vol_avg)*(glob_total_vol[i] - total_vol_avg);
        packing_density_stderr += (glob_packing_density[i] - packing_density_avg)*(glob_packing_density[i] - packing_density_avg);
        cpu_time_stderr += (glob_cpu_time[i] - cpu_time_avg)*(glob_cpu_time[i] - cpu_time_avg);
    }
    
    float correction_factor = (global_rotation_max - 1)*global_rotation_max; //to use corrected sample std dev
    
    //last step for standard error, divide by correction factor and take square root
    prot_vol_gdpts_stderr = pow(prot_vol_gdpts_stderr/correction_factor,0.5);
    mv_gdpts_stderr = pow(mv_gdpts_stderr/correction_factor,0.5);
    cav_gdpts_stderr = pow(cav_gdpts_stderr/correction_factor,0.5);
    solv_gdpts_stderr = pow(solv_gdpts_stderr/correction_factor,0.5);
    mv_clusters_stderr = pow(mv_clusters_stderr/correction_factor,0.5);
    cav_clusters_stderr = pow(cav_clusters_stderr/correction_factor,0.5);
    cav_cluster_error_stderr = pow(cav_cluster_error_stderr/correction_factor,0.5);
    largest_mv_gdpts_stderr = pow(largest_mv_gdpts_stderr/correction_factor,0.5);
    largest_cav_gdpts_stderr = pow(largest_cav_gdpts_stderr/correction_factor,0.5);
    gdpts_total_stderr = pow(gdpts_total_stderr/correction_factor,0.5);
    prot_vol_stderr = pow(prot_vol_stderr/correction_factor,0.5);
    mv_vol_stderr = pow(mv_vol_stderr/correction_factor,0.5);
    cav_vol_stderr = pow(cav_vol_stderr/correction_factor,0.5);
    solv_vol_stderr = pow(solv_vol_stderr/correction_factor,0.5);
    total_vol_stderr = pow(total_vol_stderr/correction_factor,0.5);
    packing_density_stderr = pow(packing_density_stderr/correction_factor,0.5);
    cpu_time_stderr = pow(cpu_time_stderr/correction_factor,0.5);
    
    //print off our data to the summary file
    summaryFile << a_grid << "\n";
    summaryFile << R_probe << "\n";
    summaryFile << prot_vol_gdpts_avg << " " << prot_vol_gdpts_stderr << "\n";
    summaryFile << mv_gdpts_avg << " " << mv_gdpts_stderr << "\n";
    summaryFile << cav_gdpts_avg << " " << cav_gdpts_stderr << "\n";
    summaryFile << solv_gdpts_avg << " " << solv_gdpts_stderr << "\n";
    summaryFile << mv_clusters_avg << " " << mv_clusters_stderr << "\n";
    summaryFile << cav_clusters_avg << " " << cav_clusters_stderr << "\n";
    summaryFile << cav_cluster_error_avg << " " << cav_cluster_error_stderr << "\n";
    summaryFile << largest_mv_gdpts_avg << " " << largest_mv_gdpts_stderr << "\n";
    summaryFile << largest_cav_gdpts_avg << " " << largest_cav_gdpts_stderr << "\n";
    summaryFile << gdpts_total_avg << " " << gdpts_total_stderr << "\n";
    summaryFile << prot_vol_avg << " " << prot_vol_stderr << "\n";
    summaryFile << mv_vol_avg << " " << mv_vol_stderr << "\n";
    summaryFile << cav_vol_avg << " " << cav_vol_stderr << "\n";
    summaryFile << solv_vol_avg << " " << solv_vol_stderr << "\n";
    summaryFile << total_vol_avg << " " << total_vol_stderr << "\n";
    summaryFile << packing_density_avg << " " << packing_density_stderr << "\n";
    summaryFile << cpu_time_avg << " " << cpu_time_stderr << "\n";

    //finally we close our files
    microvoidHistFile.close();
    voidHistFile.close();
    summaryFile.close();
    largestMVFile.close();
}

//runs 100 random orientations of the protein with default Rprobe and a_grid
//generates cluster statistics based on the 100 runs
void mv_clustersRun(){
    
    //want to generate results from several Rprobes and several a_grids
    //need to find out what subroutines we need to re run after each Rprobe change and after each a_grid change..
    //want to try a_grid = 0.2, 0.5, and 1.0
    //try Rprobe of 1.4, 0.5, 1.0, 2.0? 
    //exclude 1qfmH... (since it breaks at 2.0?)
    //this will finish off the story for microvoid and then we go on to processing void
    
    clock_t algorithmTime;
    generateHash();
    algorithmTime = clock();
    generateHash();
    hoshenKopelman(); //propagate through our grid - percolation
    clusterSizeCounter(); //process void and microvoid statistics
    algorithmTime = clock() - algorithmTime;
    rotateCoords(); //rotate for the next time through to run again with new values
    destroyHash(); //since we will be making it again, we need to destroy the memory first
    histogramGenerator();
    
    totalVolumeRun(); //get our volume file as well
}

void connectivityCompareRun(){
    ofstream cluster_compare;
    ostringstream gridString; gridString << setprecision(2) << a_grid; 
    ostringstream probeString; probeString << setprecision(2) << R_probe;
    cluster_compare.open((pdb_file_name+"-a_grid"+gridString.str()+"-Rprobe"+probeString.str()+"-new_connectivity.txt").c_str());
    float clusters = 0;
    float total = 0;
    
    for(int probeRep = 0; probeRep<=15;probeRep++){ //from R_probe=0.5 to 2.0
        for(int rep=0; rep<100; rep++){ //across 100 random rotations
            R_probe = (probeRep+5) * 0.1;
            generateHash();
            hoshenKopelman(); //propagate through our grid - percolation
            clusterSizeCounter(); //process void and microvoid statistics
            rotateCoords(); //rotate for the next time through to run again with new values
            destroyHash(); //since we will be making it again, we need to destroy the memory first
            clusters = clusters + global_mv_clusters;
            total = total + global_microvoid;
        }
        clusters = clusters / 100.0;
        total = total / 100.0;
        cluster_compare << R_probe << " " << total << " " << clusters << "\n"; //no need to put in volume as well since this is all for fixed a_grid
        clusters = 0;
        total = 0;
    }
}

/*
//generates void total over random orientations
void percThreshRun(){
    //REMEMBER TO CHANGE CLUSTER SIZE COUNTER BACK TO REGULAR FLOAT AFTER USING THIS
    //we will store in each file a line that has a_grid, R_probe, the void total, and the standard deviation
    ostringstream gridString; gridString << setprecision(2) << a_grid; 
    ostringstream iterString; iterString << setprecision(2) << statRunNum;
    ofstream voidTotalFile;
    voidTotalFile.open((pdb_file_name+"-a_grid"+gridString.str()+"-voidTotal."+iterString.str()+".txt").c_str());
    int rotations = 35;
    float * totals = new float[3];
    for(int rotCount = 0; rotCount <rotations; rotCount++){ //we will do our a_grid convergence test over 100 random orientations PER protein
        generateHash();
        R_probe = 0.01; targetR_probe = 0.0; int probeCount = 0;
        while(R_probe < 2.01){
            hoshenKopelman(); //propagate through our grid - percolation
            totals = clusterSizeCounter(); //outputs microvoid, void, and protein volume
            voidTotalFile << R_probe << " " << totals[0] << " " << totals[1] << " " << totals[2] << "\n";
            R_probe = R_probe + 0.01; targetR_probe = R_probe * 1.10;
        }
        rotateCoords();
        destroyHash();
    }
    voidTotalFile.close(); 
}*/

void basicRun(){
    clock_t algorithmTime;
    //basic algorithm
    algorithmTime = clock(); //begin timing our algorithm
    generateHash(); //make our hash coded map of atoms
    hoshenKopelman(); //propagate through our grid - percolation
    clusterSizeCounter(); //process void and microvoid statistics
    algorithmTime = clock() - algorithmTime;
    printf ("Algorithm Processing time: %d clicks (%f seconds).\n",(int)algorithmTime,((float)algorithmTime)/CLOCKS_PER_SEC);
    destroyHash();
}

void gridFactorRun(){
    clock_t algorithmTime;
    ofstream gridFactor_compare;
    float packing;
    gridFactor_compare.open((pdb_file_name+"-gridFactor_check.txt").c_str());
    for(int rep=0; rep<=100; rep++){
        //basic algorithm
        grid_factor = -1+((rep / 100.0)*2);
        algorithmTime = clock(); //begin timing our algorithm
        generateHash(); //make our hash coded map of atoms
        hoshenKopelman(); //propagate through our grid - percolation
        packing = clusterSizeCounter(); //process void and microvoid statistics
        algorithmTime = clock() - algorithmTime;
        printf ("Algorithm Processing time: %d clicks (%f seconds).\n",(int)algorithmTime,((float)algorithmTime)/CLOCKS_PER_SEC);
        destroyHash();
        destroyFreeHash();
        gridFactor_compare << grid_factor << " " << packing << "\n";
    }
    gridFactor_compare.close();
}

void visualizationRun(){
    //start following prescription based on hollow
    //load exactly like this but with x, y, and z shifted by xshift, yshift, zshift
    //ATOM      1  O   HOH A   1      16.447  60.738   8.160  0.00  0.00
    //prepare our Atom object      
    strncpy(voidAtom->atom_name, "O  ", 3);     
    strncpy(voidAtom->res_sname, "HOH", 3);
    char chain[2] = "A";
    voidAtom->chain_id = chain[0];
    strncpy(microvoidAtom->atom_name, "O  ", 3);     
    strncpy(microvoidAtom->res_sname, "HOH", 3);
    microvoidAtom->chain_id = chain[0];
    
  
    clock_t algorithmTime;

    voidVis.open((pdb_file_name+"-void.pdb").c_str());
    microvoidVis.open((pdb_file_name+"-microvoid.pdb").c_str());
    //basic algorithm
    algorithmTime = clock(); //begin timing our algorithm
    generateHash(); //make our hash coded map of atoms
    hoshenKopelman(); //propagate through our grid - percolation
    clusterSizeCounter(); //process void and microvoid statistics
    algorithmTime = clock() - algorithmTime;
    printf ("Algorithm Processing time: %d clicks (%f seconds).\n",(int)algorithmTime,((float)algorithmTime)/CLOCKS_PER_SEC);
    destroyHash();
    voidVis.close(); microvoidVis.close();
}

int main(int argc, char** argv) {

    
    clock_t importTime; 
    importTime = clock();
    
    //assign values from user input
    userInput(argc, argv);
    
   
    
    //produce atom coordinate array (global variable)
    getAtomCoords(pdb_file_name);
    firstQuadCoordShift();
    importTime = clock() - importTime;
    printf ("Protein import took %d clicks (%f seconds).\n",(int)importTime,((float)importTime)/CLOCKS_PER_SEC);
    
    //reportVariables();
    
    
    //which data analysis test we will be doing. regular runs use basicRun()
    //atomsVSTimeRun(((float)importTime)/CLOCKS_PER_SEC);
    //gridSizeRun();
    //orientationRun();
    //packDensVSRprobeRun();
    //inertiaRotateRun();
    //gridVolumeVStimeRun(); //return hoshenKopelman to a void function after running this
    //vdwRadiiRun();
    //clusterStatisticsRun();
    //voidVSa_gridRun();
    //probeSizeVSorientationRun();
    //totalVolumeRun();
    //a_grid_maxTestRun();
    //visualizationRun();
    //basicRun();
    //gridFactorRun();
    //voidTotalRun();
    //percThreshRun();
    //mv_clustersRun();
    //connectivityCompareRun();
    finalDataCollectionRun();
    
    return 0;
}
