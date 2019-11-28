#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <string.h>
#include <random>

#include "mol.h"
#include "flex_dock.h"
#if !defined(REMC_H)
    #define REMC_H 1

//bool REMC(AMBER_TYPER & amber, vector <DOCKMol> & mol,Energy_Score & c_nrg, string filevdw, string fileflex,string fileflex_drive_tbl,string grid_file,int cutoff,vector<vector <float> > &x,vector<vector <float> > &y, vector<vector <float> > &z,vector<float> &,float,float,string,mt19937 &);
bool sortmol(vector<DOCKMol> &, vector<int> &);

vector<int> randomgeneratenum(int num,int maxvalue);

bool REMC(AMBER_TYPER & amber, vector <DOCKMol> & REmolvec,Energy_Score & c_nrg, string filevdw, string fileflex,string fileflex_drive_tbl,string grid_file,int cutoff,vector<vector <float> > &x,vector<vector <float> > &y, vector<vector <float> > &z,vector<float> &accenergy,float Tmin, float Tmax,string outfile,mt19937 &gen, vector<TORSION> &torsions, bool flexible_flag,vector<float>sphx,vector<float>sphy,vector<float> sphz);


using namespace std;
#endif
