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
#if !defined(INITIALMOL_H)
    #define INITIALMOL_H 1
bool intialmol(vector <DOCKMol> &molvec, DOCKMol mol,float spherecenterx,float spherecentery,float spherecenterz,int initialnum,mt19937 &gen,AMBER_TYPER & amber,Energy_Score & c_nrg, string filevdw, string fileflex,string fileflex_drive_tbl,string grid_file,vector<TORSION> &torsions, bool flexible_flag);
bool intialmol2(vector<vector<float> > &ini_x, vector<vector<float> > &ini_y, vector<vector<float> > &ini_z, int num_molvec,DOCKMol mol,float spherecenterx,float spherecentery,float spherecenterz,int initialnum,mt19937 &gen);
using namespace std;
#endif
