#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <string.h>

#include "mol.h"
#if !defined(CLUSTER_H)
    #define CLUSTER_H 1
using namespace std;

bool cluster(vector<vector<float> > &x,vector<vector<float> > &y,vector<vector<float> > &z,vector<float> &energy,DOCKMol &mol,float cutoff, string output);
//float RMSD(DOCKMol mol1,DOCKMol mol2);
float RMSD(vector<float> &, vector<float> &,vector<float> &, vector<float> &,vector<float> &, vector<float> &);
int center(vector<int> &clusterflag,vector<vector<float> > &rmsdmatrix,int cluster);
bool combo(int cluster,vector<int> &clusterflag,vector<vector<float> > &x,vector<vector<float> > &y,vector<vector<float> > &z, vector<float> &combomolx,vector<float> &combomoly,vector<float> &combomolz);
int closec(vector<int> & clusterflag,int cluster,vector<float> &closcx,vector<float> &closcy,vector<float> &closcz,vector<vector<float> > &x,vector<vector<float> > &y,vector<vector<float> > &z);
float combo_rmsd(vector<float> &combomolx,vector<float> &combomoly,vector<float> &combomolz,vector<float> &x,vector<float> &y,vector<float> &z);
#endif
