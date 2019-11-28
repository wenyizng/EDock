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

#if !defined(PRECON_H)
    #define PRECON_H 1

bool readconformations(const char *filein,vector<vector<float> > &x,vector<vector<float> > &y,vector<vector<float> > &z, vector<float> &energy, DOCKMol &mol,vector<string> &namelog);
#endif
