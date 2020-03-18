#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <string.h>

#if !defined(OUTLIG_H)
    #define OUTLIG_H 1
void outlig(DOCKMol &);
void outpdb(DOCKMol &, ofstream &,const int );
bool Write_Mol2(DOCKMol &, std::ofstream &, string);
bool Write_Mol2_2(DOCKMol &mol, vector<float> &x,vector<float> &y, vector<float> &z,ofstream & ofs);
#endif 
