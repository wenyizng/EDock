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

#if !defined(MC_SIMULATION_H)
    #define MC_SIMULATION_H 1

//bool MC(AMBER_TYPER & amber, DOCKMol & mol,Energy_Score & c_nrg, string filevdw, string fileflex,string fileflex_drive_tbl,string grid_file,int cutoff,RotatemolVec &temprotatemols,float temprature);
//bool MC(AMBER_TYPER & amber, DOCKMol &mol,Energy_Score & c_nrg, string filevdw, string fileflex,string fileflex_drive_tbl,string grid_file,int cutoff,vector<vector <float> > &,vector<vector <float> > &,vector<vector <float> > &,vector <float> &, float temprature,float angel,mt19937 &gen,vector<int> &MCsteps);
bool MC(AMBER_TYPER & amber, DOCKMol &mol,DOCKMol & receptor, Energy_Score & c_nrg, string filevdw, string fileflex,string fileflex_drive_tbl,string grid_file,float cutoff,vector<vector <float> > &x,vector<vector <float> > &y,vector<vector <float> > &z,vector <float> &tmpenergy,float temprature,mt19937 &gen,vector<float> &replica_MC_energy, bool flexible_flag,vector<vector<int>> biolip_matrix,int MC_steps,int & acceprate, vector<float> sphx,vector<float> sphy,vector<float> sphz,float bindsite_weight,vector<int> protein_atom_index,vector<vector <float> > ave, vector<vector <float> > std);
bool MC_flexible_ligand(AMBER_TYPER & amber, DOCKMol &mol, float temprature,mt19937 &gen,vector<vector<int>> biolip_matrix);
bool MC_rigid_ligand(AMBER_TYPER & amber, DOCKMol &mol, DOCKMol &receptor, Energy_Score & c_nrg, string filevdw, string fileflex,string fileflex_drive_tbl,string grid_file,float cutoff, float temprature,mt19937 &gen,vector<vector<int>> biolip_matrix,vector<float> sphx,vector<float> sphy,vector<float> sphz,float bindsite_weight,vector<int> protein_atom_index,vector<vector <float> > ave, vector<vector <float> > std);
vector<int> sort_index(RotatemolVec);
//bool Tmol(DOCKMol &, DOCKMol &);
//bool Rmol(DOCKMol &, DOCKMol &);
//bool generatationmol(DOCKMol &, DOCKMol & ,float,mt19937 &);
bool generatationmol(DOCKMol & original,DOCKMol & endmol,int method,mt19937 &gen);
//bool flexible_mol(DOCKMol & mol,mt19937 &gen,TORSION &torsion)
float flexible_mol(DOCKMol & mol, TORSION &torsion, float angle, AMBER_TYPER & amber,vector<vector<int>> biolip_matrix);
//bool generatationmol1(DOCKMol & original,DOCKMol & endmol,float angel,mt19937 &gen);
//bool generatationmol2(DOCKMol & original,DOCKMol & endmol,float angel,mt19937 &gen);
//void swap_conformations(DOCKMol &, DOCKMol &);
void swapmol(DOCKMol &, DOCKMol &);
using namespace std;
#endif
