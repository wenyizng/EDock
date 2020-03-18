#ifndef TRANSLATE_H
#define TRANSLATE_H 

#include <vector>
#include <cmath>
#include <iostream>
#include <cstdlib> 
#include "mol.h"

const float PI=3.141592653589793;
const float Extra=1.0e-4;
const float UpMax=1.0e+10;

typedef vector<vector<float> > matrix;

using namespace std;

inline float deg2rad(float deg);

void translate_rotate(DOCKMol & mol);

bool GroupRotation1(const vector<float> &axisA, const vector<float> &axisB, float angle, vector<vector<float> > &pointB);

inline bool RotationMatrixB(const vector<float> &axis, float angle, matrix &romtx);

inline bool norm(const vector<float> &c, vector<float> &cc);

void ShowMyvector(const vector<float> &cc);

inline void SetMatrix(matrix &sm, int m, int n);

inline void subtract(const vector<float> &c1, const vector<float> &c2, vector<float> &cc);
inline void crossproduct(const vector<float> &c1, const vector<float> &c2, vector<float> &cc);
inline float innerproduct(const vector<float> &c1, const vector<float> &c2);

bool GroupRotation(const vector<float> &, const vector<float> &, float, DOCKMol &);

float Points2Dihedral(const vector<float> &c1, const vector<float> &c2, const vector<float> &c3, const vector<float> &c4);

bool GroupRotation(const vector<float> &axisA, const vector<float> &axisB, float angle, DOCKMol & mol, vector<bool> &flexibleatoms);
#endif
