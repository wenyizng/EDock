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

#if !defined(RAND_H)
    #define RAND_H 1


using namespace std;

double largeangle(mt19937 &gen);
double randgeneration(mt19937 &gen);
double rand0to1(mt19937 &gen);
int randint(mt19937 &gen, int endnum);
#endif 
