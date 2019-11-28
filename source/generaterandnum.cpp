#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <string.h>

#include "generaterandnum.h"


double largeangle(mt19937 &gen)
{
    //random_device rd;  //Will be used to obtain a seed for the random number engine
    //mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    uniform_real_distribution<> dis(90.0, 180.0);
    return dis(gen);
}

double randgeneration(mt19937 &gen)
{
    uniform_real_distribution<> dis(-1.0, 1.0);
    return dis(gen);
};

double rand0to1(mt19937 &gen)
{
    uniform_real_distribution<> dis(0.0, 1.0);
    return dis(gen);
}

int randint(mt19937 &gen, int endnum)
{
    if(endnum >=1)
    {
       uniform_int_distribution <> dis(0,endnum-1);
       return dis(gen);
    }
    else
    return 0;	
    
}


