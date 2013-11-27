#ifndef SIMULATION_H
#define SIMULATION_H

#include <stdlib.h>
#include <string>
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <assert.h>
#include <vector>
#include <set>
#include <map>
#include <fstream>
#include <climits>
#include <bitset>
#include <time.h>
#include <inttypes.h>
#include "operator.h"
#include "lattice.h"


using namespace std;


class Simulation
{
public:
    vector<int> XSize;
    int iterations;
    double pMin, pMax, pStep;

    Simulation (void);
    ~Simulation (void) { }
    void run (void);
};


#endif
