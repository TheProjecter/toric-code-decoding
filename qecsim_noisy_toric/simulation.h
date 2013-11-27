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
#include <iomanip>      // std::setprecision
#include "operator.h"
#include "lattice.h"



using namespace std;


class Simulation
{
public:
    vector<int> XSize;
    int iterations;
    double pMin, pMax, pStep, q_factor;

    Simulation (void);
    ~Simulation (void) { }
    void run (void);

//for the error in time
    double faculty (double n);
    double ncr (double n, double k);
    double determine_p_per_step (double p_G, int T_in);
    double determine_p_effective (double p_G, int T_in);



};


#endif
