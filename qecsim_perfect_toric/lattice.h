#ifndef LATTICE_H
#define LATTICE_H

#define MAX_NEWRAND 65535

#include <stdlib.h>
#include <string>
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <assert.h>
#include <vector>
#include <algorithm>
#include <inttypes.h>
#include "operator.h"
#include "blossom5-v2.02.src/PerfectMatching.h"
#include "blossom5-v2.02.src/GEOM/GeomPerfectMatching.h"


using namespace std;


class Cell
{
public:
    Operator *ErrInit;
    Operator *ErrCurrent;
    Operator *ErrGuess;

    Cell (double & p);
    ~Cell (void);
};


class Vertex
{
public:
    int x;
    int y;

    Vertex (int newX, int newY);
    ~Vertex (void) { }
};


class Lattice
{
public:
    int xSize;
    int ySize;
    Cell ***cells;


    Lattice (int newXSize, double p);
    ~Lattice (void);
    bool syndromeXXXX (int xLoc, int yLoc);
    bool syndromeZZZZ (int xLoc, int yLoc);
    vector<bool> getSignature (void);
    void correct (int xLoc, int yLoc);
    void correctMesh ();
    void correctLine ( int x1Loc, int y1Loc, int x2Loc, int y2Loc);
    Operator* getLogical (string whichOp);
    bool success (Operator *p);
    void printState (void);
    // New random generator
    unsigned int randomGeneratorState;
    unsigned int NewRand(void);
    void NewSrand(unsigned int seed);

};


#endif
