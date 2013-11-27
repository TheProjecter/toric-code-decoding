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

    //Cell (double & p);
    Cell ();
    ~Cell (void);
};


class Vertex
{
public:
    int x;
    int y;
    int t;

    Vertex (int newX, int newY, int newT);
    ~Vertex (void) { }
};


class Lattice
{
public:
    int xSize;
    int ySize;
    int T;
    Cell ***cells;


    Lattice (int newXSize, int newT);
    ~Lattice (void);
    bool syndromeZZZZ (int xLoc, int yLoc, double q_err);
    vector<bool> getSignature (void);
    void correct (int xLoc, int yLoc);
    void evolute_and_correct (double p_err, double q_err);
    void evolute(double p_err);
    void perfect_correction(void);
    void correctLine ( int x1Loc, int y1Loc, int x2Loc, int y2Loc);
    Operator* getLogical (string whichOp);
    bool success (void);
    void printState (void);
    void print2DState(void);
    void print2DInit(void);

};


#endif
