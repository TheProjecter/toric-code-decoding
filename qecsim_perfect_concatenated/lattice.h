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

// which way to go for the correction line
//typedef enabled by default, not to yield compiler warnings
//typedef enum vertical_direction {up, down, neither, wrong_direction};
enum vertical_direction {up, down, neither, wrong_direction};
// is my cell a square or an octagon?
//typedef enabled by default
//typedef enum celltype {sqr, octo, wrong_type};
enum celltype {sqr, octo, wrong_type};
//typedef enabled by default
//what kind of correction shape is it?
//typedef enum correctiontype {four, before, next};
enum correctiontype {four, before, next,wrong_correctiontype};
//what kind of correction Qubit is it?
//typedef enum correctiontype {four, before, next};
enum correctiontypeQB {four_up, four_down, before_up,before_down, next_up, next_down,wrong_correctiontypeQB};





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
    int zSize;
    Cell ****cells;


    Lattice (int newXSize, double p);
    ~Lattice (void);
    bool OCTsyndromeZZZZZZZZ (int xLoc, int yLoc);
    bool SQsyndromeZZZZ (int xLoc, int yLoc);
    vector<bool> getSignature (void);
    void correct (int xLoc, int yLoc , int QubitNr);
    void correctMesh ();
    void correctLine ( int x1Loc, int y1Loc, int x2Loc, int y2Loc);
    Operator* getLogical (string whichOp);
    vertical_direction v_where_to_go(int x1, int x2);
    int delta_xy(int x1, int x2, int Size);
    celltype which_cell_type (int x, int y);
    bool is_line();
    bool is_corrected();
    void generateline();
    bool success (Operator *p);
    void printState (void);
    void printStateQB(void);

};


#endif
