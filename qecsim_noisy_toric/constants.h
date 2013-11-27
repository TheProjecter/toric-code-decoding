#ifndef MYCONST_H
#define MYCONST_H


static const bool debug1 = false;
static const bool debug2 = false;
//check where at the program we are debug
static const bool debug3 = false;
//debug 4 makes the simulation 2Dimensional with one in time and one in space
static const bool debug4 = false;

//debug 5 prints progressbar
static const bool debug5 = false;



//Pmax is the maximum p determine_p_per_step is looking for
static const double PMAX = 0.5;
//PSTEP and DELTAPWUNSCH determine the precision of the
static const double PSTEP = 0.00001;
static const double DELTAPWUNSCH = 0.00005;


#endif
