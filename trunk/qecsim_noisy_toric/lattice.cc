#include "lattice.h"
#include "constants.h"


// Constructor generates a unit cell. In case of the surface code this is a single qubit.

//Cell::Cell (double & p) {
Cell::Cell ()
{
    ErrInit    = new Operator (1);
// ErrInit->generateRandom (p);
    // ErrCurrent is a copy of the random ErrInit
    //ErrCurrent = new Operator (*ErrInit);
    ErrCurrent = new Operator (1);
    ErrGuess   = new Operator (1);

}


// Destructor for a unit cell.
Cell::~Cell (void)
{
    delete ErrInit;
    delete ErrCurrent;
    delete ErrGuess;
}


// Constructor of a vertex representing an operator in the transformed lattice.
Vertex::Vertex (int newX, int newY, int newT)
{
    x  = newX;
    y  = newY;
    t  = newT;
}


// Constructor generates the lattice.
//Lattice::Lattice(int newXSize, double p) {
Lattice::Lattice(int newXSize, int newT)
{

    if (debug3)
    {
        cout<<"Create Lattice"<<endl;
    }
    xSize = 2*newXSize;
    ySize = 2*newXSize;
    T	= newT;



    cells = new Cell**[xSize];
    for (int i = 0; i < xSize; i++)
    {
        cells[i] = new Cell*[ySize];
        for (int j = 0; j < ySize; j++)
        {
//if i,j both odd     OR    both even DO NOT generate a new cell (qubit)
            if ((i%2 && j%2) || (!(i%2) && !(j%2)))
            {
                continue;
            }
//otherwise, if only one of i,j is odd generate a new cell (qubit)
            cells[i][j] = new Cell ();
        }
    }

    if (debug3)
    {
        cout<<"Lattice created"<<endl;
    }
}


// Destructor of the lattice.
Lattice::~Lattice(void)
{

    for (int i = 0; i < xSize; i++)
    {
        for (int j = 0; j < ySize; j++)
        {
            if ((i%2 && j%2) || (!(i%2) && !(j%2)))
                continue;
            delete cells[i][j];
        }
        delete[] cells[i];
    }
    delete cells;

}




// Calculates plaquette syndrome at (around) the specified coordinates.
bool Lattice::syndromeZZZZ (int xLoc, int yLoc, double q_err)
{

    assert (xLoc < xSize && yLoc < ySize && (xLoc % 2) && (yLoc % 2));

    //q could be modeled here as well
    double q = q_err;

    bool result = false;

    Operator *P = new Operator(4);
    P->ops[0]  = Z;
    P->ops[1]  = Z;
    P->ops[2]  = Z;
    P->ops[3]  = Z;

    Operator *plaquette = new Operator(4);
    plaquette->ops[0]  = cells[xLoc][yLoc-1]->ErrCurrent->ops[0];
    plaquette->ops[1]  = cells[(xLoc+1)%xSize][yLoc]->ErrCurrent->ops[0];
    plaquette->ops[2]  = cells[xLoc][(yLoc+1)%ySize]->ErrCurrent->ops[0];
    plaquette->ops[3]  = cells[xLoc-1][yLoc]->ErrCurrent->ops[0];

    result = !P->commute(*plaquette);




    //generate noisy syndrome measurement with propapility q
    double r = (double)rand() / ((double)RAND_MAX + 1);
    //if it shall have noise, then flip it
    if (r < q)
    {
        result = !result;
    }


    delete P;
    delete plaquette;
    return result;
}


// Calculates a signature which is a concatentation of all site and plaquette syndromes in the lattice.
vector<bool> Lattice::getSignature (void)
{

    // The signature starts out empty
    vector<bool> signature;



    // Add ZZZZ syndromes
    for (int i = 1; i < xSize; i+=2)
    {
        for (int j = 1; j < ySize; j+=2)
        {
            //signature.push_back(syndromeZZZZ(i,j));
        }
    }

    return signature;
}


//  Transform X or Z errors at the specified location.
void Lattice::correct ( int xLoc, int yLoc)
{

    assert (xLoc < xSize && yLoc < ySize);
//X*X=I
    pauli p = X;

    // Correct the error
    cells[xLoc][yLoc]->ErrCurrent->ops[0] = cells[xLoc][yLoc]->ErrCurrent->ops[0] * p;
    cells[xLoc][yLoc]->ErrGuess->ops[0] = cells[xLoc][yLoc]->ErrGuess->ops[0] * p;

    if (debug1) cout << " -> correcting " << p << " at x=" << xLoc << " y=" << yLoc << endl;
}


//evtl void Lattice::evolute(double & p_err){
void Lattice::evolute(double p_err)
{
    if (debug3)
    {
        cout<<"Evolute "<<endl;
    }
    //qubit error
    double p = p_err;

    for (int i = 0; i < xSize; i++)
    {
        for (int j = 0; j < ySize; j++)
        {

            if (debug4)
            {

                //only second line nr 1 gets errors
                if (i==1 && !(j%2))  //j=1,3,5is syndrome,0,2,4 is qubit
                {
                    cells[i][j]->ErrInit->generate_random_evolution(p);
                    //and copy into Current
                    cells[i][j]->ErrCurrent->ops[0] = cells[i][j]->ErrInit->ops[0];
                }

            }
            else
            {

                //if i,j both odd     OR    both even do nothing
                if ((i%2 && j%2) || (!(i%2) && !(j%2)))
                {
                    continue;
                }
                //otherwise, throw error at qubit
                cells[i][j]->ErrInit->generate_random_evolution(p);
                //and copy into Current
                cells[i][j]->ErrCurrent->ops[0] = cells[i][j]->ErrInit->ops[0];

                //if that doesnt work, we need a copy method for operators (seems to work:-)

            }



        }
    }



}






// Calls the mathcing algorithm to find pairs of squares with nontrivial syndromes,
// and calls the correction subroutine.
void Lattice::evolute_and_correct(double p_err, double q_err)
{

    if (debug3)
    {
        cout<<"Evolute and Correct "<<endl;
    }

    //proper algorithm for q
    double p = p_err;
    double q = q_err;


    /*
    find proper algorithm for q, (is maybe done in simulation.cc)

    create history

    loop:
    evolute,measure,save

    create vertices

    match vertices

    correct


    */

    //Create 3d
    bool*** history;


    history = new bool** [xSize];

    //fill history with false;
    for (int i = 0; i < xSize; i++)
    {
        history[i] = new bool *[ySize];
        for (int j = 0; j < ySize; j++)
        {
            history[i][j] = new bool [T+1];
            for (int t = 0; t < T+1; t++)
            {
                history[i][j][t] = false;
            }
        }
    }


    //now evolute the lattice, measure after each time step and save the syndromes in the history.


    if (debug4)
    {
        cout<<"Evolution"<<endl;
    }

    for (int t = 0; t < T; t++)
    {
        //evolute

        this->evolute(p);

        if (debug4)
        {
            cout<<"T="<<t<<": ";
            print2DState();
        }


        //measure and save

        for (int i = 0; i < xSize; i++)
        {
            for (int j = 0; j < ySize; j++)
            {


                // if i,j both odd and mesh id 1 do Zsyndrome (detects X Errors)
                if (i%2 && j%2 && syndromeZZZZ(i,j,q))
                {

                    //save in history
                    history[i][j][t] = true;
                    if (debug3)
                    {
                        cout<<"history i,j,t = "<<i<<" "<<j<<" "<<t<<" = "<<history[i][j][t]<<endl;
                    }
                }


            }//end of j
        }//end of i
    }//end of for t evolute and measure

    if (debug3 || debug4)
    {
        cout<<"Finished Measurement and Saving History"<<endl;
    }

    //Print history
    if (debug4)
    {
        for (int t = 0; t < T+1; t++)
        {

            cout<<"T="<<t<<": ";
            for (int i = 0; i < xSize; i++)
            {
                for (int j = 0; j < ySize; j++)
                {


                    if (i==1 && !(j%2))  //j=1,3,5is syndrome,0,2,4 is qubit
                    {
                        cout<< " ";
                    }
                    else if (i==1 && j%2)
                    {
                        cout<<history[i][j][t];
                    }

                }
                if (i==1)
                {
                    cout<<endl;
                }
            }//end of i

        }//endof t
        cout<<endl;
    }



    // Identify vertices that represent squares with change in syndrome
    vector<Vertex*> vertices;
    if (debug4)
    {
        cout<<"Vertices"<<endl;
    }

    for (int t = 0; t < T+1; t++)
    {
        if (debug4)
        {
            cout<<"T="<<t<<": ";
        }
        for (int i = 0; i < xSize; i++)
        {
            for (int j = 0; j < ySize; j++)
            {

                bool is_vertex = false;
                // if i,j both odd check if vertex
                if (i%2 && j%2 )
                {
                    /*
                    in the first and last slice, every syndrome is a vertex
                    otherwise, check the current t and the one before
                    if the syndrome changed, it is considered a vertex
                    */
                    if((t==0 && history[i][j][t]) || (t!=0 && history[i][j][t]!=history[i][j][t-1]))
                    {
                        is_vertex = true;
                    }//end if is vertex
                }//end of if Plaquette (i,j)

                //Print vertex history
                if (debug4)
                {


                    if (i==1 && !(j%2))  //j=1,3,5is syndrome,0,2,4 is qubit
                    {
                        cout<< " ";
                    }
                    else if (i==1 && j%2)
                    {
                        if(is_vertex)
                        {
                            cout<<vertices.size();
                        }
                        else
                        {
                            cout<<".";
                        }

                    }


                }






                //if it is a vertex, save it
                if (is_vertex)
                {
                    if (debug3) cout << "Is_Vertex at: "  << " i=" << i << " j= " << j << " t= " << t<< endl;
                    Vertex* v = new Vertex (i, j, t);
                    vertices.push_back(v);

                }

            }//end of j

        }//end of i
        if (debug4)
        {
            cout<<endl;
        }
    }//end of t



    if (debug3)
    {
        cout<<"Finished Creating Vertices "<<endl;
    }




    int numVert = vertices.size();
    int edges = 0;

    double** edgeLengths;
    edgeLengths = new double*[numVert];
    for (int i = 0; i < numVert; i++)
    {
        edgeLengths[i] = new double[numVert];
        for (int j = 0; j < numVert; j++)
        {
            edgeLengths[i][j] = -1.0;
        }
    }



    // Add all weighted edges (u,v) with weight w in complete graph G
    for (int i = 0; i < numVert; i++)
    {
        for (int j = 0; j < numVert; j++)
        {
            if (i == j)
            {
                continue;
            }





//this is the distance in x y t
            double deltaX = min (abs (vertices[i]->x - vertices[j]->x), xSize - abs (vertices[i]->x - vertices[j]->x));
            double deltaY = min (abs (vertices[i]->y - vertices[j]->y), ySize - abs (vertices[i]->y - vertices[j]->y));
            double deltaT = abs (vertices[i]->t - vertices[j]->t);


            /*
            //here, some fancy algorithm for the distance is needed (later on for q > 0)
            if(distance_tobi){
            // Tobias' algorithm


            	//First, make shure those guys are connected through a finite distance
            	bool is_infinite = false;
            	//if prob is 0 but there is a distance they shouldn be connected
            	if ((p==0 && (deltaY>0 || deltaX>0)) ||  (q==0 && deltaT>0)){
            		is_infinite = true;
            	}

            	//so if the distance is infinite, dont add the edge
            	if(is_infinite){continue;}


            	//otherwise, weight the edge and add it:



            	//max delta is the most unprobable distance on the lattice
            	//XSize/2 is the real lattice, /2 again, because it is a toric code and max dist is half the lattice
            	double min_p_factor = pow(p,(xSize/4));
            	double max_delta = 1/min_p_factor;




            	if(p>0) {
            		double X_p_factor = p;//pow(p,deltaX);
            		double Y_p_factor = p;//pow(p,deltaY);

            		//deltaX *=(1/X_p_factor);
            		//deltaY *=(1/Y_p_factor);
            		deltaX *=(1-X_p_factor);
            		deltaY *=(1-Y_p_factor);
            	}

            	if(q>0) {
            		double T_q_factor = q;//pow(q,deltaT);
            		//deltaT *=(1/T_q_factor);
            		deltaT *=(1-T_q_factor);
            	}

            //ENd of Tobias' algorithm
            */


//SUCHARAS' algorithm #######################################

// Calculate relative edge weights of time and space edges
            double ws = 1000;
            double wt = 1000;
            //double wf = 1000;

            if (p!=0) ws = log ((1 - p) / p);
            if (q!=0) wt = log ((1 - q) / q);
            //if (p!=0) wf = log ((1 - p) / p);

            /* if (virtualRound && (vertices[i]->t == virtualRound || vertices[j]->t == virtualRound)) {
             	    deltaX *= wf;
               deltaY *= wf;
             } else if (virtualRound) {
               deltaX *= ws;
               deltaY *= ws;
             }*/

            deltaX *= ws;
            deltaY *= ws;
            deltaT *= wt;


//END SUCHARAS' algorithm ############################################


            //half the distance because of the lattice shape vs. the array shape
            double dist   = ((deltaX + deltaY)/2)  + (double) deltaT;
            edgeLengths[i][j] = dist;
            edges++;
            if (debug1) cout << "Matching: i=" << i << " j= " << j << " dist=" << dist << endl;
        }
    }



// Call minimum weight perfect matching
    struct PerfectMatching::Options options;
    struct GeomPerfectMatching::GPMOptions gpm_options;
    options.verbose = false;
    PerfectMatching *pm = new PerfectMatching(numVert,edges);

    for (int i = 0; i < numVert; i++)
    {
        for (int j = 0; j < numVert; j++)
        {
            double dist = edgeLengths[i][j];
            if (dist >= 0)
                pm->AddEdge(i,j,dist);
        }
    }



    pm->options = options;
    pm->Solve();

    if (debug4)
    {
        for (int i = 0; i < numVert; i++ )
        {
            int j = pm->GetMatch(i);
            if (i<j)
            {
                cout<<"Matched: "<<i<<" and "<<j<<endl;
            }
        }
    }

    // Call the correction subroutine for all matched vertices
    for (int i = 0; i < numVert; i++ )
    {
        // i and j are matched
        int j = pm->GetMatch(i);
        //only correct:
        //-once
        // -those vertices whose connections have more dimensions than time == d_x OR d_y >0
        //-the last history slice is special:
        // if there is a vertex, it is either in time and does not need
        //correction OR it is in space, because the last slice is without
        //error by default SO: only correct which matching partner  vertex is not in the last slice == correct in space AND time


        if ( !(vertices[i]->t == T && vertices[j]->t == T) && (i < j) && (j < numVert) && ( (vertices[i]->x != vertices[j]->x) || (vertices[i]->y != vertices[j]->y) )   )
        {
            if (debug4)
            {
                cout<<"Corrected: "<<i<<" to "<<j<<endl;
            }
            correctLine (vertices[i]->x, vertices[i]->y, vertices[j]->x, vertices[j]->y);


        }
    }

    for (int i = 0; i < numVert; i++)
    {
        delete vertices[i];
    }

    delete pm;

    for (int i = 0; i < numVert; i++)
    {
        delete edgeLengths[i];
    }
    delete edgeLengths;


// Delete syndrome history
    for (int i = 0; i < xSize; i++)
    {
        for (int j = 0; j < ySize; j++)
        {
            delete history[i][j];
        }
        delete history[i];
    }
    delete history;

    if (debug3)
    {
        cout<<"Finished Correction "<<endl;
    }
}





// Does one round of perfect syndrome measurement and corrects the lattice
// Without it, it might not be part of the stabilizer
//(errors could still be in it)
void Lattice::perfect_correction(void)
{

    if (debug3)
    {
        cout<<"Evolute and Correct "<<endl;
    }

    double q = 0;



    //Create 1d
    bool*** history;


    history = new bool** [xSize];

    //fill history with false;
    for (int i = 0; i < xSize; i++)
    {
        history[i] = new bool *[ySize];
        for (int j = 0; j < ySize; j++)
        {
            history[i][j] = new bool [T+1];
            for (int t = 0; t < 1; t++)
            {
                history[i][j][t] = false;
            }
        }
    }


    //measure and save
    for (int t = 0; t < 1; t++)
    {
        for (int i = 0; i < xSize; i++)
        {
            for (int j = 0; j < ySize; j++)
            {


                // if i,j both odd and mesh id 1 do Zsyndrome (detects X Errors)
                if (i%2 && j%2 && syndromeZZZZ(i,j,q))
                {

                    //save in history
                    history[i][j][t] = true;
                }

            }//end of j
        }//end of i
    }//end of t




    // Identify vertices that represent squares with change in syndrome
    vector<Vertex*> vertices;
    if (debug4)
    {
        cout<<"Vertices"<<endl;
    }

    for (int t = 0; t < 1; t++)
    {
        for (int i = 0; i < xSize; i++)
        {
            for (int j = 0; j < ySize; j++)
            {

                bool is_vertex = false;
                // if i,j both odd check if vertex
                if (i%2 && j%2 && history[i][j][t])
                {
                    is_vertex = true;

                }//end of is vertex






                //if it is a vertex, save it
                if (is_vertex)
                {
                    if (debug3) cout << "Is_Vertex at: "  << " i=" << i << " j= " << j << " t= " << t<< endl;
                    Vertex* v = new Vertex (i, j, t);
                    vertices.push_back(v);

                }

            }//end of j

        }//end of i
    }//end of t







    int numVert = vertices.size();
    int edges = 0;

    double** edgeLengths;
    edgeLengths = new double*[numVert];
    for (int i = 0; i < numVert; i++)
    {
        edgeLengths[i] = new double[numVert];
        for (int j = 0; j < numVert; j++)
        {
            edgeLengths[i][j] = -1.0;
        }
    }



    // Add all weighted edges (u,v) with weight w in complete graph G
    for (int i = 0; i < numVert; i++)
    {
        for (int j = 0; j < numVert; j++)
        {
            if (i == j)
            {
                continue;
            }



//this is the distance in x y t
            double deltaX = min (abs (vertices[i]->x - vertices[j]->x), xSize - abs (vertices[i]->x - vertices[j]->x));
            double deltaY = min (abs (vertices[i]->y - vertices[j]->y), ySize - abs (vertices[i]->y - vertices[j]->y));
            double deltaT = abs (vertices[i]->t - vertices[j]->t);

            //half the distance because of the lattice shape vs. the array shape
            double dist   = ((deltaX + deltaY)/2)  + (double) deltaT;
            edgeLengths[i][j] = dist;
            edges++;
            if (debug1) cout << "Matching: i=" << i << " j= " << j << " dist=" << dist << endl;
        }
    }



// Call minimum weight perfect matching
    struct PerfectMatching::Options options;
    struct GeomPerfectMatching::GPMOptions gpm_options;
    options.verbose = false;
    PerfectMatching *pm = new PerfectMatching(numVert,edges);

//cout<<"numVert: "<<numVert<<endl;

    for (int i = 0; i < numVert; i++)
    {
        for (int j = 0; j < numVert; j++)
        {
            double dist = edgeLengths[i][j];
            if (dist >= 0)
                pm->AddEdge(i,j,dist);
        }
    }



    pm->options = options;
    pm->Solve();

    if (debug4)
    {
        for (int i = 0; i < numVert; i++ )
        {
            int j = pm->GetMatch(i);
            if (i<j)
            {
                cout<<"Matched: "<<i<<" and "<<j<<endl;
            }
        }
    }

    // Call the correction subroutine for all matched vertices
    for (int i = 0; i < numVert; i++ )
    {
        // i and j are matched
        int j = pm->GetMatch(i);

        if (  (i < j) && (j < numVert) && ( (vertices[i]->x != vertices[j]->x) || (vertices[i]->y != vertices[j]->y) )   )
        {
            if (debug4)
            {
                cout<<"Corrected: "<<i<<" to "<<j<<endl;
            }
            correctLine (vertices[i]->x, vertices[i]->y, vertices[j]->x, vertices[j]->y);

        }
    }

    for (int i = 0; i < numVert; i++)
    {
        delete vertices[i];
    }

    delete pm;

    for (int i = 0; i < numVert; i++)
    {
        delete edgeLengths[i];
    }
    delete edgeLengths;


// Delete syndrome history
    for (int i = 0; i < xSize; i++)
    {
        for (int j = 0; j < ySize; j++)
        {
            delete history[i][j];
        }
        delete history[i];
    }
    delete history;

    if (debug3)
    {
        cout<<"Finished Correction "<<endl;
    }
}















// Given the locations of two syndromes with a non-trivial syndrome, correct errors on a line
// connecting the two syndromes.
void Lattice::correctLine ( int x1Loc, int y1Loc, int x2Loc, int y2Loc)
{

    // Initialize x1, y1, x2, y2 so that we can correct errors on a line in
    // the SE or SW direction coming from coordinates 1 to coordinates 2
    int x1, y1, x2, y2;
    if ( (y1Loc < y2Loc && y2Loc - y1Loc < ySize - y2Loc + y1Loc) ||
            (y2Loc < y1Loc && y1Loc - y2Loc > ySize - y1Loc + y2Loc) )
    {
        x1 = x1Loc;
        y1 = y1Loc;
        x2 = x2Loc;
        y2 = y2Loc;
    }
    else
    {
        x1 = x2Loc;
        y1 = y2Loc;
        x2 = x1Loc;
        y2 = y1Loc;
    }

    // Determine if we need to go to SE or SW
    bool east;
    if ( (x1 < x2 && x2 - x1 <  xSize - x2 + x1) ||
            (x2 < x1 && x1 - x2 >  xSize - x1 + x2) )
        east = true;
    else
        east = false;

    if (debug1) cout << "Correcting line x1=" << x1 << " y1=" << y1 << " x2=" << x2 << " y2=" << y2 << " east=" << east << ":" << endl;

    // correct in the southbound direction first
    while (y1 != y2)
    {
        y1 = (y1 + 2) % (ySize);
        correct ( x1, (y1-1+ySize) % ySize);
    }

    // if correction in SE direction
    if (east)
    {
        while (x1 != x2)
        {
            x1 = (x1 + 2) % (xSize);
            correct ( (x1-1+xSize) % xSize, y1);
        }
        // if correction in SW direction
    }
    else
    {
        while (x1 != x2)
        {
            x1 = (x1 -2 + xSize) % (xSize);
            correct ( (x1+1) % xSize, y1);
        }
    }
}


// Generates logical operator  Z1, or Z2 on this lattice which do anticommute with nontrivil Xloops
Operator* Lattice::getLogical(string whichOp)
{


    Operator *myOp = new Operator (xSize*ySize/2);

    for (int i = 0, k = 0; i < xSize; i++)
    {
        for (int j = 0; j < ySize; j++)
        {
            if ((i%2 && j%2) || (!(i%2) && !(j%2)))
                continue;
            if (whichOp == "Z1")
            {
                if (i == 0)
                {
                    myOp->ops[k] = Z;
                }
            }
            else if (whichOp == "Z2")
            {
                if (j == 0)
                {
                    myOp->ops[k] = Z;
                }
            }
            else
            {
                assert (0);
            }
            k++;
        }
    }

    return myOp;
}


// Examimes ErrInit and ErrGuess to determine if error correction succeeds.
// Multiplies ErrGuess with operator p.
bool Lattice::success(void)
{

    // Calculate O = E * E_guessed
    Operator E(0), EGuessP(0);

    for (int i = 0; i < xSize; i++)
    {
        for (int j = 0; j < ySize; j++)
        {
            if ((i%2 && j%2) || (!(i%2) && !(j%2)))
                continue;

            E.pushBack (cells[i][j]->ErrInit->ops[0]);
            EGuessP.pushBack (cells[i][j]->ErrGuess->ops[0]);
            //ECurrent.pushBack (cells[i][j]->ErrCurrent->ops[0]);

        }
    }
    Operator O = E * EGuessP;
    /*
    	Debug

      Operator  ECurrent(0);
      //O and ECurrent should be the same
    	cout<<endl<<"E: ";
    	E.printState();
    	cout<<endl<<"G: ";
    	EGuessP.printState();
    	cout<<endl<<"O: ";
    	O.printState();
    	cout<<endl<<"C: ";
    	ECurrent.printState();
    	cout<<endl<<endl;
    */





    // Generate logical operators  Z1, Z2

    Operator *Z1 = getLogical("Z1");
    Operator *Z2 = getLogical("Z2");

    assert ( Z1->commute(*Z2));

    // Determine result of error correction
    bool result = true;
    //if ( !O.commute(*Z1) || !O.commute(*Z2))
    if ( !O.commute(*Z1) )
        result = false;

    delete Z1;
    delete Z2;
    if (debug3)
    {
        cout<<"Return Result Success =  "<<result<<endl;
    }
    return result;
}


// Print the current state of the lattice. The content of ErrCurrent for all qubits is printed.
void Lattice::printState(void)
{

    for (int i = 0; i < xSize; i++)
    {
        for (int j = 0; j < ySize; j++)
        {



            if ((i%2 && j%2) || (!(i%2) && !(j%2)))
            {
                continue;
            }


            cells[i][j]->ErrCurrent->printState();



        }
    }
}

// Print the current 2Dstate of the lattice. The content of ErrCurrent and ZSyndrome for all qubits in i==1 is printed.
void Lattice::print2DState(void)
{

    for (int i = 0; i < xSize; i++)
    {
        for (int j = 0; j < ySize; j++)
        {



            if (i==1 && !(j%2))  //j=1,3,5is syndrome,0,2,4 is qubit
            {
                cells[i][j]->ErrCurrent->printState();
            }
            else if (i==1 && j%2)
            {
                if (syndromeZZZZ (i, j, 0))
                {
                    cout<<"-";
                }
                else
                {
                    cout<<"+";
                }
            }
        }
    }
    cout<<endl;
}



// Print the current 2Dstate of the Initialized Error of the lattice.
void Lattice::print2DInit(void)
{

    for (int i = 0; i < xSize; i++)
    {
        for (int j = 0; j < ySize; j++)
        {



            if (i==1 && !(j%2))  //j=1,3,5is syndrome,0,2,4 is qubit
            {
                cells[i][j]->ErrInit->printState();
            }
            else if (i==1 && j%2)
            {
                cout<<"O";
            }
        }
    }
    cout<<endl;
}

//end of file lattice.cc
