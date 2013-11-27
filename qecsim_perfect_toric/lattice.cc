#include "lattice.h"
#include "constants.h"


// Constructor generates a unit cell. In case of the surface code this is a single qubit.
Cell::Cell (double & p)
{
    ErrInit    = new Operator (1);
    ErrInit->generateRandom (p);
    // ErrCurrent is a copy of the random ErrInit
    ErrCurrent = new Operator (*ErrInit);
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
Vertex::Vertex (int newX, int newY)
{
    x  = newX;
    y  = newY;
}


// Constructor generates the lattice.
Lattice::Lattice(int newXSize, double p)
{

    xSize = 2*newXSize;
    ySize = 2*newXSize;

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
            cells[i][j] = new Cell (p);
        }
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
bool Lattice::syndromeZZZZ (int xLoc, int yLoc)
{

    assert (xLoc < xSize && yLoc < ySize && (xLoc % 2) && (yLoc % 2));

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
            signature.push_back(syndromeZZZZ(i,j));
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


// Calls the mathcing algorithm to find pairs of squares with nontrivial syndromes,
// and calls the correction subroutine.
void Lattice::correctMesh ()
{



    // Identify vertices that represent squares with non-trivial syndrome
    vector<Vertex*> vertices;

    for (int i = 0; i < xSize; i++)
    {
        for (int j = 0; j < ySize; j++)
        {
            bool syndrome = false;
// if i,j both odd and mesh id 1 do Zsyndrome (detects X Errors)
            if (i%2 && j%2 )
                syndrome = syndromeZZZZ(i,j);

            if (syndrome)
            {
                if (debug1) cout << "Syndrome:"  << " i=" << i << " j= " << j << endl;
                Vertex* v = new Vertex (i, j);
                vertices.push_back(v);
            }
        }
    }

    // Call minimum weight perfect matching
    struct PerfectMatching::Options options;
    struct GeomPerfectMatching::GPMOptions gpm_options;
    options.verbose = false;

    int numVert = vertices.size();
    int edges = (numVert*(numVert-1));
    PerfectMatching *pm = new PerfectMatching(numVert,edges);

    // Add all weighted edges (u,v) with weight w in complete graph G
    for (int i = 0; i < numVert; i++)
    {
        for (int j = 0; j < numVert; j++)
        {
            if (i == j)
                continue;
            int deltaX = min (abs (vertices[i]->x - vertices[j]->x), xSize - abs (vertices[i]->x - vertices[j]->x));
            int deltaY = min (abs (vertices[i]->y - vertices[j]->y), ySize - abs (vertices[i]->y - vertices[j]->y));
            int dist   = (deltaX + deltaY) / 2;
            pm->AddEdge(i,j,dist);
            if (debug1) cout << "Matching: i=" << i << " j= " << j << " dist=" << dist << endl;
        }
    }

    pm->options = options;
    pm->Solve();

    // Call the correction subroutine for all matched vertices
    for (int i = 0; i < numVert; i++ )
    {
        // i and j are matched
        int j = pm->GetMatch(i);
        if (i < j)
            correctLine (vertices[i]->x, vertices[i]->y, vertices[j]->x, vertices[j]->y);
    }

    for (int i = 0; i < numVert; i++)
    {
        delete vertices[i];
    }
    delete pm;

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
bool Lattice::success(Operator *p)
{

    // Calculate O = E * E_guessed
    Operator E(0), EGuessP(0);
    for (int i = 0, l = 0; i < xSize; i++)
    {
        for (int j = 0; j < ySize; j++)
        {
            if ((i%2 && j%2) || (!(i%2) && !(j%2)))
                continue;
            E.pushBack (cells[i][j]->ErrInit->ops[0]);
            if (p)
            {
                pauli pErrGuessP = cells[i][j]->ErrGuess->ops[0] * p->ops[l++];
                EGuessP.pushBack (pErrGuessP);
            }
            else
            {
                EGuessP.pushBack (cells[i][j]->ErrGuess->ops[0]);
            }
        }
    }
    Operator O = E * EGuessP;

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
                continue;
            cout << "Error at location X=" << i << " Y=" << j << ": ";
            cells[i][j]->ErrCurrent->printState();
        }
    }
}


