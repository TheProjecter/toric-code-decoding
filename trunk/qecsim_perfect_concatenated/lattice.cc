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

// zSize: 0 is the syndrome level, 1,2,3,4 are the qubit indices
    zSize = 5;

    cells = new Cell***[xSize];
    for (int i = 0; i < xSize; i++)
    {
        cells[i] = new Cell**[ySize];
        for (int j = 0; j < ySize; j++)
        {
            cells[i][j] = new Cell *[zSize];
            for ( int k = 0; k < zSize; k++)
            {
                //if k= 0     OR  i,j both odd OR    both even DO NOT generate a new cell (qubit)
                if ( k == 0 || (i%2 && j%2)  || (!(i%2) && !(j%2)))
                {
                    continue;
                }

                //otherwise, if only one of i,j is odd  and k > 0 generate a new cell (qubit)
                if(debug6)
                {
                    cout<< i<<j<<k<<" ";
                }
                cells[i][j][k] = new Cell (p);

            }//end of k
        }//end of j
    }// end of i


}


// Destructor of the lattice.
Lattice::~Lattice(void)
{

    for (int i = 0; i < xSize; i++)
    {
        for (int j = 0; j < ySize; j++)
        {
            for (int k = 0; k < zSize; k++)
            {
                if ((k ==0) || (i%2 && j%2) || (!(i%2) && !(j%2)))
                    continue;
                delete cells[i][j][k];
            }
            delete[] cells[i][j];
        }
        delete[] cells[i];
    }
    delete cells;

}




// Calculates Octaeder plaquette syndrome at (around) the specified coordinates.
bool Lattice::OCTsyndromeZZZZZZZZ (int xLoc, int yLoc)
{
//both i,j odd

    /* Debug
    cout<<" xLoc "<< xLoc <<"  xSize "<< xSize <<" yLoc  " <<yLoc <<" ySize  "<< ySize << " (xLoc % 2)  "<< (xLoc % 2)<<"   (yLoc % 2) "<< (yLoc % 2)<<endl;

    */
    assert (xLoc < xSize && yLoc < ySize && (xLoc % 2) && (yLoc % 2));

    bool result = false;

    Operator *P = new Operator(8);
    P->ops[0]  = Z;
    P->ops[1]  = Z;
    P->ops[2]  = Z;
    P->ops[3]  = Z;
    P->ops[4]  = Z;
    P->ops[5]  = Z;
    P->ops[6]  = Z;
    P->ops[7]  = Z;

    Operator *plaquette = new Operator(8);

//Remember QB are
//1--2
//|  |
//3--4
    //One step left (if i is line, y is column) to the square, get QB 2,4
    plaquette->ops[0]  = cells[xLoc][(yLoc-1+ySize)%ySize][2]->ErrCurrent->ops[0];
    plaquette->ops[1]  = cells[xLoc][(yLoc-1+ySize)%ySize][4]->ErrCurrent->ops[0];
    //One step down get QB 1,2
    plaquette->ops[2]  = cells[(xLoc+1)%xSize][yLoc][1]->ErrCurrent->ops[0];
    plaquette->ops[3]  = cells[(xLoc+1)%xSize][yLoc][2]->ErrCurrent->ops[0];
// One step right, get QB 1,3
    plaquette->ops[4]  = cells[xLoc][(yLoc+1)%ySize][1]->ErrCurrent->ops[0];
    plaquette->ops[5]  = cells[xLoc][(yLoc+1)%ySize][3]->ErrCurrent->ops[0];
//One step up get QB 3,4
    plaquette->ops[6]  = cells[(xLoc-1+xSize)%xSize][yLoc][3]->ErrCurrent->ops[0];
    plaquette->ops[7]  = cells[(xLoc-1+xSize)%xSize][yLoc][4]->ErrCurrent->ops[0];

    /* DEBUG
    //test where the plaquette tests the cells
    cout << "xLoc  " << xLoc << "    yLoc-1  " << yLoc-1 << endl;
    cout << "xLoc+1modxSiyze  " << (xLoc+1)%xSize << "    yLoc  " << yLoc << endl;
    cout << "xLoc  " << xLoc << "    yLoc+1modySize  " << (yLoc+1)%ySize << endl;
    cout << "xLoc-1  " << xLoc-1 <<"  yLoc  " << yLoc <<  endl;

    */

    result = !P->commute(*plaquette);

    delete P;
    delete plaquette;
    return result;
}

// Calculates SQUARE plaquette syndrome at (around) the specified coordinates.
bool Lattice::SQsyndromeZZZZ (int xLoc, int yLoc)
{
//exactly one of i,j odd and one even
    assert (xLoc < xSize && yLoc < ySize && ( ((xLoc % 2) && !(yLoc % 2)) || ( !(xLoc % 2) && (yLoc % 2))) );

    bool result = false;

    Operator *P = new Operator(4);
    P->ops[0]  = Z;
    P->ops[1]  = Z;
    P->ops[2]  = Z;
    P->ops[3]  = Z;


    Operator *plaquette = new Operator(4);

//Remember QB are
//1--2
//|ij|
//3--4
    plaquette->ops[0]  = cells[xLoc][yLoc][1]->ErrCurrent->ops[0];
    plaquette->ops[1]  = cells[xLoc][yLoc][2]->ErrCurrent->ops[0];
    plaquette->ops[2]  = cells[xLoc][yLoc][3]->ErrCurrent->ops[0];
    plaquette->ops[3]  = cells[xLoc][yLoc][4]->ErrCurrent->ops[0];



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
            signature.push_back(SQsyndromeZZZZ(i,j));
        }
    }

    return signature;
}


//  Transform X or Z errors at the specified location.
void Lattice::correct ( int xLoc, int yLoc, int QubitNr)
{
    if (debug4)
    {
        cout<< "xsize,ysixze,ysize "<<xSize<<ySize<<zSize<<endl;
        cout<< "xLoc,yLoc,QBNr "<<xLoc<<yLoc<<QubitNr<<endl;
    }
    assert (xLoc < xSize && yLoc < ySize && QubitNr < zSize);

//X*X=I
    pauli p = X;

    // Correct the error
    cells[xLoc][yLoc][QubitNr]->ErrCurrent->ops[0] = cells[xLoc][yLoc][QubitNr]->ErrCurrent->ops[0] * p;
    cells[xLoc][yLoc][QubitNr]->ErrGuess->ops[0] = cells[xLoc][yLoc][QubitNr]->ErrGuess->ops[0] * p;

    if (debug1)
    {
        cout << " -> correcting " << p << " at x=" << xLoc << " y=" << yLoc << " QB_Nr="<<QubitNr<<endl;
    }
    if (debug3)
    {
        cout << " -> correcting " << p << " at x=" << xLoc << " y=" << yLoc << " QB_Nr="<<QubitNr<<endl;
    }

}//end of correct()


// Calls the mathcing algorithm to find pairs of squares with nontrivial syndromes,
// and calls the correction subroutine.
void Lattice::correctMesh ()
{


    // Identify vertices that represent squares/octagons with non-trivial syndrome
    vector<Vertex*> vertices;

    for (int i = 0; i < xSize; i++)
    {
        for (int j = 0; j < ySize; j++)
        {
            bool syndrome = false;
// if i,j both odd do OCT-Z-syndrome (detects X Errors)

            if ( ((i%2) && (j%2)) )
            {
                syndrome = OCTsyndromeZZZZZZZZ(i,j);
            }
//if exactly  one of i and j is odd and one is even do SQ-Z-Syndrome

            if (  (i%2 && !(j%2))   ||   (!(i%2) && j%2) )
            {
                syndrome = SQsyndromeZZZZ(i,j);
            }

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
            //if both syndromes are both on the same odd column the distance in x increases by2
            //just a fancy specialty of this metrik (taxi-cab with holes)
            if ( (vertices[i]->y == vertices[j]->y) && (!(vertices[i]->y % 2)) )
            {
                deltaX = deltaX + 2;
            }

            int deltaY = min (abs (vertices[i]->y - vertices[j]->y), ySize - abs (vertices[i]->y - vertices[j]->y));
            //if the syndromes are both on the same odd line the distance in y increases by 2
            //just a fancy specialty of this metrik (taxi-cab with holes)
            if ( (vertices[i]->x == vertices[j]->x) && (!(vertices[i]->x % 2)) )
            {
                deltaY = deltaY + 2;
            }
            int dist   = (deltaX + deltaY);
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
            //cout<<"correctline" <<vertices[i]->x<<" "<< vertices[i]->y<<" to "<< vertices[j]->x<< " "<<vertices[j]->y<<endl;
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

//x y-->
//I
//v    EAST
//  00 01 02
//N 10 11 13 South
//  20 21 22
//  West

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

// go from x1|y1 -> x2|y2
// y2 = y1+(...)%ySize

//DIRECTION
    // Determine if we need to go up, down or neither
    vertical_direction v_dir = v_where_to_go(x1, x2);
//up,down, neither, wrong_direction

//HORIZONTAL Y-ACHSIS j
    // Determine if we need to start from sqr or oct
    celltype h_start_type = which_cell_type ( x1, y1);
//sqr,oct,wrong_type
    // Determine if we need to start from sqr or oct
    celltype h_end_type = which_cell_type ( x1, y2);

//VERTICAL X-ACHSIS i
// Determine if we need to start from sqr or oct
    celltype v_start_type = h_end_type;
    // Determine if we need to start from sqr or oct
    celltype v_end_type = which_cell_type ( x2, y2);

//Determine distance between syndromes

    int delta_x = delta_xy(x2,x1,xSize);
    int delta_y = delta_xy(y2,y1,ySize);








    /*old debugging
      if (debug1) cout << "Correcting line x1=" << x1 << " y1=" << y1 << " x2=" << x2 << " y2=" << y2 << " east=" << east << ":" << endl;
    */

    if (debug3)
    {
        cout <<endl<< "Correcting line x1=" << x1 << " y1=" << y1 << " x2=" << x2 << " y2=" << y2
             << " direction= ";
        if (v_dir==up)
        {
            cout<< " up "<<endl;
        }
        else if (v_dir==down)
        {
            cout<< " down "<<endl;
        }
        else if (v_dir==neither)
        {
            cout<< " neither "<<endl;
        }
        cout<< "Start Horizontal= ";
        if (h_start_type==sqr)
        {
            cout<< " square "<<endl;
        }
        else if (h_start_type==octo)
        {
            cout<< " 0ctogon "<<endl;
        }
        cout<< "End Horizontal= ";
        if (h_end_type==sqr)
        {
            cout<< " square "<<endl;
        }
        else if (h_end_type==octo)
        {
            cout<< " 0ctogon "<<endl;
        }
        cout<< "Start Vertical= ";
        if (v_start_type==sqr)
        {
            cout<< " square "<<endl;
        }
        else if (v_start_type==octo)
        {
            cout<< " 0ctogon "<<endl;
        }
        cout<< "End Vertical= ";
        if (v_end_type==sqr)
        {
            cout<< " square "<<endl;
        }
        else if (v_end_type==octo)
        {
            cout<< " 0ctogon "<<endl;
        }
        cout<< "Delta_y= "<<delta_xy(y2,y1,ySize)<<endl;
        cout<< "Delta_x= "<<delta_xy(x2,x1,xSize)<<endl;
    }
/////////////////////////////////////////////HORIZONTAL
    // correct in the horizontal direction first
//Qubits to correct in horizontal y Achsis
//y1 ++ to y2


// use correct( x1, y1, qubitcorr)

    //do we need to correct in horizontal direction?
    if (delta_y != 0)
    {

        if (debug4)
        {
            cout<<"delta_y != 0"<<endl;
        }
        //then start correction

        //exitvariable
        bool exit = false;
        //exit next time variable
        bool exitnextloop = false;


        //laufvariable
        int Y = y1;
        //countingvariable counting the amount of corrections
        int l = 1;

        //while exit is false
        while (!exit)
        {

            if (debug4)
            {
                cout<<"Y, y2: "<<Y<<" "<<y2<<endl;
            }


            //only if on an octagon
            if ( which_cell_type(x1,Y) == octo)
            {

                if (debug4)
                {
                    cout<<" is octagon "<<endl;
                }

                //initialize the stuff!
                bool first = false;
                bool last = false;
                correctiontype whatcorrection = wrong_correctiontype;
                correctiontypeQB whatcorrectionqubit = wrong_correctiontypeQB;


                //check if first
                if ( l == 1)
                {
                    first = true;
                }

                //check if last
                if ( delta_xy(Y,y2,ySize) <= 1 )
                {
                    last = true;
                }

                //Determine correctiontype
                if ( last && h_end_type == octo )
                {
                    whatcorrection = before;
                }
                else if (first && h_start_type)
                {
                    whatcorrection = next;
                }
                else
                {
                    whatcorrection = four;
                }

                //Determine correctiontypequbit

                if ( (whatcorrection == before) && ( v_dir == up) )
                {
                    whatcorrectionqubit = before_up;
                }
                else if ( (whatcorrection == before) && ( v_dir == down  || v_dir == neither) )
                {
                    whatcorrectionqubit = before_down;
                }
                else if ( (whatcorrection == next) && ( v_dir == up) )
                {
                    whatcorrectionqubit = next_up;
                }
                else if ( (whatcorrection == next) && ( v_dir == down  || v_dir == neither) )
                {
                    whatcorrectionqubit = next_down;
                }
                else if ( (whatcorrection == four) && ( v_dir == up) )
                {
                    whatcorrectionqubit = four_up;
                }
                else if ( (whatcorrection == four) && ( v_dir == down  || v_dir == neither) )
                {
                    whatcorrectionqubit = four_down;
                }
                else
                {
                    assert(0);
                }

                if(debug4)
                {
                    if (first)
                    {
                        cout<< " first ";
                    }
                    else
                    {
                        cout<<" not first ";
                    }
                    if (last)
                    {
                        cout<< " last ";
                    }
                    else
                    {
                        cout<<" not last ";
                    }

                    if (whatcorrectionqubit == before_up)
                    {
                        cout<< " before_up ";
                    }
                    else if (whatcorrectionqubit == before_down)
                    {
                        cout<< " before_down ";
                    }
                    else if(whatcorrectionqubit == next_up)
                    {
                        cout<< " next_up ";
                    }
                    else if(whatcorrectionqubit == next_down)
                    {
                        cout<< " next_down ";
                    }
                    else if(whatcorrectionqubit == four_up)
                    {
                        cout<< " four_up ";
                    }
                    else if(whatcorrectionqubit == four_down)
                    {
                        cout<< " four_down ";
                    }
                    else
                    {
                        assert(0);
                    }
                }




                switch (whatcorrectionqubit)
                {
                case before_up:
                    correct(x1,(Y+ySize-1)%ySize, 2);
                    Y = (Y+1)%ySize;
                    break;


                case before_down:
                    correct(x1,(Y+ySize-1)%ySize, 4);
                    Y = (Y+1)%ySize;
                    break;
                case next_up:
                    correct(x1,(Y+1)%ySize, 1);
                    Y = (Y+1)%ySize;
                    break;
                case next_down:
                    correct(x1,(Y+1)%ySize, 3);
                    Y = (Y+1)%ySize;
                    break;
                case four_down:
                    correct(x1,(Y+ySize-1)%ySize, 4);
                    correct((x1+1)%xSize,Y, 1);
                    correct((x1+1)%xSize,Y, 2);
                    correct(x1,(Y+1)%ySize, 3);
                    Y = (Y+1)%ySize;
                    break;
                case four_up:
                    correct(x1,(Y+ySize-1)%ySize, 2);
                    correct((x1+xSize-1)%ySize,Y, 3);
                    correct((x1+xSize-1)%ySize,Y, 4);
                    correct(x1,(Y+1)%ySize, 1);
                    Y = (Y+1)%ySize;
                    break;
                default:
                    assert(0);
                }

                //if we were at an octagon it shure must have corrected something
                // (otherwise assert(0);)-->so counter up!
                l++;

            }
            else if ( which_cell_type(x1,Y) == sqr )
            {

                Y = (Y+1)%ySize;

            }
            else
            {
                assert(0);
            }// end if octo or square

            // this makes the loop run a last time after Y hits y2
            // No, I do not want to talk about it!
            if (exitnextloop)
            {
                exit = true;
            }
            if (Y == y2)
            {
                exitnextloop = true;
            }




        } //endofwhileY
    }//endofifdeltay


/////////////////////////////////////////////VERTICAL
    // correct in the vertical direction
//Qubits to correct in vertical x Achsis
// go from x1|y2 to x2|y2
// x1 --> x2


// use correct( x1, y1, qubitcorr)

    //do we even need to correct in vertical direction?
    if (delta_x != 0)
    {
        //then start correction

        //exitvariable
        bool exit = false;
        //exit next time variable
        bool exitnextloop = false;

        //laufvariable
        int X = x1;
        //countingvariable counting the amount of corrections
        int l = 1;

        // while exit is false
        while (!exit)
        {

            //only if on an octagon
            if ( which_cell_type(X,y2) == octo)
            {


                //initialize that stuff
                bool first = false;
                bool last = false;
                correctiontype whatcorrection = wrong_correctiontype;
                correctiontypeQB whatcorrectionqubit = wrong_correctiontypeQB;



                //check if first
                if ( l == 1)
                {
                    first = true;
                }
                //check if last
                if ( delta_xy(X,x2,xSize) <= 1 )
                {
                    last = true;
                }

                //Determine correctiontype
                if ( last && v_end_type == octo )
                {
                    whatcorrection = before;
                }
                else if (first && v_start_type)
                {
                    whatcorrection = next;
                }
                else
                {
                    whatcorrection = four;
                }

                //Determine correctiontypequbit

                if ( (whatcorrection == before) && ( v_dir == up) )
                {
                    whatcorrectionqubit = before_up;
                }
                else if ( (whatcorrection == before) && ( v_dir == down ) )
                {
                    whatcorrectionqubit = before_down;
                }
                else if ( (whatcorrection == next) && ( v_dir == up) )
                {
                    whatcorrectionqubit = next_up;
                }
                else if ( (whatcorrection == next) && ( v_dir == down ) )
                {
                    whatcorrectionqubit = next_down;
                }
                else if ( (whatcorrection == four) && ( v_dir == up) )
                {
                    whatcorrectionqubit = four_up;
                }
                else if ( (whatcorrection == four) && ( v_dir == down  ) )
                {
                    whatcorrectionqubit = four_down;
                }
                else
                {
                    assert(0);
                }





                switch (whatcorrectionqubit)
                {
                case before_up:
                    correct((X+xSize+1)%xSize,y2, 1);
                    X = (X+xSize-1)%xSize;
                    break;
                case before_down:
                    correct((X+xSize-1)%xSize,y2, 3);
                    X = (X+1)%xSize;
                    break;
                case next_up:
                    correct((X+xSize-1)%xSize,y2, 3);
                    X = (X+xSize-1)%xSize;
                    break;
                case next_down:
                    correct((X+xSize+1)%xSize,y2, 1);
                    X = (X+1)%xSize;
                    break;
                case four_down:
                    correct((X+xSize-1)%xSize,y2, 3);
                    correct(X,(y2+ySize-1)%ySize, 2);
                    correct(X,(y2+ySize-1)%ySize, 4);
                    correct((X+xSize+1)%xSize,y2, 1);
                    X = (X+1)%xSize;
                    break;
                case four_up:
                    correct((X+xSize-1)%xSize,y2, 3);
                    correct(X,(y2+ySize-1)%ySize, 2);
                    correct(X,(y2+ySize-1)%ySize, 4);
                    correct((X+xSize+1)%xSize,y2, 1);
                    X = (X+xSize-1)%xSize;
                    break;
                default:
                    assert(0);
                }

                //if we were at an octagon it shure must have corrected something
                // (otherwise assert(0);)-->so counter up!
                l++;

            }
            else if ( which_cell_type(X,y2) == sqr )
            {

                switch (v_dir)
                {
                case up:
                    X = (X+xSize-1)%xSize;
                    break;
                case down:
                    X = (X+1)%xSize;
                    break;
                case neither:
                    assert(0);
                default:
                    assert(0);
                }


            }
            else
            {
                assert(0);
            }// end if octo or square


            // this makes the loop run a last time after X hits x2
            // No, I do not want to talk about it!
            if (exitnextloop)
            {
                exit = true;
            }
            if (X == x2)
            {
                exitnextloop = true;
            }




        } //endofwhileY
    }//endofifdeltax


}// end of function I guess...




// Determines if we need to go up, down or neither in the lattice
vertical_direction Lattice::v_where_to_go (int x1, int x2)
{

    vertical_direction v_dir = wrong_direction;

    if ( (x1 < x2 && x2 - x1 <  xSize - x2 + x1) ||
            (x2 < x1 && x1 - x2 >  xSize - x1 + x2) )
        v_dir = down;
    else if (x1 == x2)
        v_dir = neither;
    else
        v_dir = up;

    return v_dir;
}




// Determines if the celltype is square or an octagon
celltype Lattice::which_cell_type ( int x, int y)
{

    int i = x;
    int j = y;

    celltype type = wrong_type;

    if ( (i%2 && j%2)  || (!(i%2) && !(j%2)) )
    {
        type = octo;
    }
    else if ( (!(i%2) && j%2)  || (i%2 && !(j%2)) )
    {
        type = sqr;
    }
    else
    {
        type = wrong_type;
    }

    return type;

}//end of which_cell_type



//determines distance between coordinates in x direction
int Lattice::delta_xy(int xy1, int xy2, int Size)
{

    int d = min (abs(xy2-xy1), Size-abs(xy2-xy1) );

    return d;

}//endof delta_x



//returns true if there is no error syndrome anywhere on the lattice

bool Lattice::is_line()
{


    bool is_line_var = false;
    int syndromecounter = 0;


    // Identify vertices that represent squares/octagons with non-trivial syndrome


    for (int i = 0; i < xSize; i++)
    {
        for (int j = 0; j < ySize; j++)
        {
            bool syndrome = false;
// if i,j both odd do OCT-Z-syndrome (detects X Errors)

            if ( ((i%2) && (j%2)) )
            {
                syndrome = OCTsyndromeZZZZZZZZ(i,j);
            }
//if exactly  one of i and j is odd and one is even do SQ-Z-Syndrome

            if (  (i%2 && !(j%2))   ||   (!(i%2) && j%2) )
            {
                syndrome = SQsyndromeZZZZ(i,j);
            }


// if one of the above has a syndrome
            if (syndrome)
            {
                syndromecounter++;
            }
//otherwise go further
        }
    }

    if (syndromecounter == 2)
    {
        is_line_var = true;
    }


    return is_line_var;


}




//generates exactly 1 line of errors on the lattice

void Lattice::generateline()
{
    int x1 =  -1;
    int y1 = -1;

    int x2 =  -1;
    int y2 = -1;




    while(true)
    {

        x1 =  rand()%xSize;	//element [0,xSize]
        y1 = rand()%ySize;	//element [0,ySize]

        x2 =  rand()%xSize;	//element [0,xSize]
        y2 = rand()%ySize;	//element [0,ySize]

        //no syndrome at X Octagons erlaubt
        if ( (!(x1%2) && !(y1%2)) || (!(x2%2) && !(y2%2)))
        {
            continue;
        }



        //if we have generated two different points stop the generator
        if ((x1 != x2 || y1 != y2)    )
        {
            break;
        }

    }




    correctLine( x1, y1, x2, y2);


}//end of generateline



//returns true if there is no error syndrome anywhere on the lattice

bool Lattice::is_corrected()
{



    // Identify vertices that represent squares/octagons with non-trivial syndrome


    for (int i = 0; i < xSize; i++)
    {
        for (int j = 0; j < ySize; j++)
        {
            bool syndrome = false;
// if i,j both odd do OCT-Z-syndrome (detects X Errors)

            if ( ((i%2) && (j%2)) )
            {
                syndrome = OCTsyndromeZZZZZZZZ(i,j);
            }
//if exactly  one of i and j is odd and one is even do SQ-Z-Syndrome

            if (  (i%2 && !(j%2))   ||   (!(i%2) && j%2) )
            {
                syndrome = SQsyndromeZZZZ(i,j);
            }


// if one of the above has a syndrome
            if (syndrome)
            {
                return false;
            }
//otherwise go further
        }
    }



//if there was no syndrome anywhere
    return true;


}







// Generates logical operator  Z1, or Z2 on this lattice which do anticommute with nontrivil Xloops
Operator* Lattice::getLogical(string whichOp)
{

    //number of qubits
    Operator *myOp = new Operator (xSize*ySize/2*4);

    for (int i = 0, l = 0; i < xSize; i++)
    {
        for (int j = 0; j < ySize; j++)
        {
            for (int k = 0; k < zSize; k++)
            {

                if ( which_cell_type(i,j) == octo || (k == 0) )
                {
                    if(debug8)
                    {
                        cout<<"continue"<<endl;
                    }
                    continue;

                }

                if (whichOp == "Z1")
                {
                    //first line squares
                    if (((i == 0) && which_cell_type(i,j) == sqr && (k == 1)) || ((i == 0) && which_cell_type(i,j) == sqr && (k == 2)) )
                    {
                        if(debug8)
                        {
                            cout<<"Z"<<endl;
                        }
                        myOp->ops[l] = Z;
                    }
                    else
                    {
                        if(debug8)
                        {
                            cout<<"I"<<endl;
                        }
                        myOp->ops[l] = I;
                    }
                }
                else if (whichOp == "Z2")
                {
                    //first column squares
                    if ((j == 0 && which_cell_type(i,j) == sqr && k == 1) || (j == 0 && which_cell_type(i,j) == sqr && k == 3))
                    {
                        if(debug8)
                        {
                            cout<<"Z"<<endl;
                        }
                        myOp->ops[l] = Z;
                    }
                    else
                    {
                        if(debug8)
                        {
                            cout<<"I"<<endl;
                        }
                        myOp->ops[l] = I;
                    }
                }
                else
                {
                    assert (0);
                }
                l++;

            }//endofk
        }//endofj
    }//endofi

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
            for (int k = 0; k < zSize; k++)
            {
                if ((i%2 && j%2) || (!(i%2) && !(j%2)) || (k == 0))
                {
                    continue;
                }

                E.pushBack (cells[i][j][k]->ErrInit->ops[0]);
                if (p)
                {
                    pauli pErrGuessP = cells[i][j][k]->ErrGuess->ops[0] * p->ops[l++];
                    EGuessP.pushBack (pErrGuessP);
                }
                else
                {
                    EGuessP.pushBack (cells[i][j][k]->ErrGuess->ops[0]);
                }
            }
        }
    }
    Operator O = E * EGuessP;

    // Generate logical operators  Z1, Z2

    Operator *Z1 = getLogical("Z1");
    Operator *Z2 = getLogical("Z2");

    assert ( Z1->commute(*Z2));
    if (debug7)
    {
        cout<<"E0: ";
        E.printState();
        cout<<endl;
        cout<<"EG: ";
        EGuessP.printState();
        cout<<endl;
        cout<<"Ob: ";
        O.printState();
        cout<<endl;
        cout<<"Z1: ";
        Z1->printState();
        cout<<endl;

    }

    // Determine result of error correction
    bool result = true;
    if ( !O.commute(*Z1) || !O.commute(*Z2) )
    {
        result = false;
    }

    delete Z1;
    delete Z2;

    return result;

//der 2. Z muss noch weg dann
    return true;
}



// Print the current state of the lattice. The content of ErrCurrent for all qubits is printed.
void Lattice::printState(void)
{

    for (int i = 0; i < xSize; i++)
    {
        for (int j = 0; j < ySize; j++)
        {
            for (int k = 0; k < zSize; k++)
            {
                if ((k == 0) || (i%2 && j%2) || (!(i%2) && !(j%2)))
                    continue;
                cout << "Error at location X=" << i << " Y=" << j << " Z=" << k <<": ";
                cells[i][j][k]->ErrCurrent->printState();
                cout<< endl;
            }
        }
    }
}

// Print the current state of the QB lattice. The content of ErrCurrent for all qubits is printed.
void Lattice::printStateQB(void)
{
    cout<<"Print Lattice State"<<endl;

    for (int i = 0; i < xSize; i++)
    {
        for (int l = 0; l < 3; l++)
        {
            for (int j = 0; j < ySize; j++)
            {

                if ( (i%2 && j%2) || (!(i%2) && !(j%2)))
                {

                    if (l==1)
                    {
                        if(i%2 && j%2)
                        {
                            cout<< "    ";
                            // if there is syndrome print -
                            if(OCTsyndromeZZZZZZZZ(i,j))
                            {
                                cout<<"-";
                            }
                            else
                            {
                                cout<<"+";
                            }
                            cout<<"    ";
                            continue;
                        }
                        else
                        {
                            cout<< "         ";
                            continue;
                        }

                    }
                    else
                    {
                        cout<< "         ";
                        continue;
                    }
                }



                if (l==0)
                {
                    cout<<" ";
                    cells[i][j][1]->ErrCurrent->printState();
                    cout<<" --- ";
                    cells[i][j][2]->ErrCurrent->printState();
                    cout<<" ";


                }
                else if (l==1)
                {
                    cout<<"  | ";
                    //if there is syndrome print -
                    if(SQsyndromeZZZZ(i,j))
                    {
                        cout<<"-";
                    }
                    else
                    {
                        cout<<"+";
                    }
                    cout <<" |  ";
                }
                else if (l==2)
                {
                    cout<<" ";
                    cells[i][j][3]->ErrCurrent->printState();
                    cout<<" --- ";
                    cells[i][j][4]->ErrCurrent->printState();
                    cout<<" ";
                }










            }//end of j
            cout<<endl;
        }// end of l
        cout<<endl;



    }//end of i
}//end of print






