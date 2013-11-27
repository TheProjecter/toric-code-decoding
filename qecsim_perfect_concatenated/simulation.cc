#include "simulation.h"
#include "constants.h"
#include <omp.h>

// Constructor reads parameters from a config file.
Simulation::Simulation (void)
{

    string tmp;
    string item;
    int pos;
    ifstream f_conf("in/parameters.txt");

    getline(f_conf,tmp);
    pos = tmp.find(" = ");
    tmp = tmp.substr (pos + 3);
    stringstream ss(tmp);
    while(ss >> item)
    {
        XSize.push_back(atoi(item.c_str()));
    }

    getline(f_conf,tmp);
    pos = tmp.find(" = ");
    tmp = tmp.substr (pos + 3);
    iterations = atoi(tmp.c_str());

    getline(f_conf,tmp);
    pos = tmp.find(" = ");
    tmp = tmp.substr (pos + 3);
    pMin = atof(tmp.c_str());

    getline(f_conf,tmp);
    pos = tmp.find(" = ");
    tmp = tmp.substr (pos + 3);
    pMax = atof(tmp.c_str());

    getline(f_conf,tmp);
    pos = tmp.find(" = ");
    tmp = tmp.substr (pos + 3);
    pStep = atof(tmp.c_str());

    f_conf.close();
}


// Simulate a simple decoding algorithm. Run simulation and save result.
void Simulation::run (void)
{

    // Initialize the random number generator
    srand (time(0));
    // Initialize constant random number seed
    //srand(8844);


    //loop for different XSizes
    for (int i = 0; i < (int) XSize.size(); i++)
    {

        //data file
        ofstream f_out;
        stringstream ss;

        //add a 0 to the name if XSize is smaller than 0
        if(XSize[i] < 10)
        {
            ss << "out/results_0" << XSize[i] << ".txt";
        }
        else
        {
            ss << "out/results_" << XSize[i] << ".txt";

        }

        //open data file stream
        f_out.open (ss.str().c_str());


        //Print some infromation to screen
        cout   << "RUNNING SIMULATION:   XSize=" << XSize[i] << "   iterations=" << iterations << "   pMin=" << pMin << "   pMax=" << pMax << "   pStep=" << pStep << endl;


        //loop for different p's
        for (double p = pMin; p <= pMax + 0.00001; p += pStep)
        {
            //initialize the counter summing up the successes
            int cntSucc = 0;
            //start parallel programming
#pragma omp parallel for schedule(dynamic) default(shared) reduction(+:cntSucc)
            for (int j = 1; j <= iterations; j++ )
            {

                //There is a lot of debugging in this loop
                //important tasks are marked with ######



                if (debug4)
                {
                    cout<<"iteration: "<<j<<endl;
                }

                if (debug1) cout << "1. ERROR WAS GENERATED:" << endl;

                //###########
                //Build a lattice with errors
                Lattice myLattice (XSize[i], p);
                //###########

                if (debug1) myLattice.printState();
                if (debug3)
                {
                    cout<<"Before Correction:"<<endl;
                    myLattice.printStateQB();
                }

                if (debug1) cout << "2. Z ERRORS CORRECTED USING SITE STABILIZERS:" << endl;

                //############
                //Correct the lattice
                myLattice.correctMesh();
                //############


                if (debug1) myLattice.printState();
                if (debug3)
                {
                    cout<<"After Correction:"<<endl;
                    myLattice.printStateQB();
                }


                if (debug4)
                {
                    if ( !myLattice.is_corrected() )
                    {
                        cout<<"Not Corrected"<<endl;
                        assert(0);
                    }
                }



                //##############
                //Test, if the recovery was successful
                if (myLattice.success (0))
                {
                    cntSucc++;

                    //##############
                    if (debug1) cout << "SIMULATION REPORTS SUCCESS!!!" << endl;

                    if(debug9)
                    {
                        char in;
                        cout<< "Correction Succeded"<<endl;
                        cin>>in;
                    }
                }
                else
                {

                    if(debug7)
                    {
                        char in;
                        cout<< "Correction Failed"<<endl;
                        cin>>in;
                    }
                    if (debug1) cout << "SIMULATION REPORTS FAILURE!!!" << endl;
                }

            }


            //Print the results to the screen
            cout << "Success at error level p=" << 100 * p << "% is " << (double) (100 * cntSucc) / iterations << "%." << endl;

            //Print the resultsto thefile
            f_out << 100 * p << " " <<  (double) (100 * cntSucc) / iterations << endl;


        }

        f_out.close();

    }
}


//end of file simulation.cc
