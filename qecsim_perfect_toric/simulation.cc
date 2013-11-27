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

    // Initialize the random number generator.
    srand (time(0));


    for (int i = 0; i < (int) XSize.size(); i++)
    {

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
        f_out.open (ss.str().c_str());

        cout   << "RUNNING SIMULATION:   XSize=" << XSize[i] << "   iterations=" << iterations << "   pMin=" << pMin << "   pMax=" << pMax << "   pStep=" << pStep << endl;

        for (double p = pMin; p <= pMax + 0.00001; p += pStep)
        {
            int cntSucc = 0;
#pragma omp parallel for schedule(dynamic) default(shared) reduction(+:cntSucc)
            for (int j = 1; j <= iterations; j++ )
            {



                if (debug1) cout << "1. ERROR WAS GENERATED:" << endl;
                Lattice myLattice (XSize[i],p);
                if (debug1) myLattice.printState();

                if (debug1) cout << "2. Z ERRORS CORRECTED USING SITE STABILIZERS:" << endl;


                myLattice.correctMesh();
                if (debug1) myLattice.printState();



                if (myLattice.success (0))
                {
                    cntSucc++;
                    if (debug1) cout << "SIMULATION REPORTS SUCCESS!!!" << endl;
                }
                else
                {
                    if (debug1) cout << "SIMULATION REPORTS FAILURE!!!" << endl;
                }


            }


            cout << "Success at error level p=" << 100 * p << "% is " << (double) (100 * cntSucc) / iterations << "%." << endl;
            f_out << 100 * p << " " <<  (double) (100 * cntSucc) / iterations << endl;

        }

        f_out.close();

    }
}
