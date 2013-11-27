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

    getline(f_conf,tmp);
    pos = tmp.find(" = ");
    tmp = tmp.substr (pos + 3);
    q_factor = atof(tmp.c_str());

    f_conf.close();
}


// Simulate a simple decoding algorithm. Run simulation and save result.
void Simulation::run (void)
{

    // Initialize the random number generator.
    srand (time(0));


    for (int i = 0; i < (int) XSize.size(); i++)
    {


        //data files
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
        //Conclusion
        //ofstream conclusion_out;
        //stringstream ss_conclusion;

        //add a 0 to the name if XSize is smaller than 0
        if(XSize[i] < 10)
        {
            ss << "out/results_0" << XSize[i] << ".txt";
        }
        else
        {
            ss << "out/results_" << XSize[i] << ".txt";

        }

        //ss_conclusion<<"out/conclusion.txt";

        f_out.open (ss.str().c_str());
        //f_out.open (ss_conclusion.str().c_str());


        cout   << "RUNNING SIMULATION: XSize=" << XSize[i] << "   iterations=" << iterations << "   pMin=" << pMin << "   pMax=" << pMax << "   pStep=" << pStep << endl;

        cout << "ERROR [%]:\t\t\t\t\t\tSUCCESS [%]"<<endl<<"Level\tEffective\tper_Step\tSyndrome\tSuccess"<< endl;



        for (double p = pMin; p <= pMax + 0.00001; p += pStep)
        {

            if (debug3)
            {
                cout<<"loop p = "<<p<<endl;
            }

            //configuration
            //number of succesful iterations
            int cntSucc = 0;


            //Error probability per time step
            double p_per_step = p;//determine_p_per_step(p,XSize[i]);

            //Syndrome measurement error probability
            double q = p_per_step*q_factor;

            //effective error probability
            double p_effective = determine_p_effective(p_per_step,XSize[i]);


            //Measuerement and Error Repetitions
            int Timesteps = XSize[i];

#pragma omp parallel for schedule(dynamic) default(shared) reduction(+:cntSucc)
            for (int iterate = 1; iterate <= iterations; iterate++ )
            {
                if (debug3)
                {
                    cout<<"loop iterate = "<<iterate<<endl;
                }

                if (debug3)
                {
                    cout<<"Start Function create Lattice"<<endl;
                }

                Lattice myLattice (XSize[i],Timesteps);


                if (debug4)
                {
                    cout<<"Lattice Created"<<endl;
                    cout<<"     ";
                    myLattice.print2DState();
                    cout<<endl;
                }

                if (debug1) cout << "2. Z ERRORS CORRECTED USING SITE STABILIZERS:" << endl;







                //now evolute and correct the lattice
                myLattice.evolute_and_correct(p_per_step,q);

                if (debug4)
                {
                    int k=0;
                    cout<<endl<<"Press any key to continue:";
                    cin>>k;
                    cout<<endl;
                }
                if (debug1) myLattice.printState();

                //now do one last round of perfect correction
                myLattice.perfect_correction();

                if (debug4)
                {
                    int k=0;
                    cout<<endl<<"Press any key to continue:";
                    cin>>k;
                    cout<<endl;
                }
                if (debug1) myLattice.printState();


                if (myLattice.success ())
                {
                    cntSucc++;
                    if (debug1) cout << "SIMULATION REPORTS SUCCESS!!!" << endl;
                }
                else
                {
                    if (debug1) cout << "SIMULATION REPORTS FAILURE!!!" << endl;
                }

                if (debug3)
                {
                    cout<<"Finished Success "<<endl;
                }

                if (debug5 && !(iterations%100))
                {
                    float progress = (float)iterate/(float)iterations;
                    cout<<"\rprogress[%]: "<<setprecision(3)<<100* progress;
                    //if(iterate == iterations-1){cout<<"\r";}

                }



            }//end of loop iterations





            if(debug5)
            {
                cout<<"\r                 \r";
            }
            cout<< 100 * p <<"\t"<<setprecision(5)<<p_effective * 100<<"\t\t" <<p_per_step * 100<<  "\t\t"<<100*q<<"\t\t"<< (double) (100 * cntSucc) / iterations  << endl;
            f_out << 100 * p_per_step<< " " <<  (double) (100 * cntSucc) / iterations << endl;

        }

        f_out.close();

    }
}













double Simulation::determine_p_per_step (double p_G, int T_in)
{

    double p_step = PSTEP;
    double p_max = PMAX;
    double T = (double) T_in;
    //durchlaeufe


    //Gesamtwahrscheinlichkeit nach T-steps
    vector <double> vec_p_eff_T;
    //Effektive Einzelwahrscheinlichkeit per step (DIE WOLLEN WIR WISSEN)
    vector <double> vec_p_per_step;



    for (double p_per_step = 0.0000f; p_per_step <= p_max; p_per_step  += p_step)
    {


        //hier errechne wahrscheinlichkeit nach T schritten
        double p_eff_T = 0.0000f;

        for (int k = 0; k <= T; k++)
        {
            //only sum up odd k's
            if (k%2)
            {
                double k_loop =  k;

                double q = 1- p_per_step;
                double l = T-k_loop;

                double prod1 =ncr(T,k_loop);
                double prod2 = pow(p_per_step,k_loop);
                double prod3 = pow(q,l);

                p_eff_T += prod1*prod3*prod2;
                /*
                cout<<"kloop:"<<k_loop<<endl<<" ncr("<<T<<","<<k<<")="<<prod1<<endl<<"pow("<<p<<","<<k_loop<<")="<<prod2<<endl<<"pow("<<1-p<<","<<T-k_loop<<")="<<prod3<<endl;
                cout<<"p_eff="<<p_eff<<endl;
                */

            }//endofif
        }//endoffor

        vec_p_per_step.push_back(p_per_step);
        vec_p_eff_T.push_back(p_eff_T);

    }

    /*
    	//debug
    	cout<<"p\tp_eff"<<endl;
    	for (int i = 0; i < abs((int)prob.size()); i++) {

    		cout<<prob[i]<<"\t"<<prob_eff[i]<<endl;

    	}
    */

    //now check which p_real fits to which p_0

    double delta_p_wunsch = DELTAPWUNSCH;


    //nun suche ein p das das gewuenschte p gesamt ergibt
    while (true)
    {

        for (int i = 0; i < (int)vec_p_per_step.size(); i++)
        {

            double delta_p = p_G - vec_p_eff_T[i];
            //absolute value
            if (delta_p < 0)
            {
                delta_p *= -1;
            }

            if (delta_p < delta_p_wunsch )
            {
                return vec_p_per_step[i];
            }

        }

        //if there is no p found that fits

        delta_p_wunsch += 0.000005;


        if (delta_p_wunsch > 0.1)
        {
            cout<<"NO PROPER P FOUND!"<<endl;
            break;

        }


    }







    return p_G/ ( (double) T);




}//


double Simulation::determine_p_effective (double p_per_step, int T_in)
{



    //hier errechne wahrscheinlichkeit nach T schritten
    double p_eff_T = 0.0000f;

    for (int k = 0; k <= T_in; k++)
    {
        //only sum up odd k's
        if (k%2)
        {
            double k_loop =  k;

            double q = 1 - p_per_step;
            double l = T_in -k_loop;

            double prod1 = ncr(T_in,k_loop);
            double prod2 = pow(p_per_step,k_loop);
            double prod3 = pow(q,l);

            p_eff_T += prod1*prod3*prod2;

        }//endofif
    }//endoffor


    return p_eff_T;

}//









double Simulation::faculty (double n)
{
    if (n == 0)
    {
        return 1;
    }
    else
    {
        double n_minus_1 = n-1;
        n = n * faculty(n_minus_1);
    }

    return n;

}

double Simulation::ncr (double n, double k)
{

    double n_min_k = n - k;

    double bin =  faculty(n) / (faculty(n_min_k) * faculty(k));

    return bin;


}








//end of file simulation.cc
