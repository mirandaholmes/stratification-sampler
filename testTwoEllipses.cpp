//
//  test1.cpp
//
//  Test Stratification member functions for dealing with neighbours
//  
//
//  Created by Miranda Holmes-Cerfon on 7/10/17.
//
//

#include <stdio.h>
#include "Stratification.hpp"
#include "Point.hpp"
#include "Move.hpp"
#include "Proposals.hpp"
#include "SampleStrat.hpp"
#include <Eigen/Dense>
#include <chrono>   // for high-precision timing
#include <cmath>


typedef std::chrono::high_resolution_clock Clock;

using namespace Eigen;
using namespace std;



// Main loop
int main(int argc,char *argv[])
{
    // Sampling parameters
    int npts = 1000*100;
    int dsave = 5;
    
    // Stratification parameters
    int whichcase = 0;   // which stratification to choose
    int paramcase = 2;  // which parameters to choose
    
    // Proposal Parameters
    double sig = 1.6; //0.7;
    double sigbdy = 0.9; //0.4;   //Point::inf
    double sigtan = 1.4;
    double lamlose = 0.7;
    double lamgain = sigbdy*lamlose;  // 2 manif: want sigbdy*lamlose; 
    int seed = 0;     // = 0 for random seed (based on clock time)
    
    
    // parameters for ellipses
    double a; double b; double c; double d;
    switch(paramcase) {
        case 1:        // works
            a = 1; b=1; c=1; d=1;
            break;
        case 2:      // seems to work (ratios vary quite a lot, but in both directions)
            a = 1; b=1.5; c=0.8; d=2;
            break;
        case 3:      //  works
            a = 1; b = 1.5; c = 0.5; d = 0.8;
            break;
        case 4:      // works
            a = 1; b = 2; c = 3; d = 5;
            break;
        case 5:       // works
            a = 1; b = 2; c = 0.5; d = 5;
            break;
    }

    VectorXd params(4);
    params << a,b,c,d;
    
    // Create Stratification object
    Stratification mystrat(whichcase,params);

    Proposals myprops(sig,sigbdy,sigtan,lamgain,lamlose,seed);
    
    cout << "seed = " << myprops.seed() << endl;
    cout << "lamlose = " << lamlose << endl;
    cout << "lamgain = " << lamgain << ", should = " << sigbdy*lamlose << endl;
    cout << "lamsame = " << 1-lamlose-lamgain << endl;
    cout << "sig = " << sig << ", sigbdy = " << sigbdy << ", sigtan = " << sigtan << endl;
    cout << "Parameters = " << params.transpose() << endl;


    
     // --------   Now test SampleStrat   --------
    // Initialize sampler
    SampleStrat mysampler(&mystrat,&myprops,mystrat.p0,mystrat.L0);
    //mysampler.setIfDebug(SampleStrat::cYes);

    
    // Sample!
    mysampler.sample(npts, dsave);

    
    // Print out rejection statistics
    //mysampler.print_rejections_manifold();
    mysampler.print_rejections();
    
    // Print out statistics
    VectorXd stats = mysampler.stats();
    cout << "Statistics: \n" << stats.transpose() << endl;
    if(whichcase == 0) {
        //cout << "Ratios: \n" << stats.transpose()/stats(mysampler.nstats()-1) << endl;
        cout << "Ratios: " << stats.transpose()/stats(1)  << endl;
        double theory = 0;
        switch(paramcase) {
            case 1:
                theory = 2; break;
            case 2:
                theory = 1.8313; break;
            case 3:
                theory = 3.2934; break;
            case 4:
                theory = 0.5417; break;
            case 5:
                theory = 2.5872; break;
        }
        cout << "Theory: " << theory << "\n" << endl;
    }
    if(whichcase == 1) {
        cout << "Ratios: " << stats.transpose()/stats(1)  << endl;
        cout << "Theory: " << 2.0/Proposals::pi*a << "\n" << endl;
    }
    
    
    // Save points and flags to file
    cout << "writing pts, labels to file " << endl;
    auto t1 = Clock::now();
    mysampler.write_pts("Tests/twoellipses_pts.txt");
    mysampler.write_write_labels("Tests/twoellipses_labels.txt");
    auto t2 = Clock::now();
    cout << "Writing to file took "
    << chrono::duration_cast<chrono::nanoseconds>(t2 - t1).count() / 1e9
    << " seconds" << std::endl;

    
    return 0;
}





