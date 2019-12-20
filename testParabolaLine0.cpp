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
    
     /*  ------  Set up parameters  ------ */
    // Sampling parameters
    int npts = 1e7;   // how many points to generate
    int dsave = 10;         // how often to save points
    
    // Stratification parameters
    int whichcase = 0;
    
    // Proposal parameters
    double sig = 0.9;
    double sigbdy = 0.3;   //Point::inf
    double sigtan = 0.6;
    double lamlose = 0.7;
    double lamgain = sigbdy*lamlose;  // 2 manif: want sigbdy*lamlose; 
    int seed = 0;     // = 0 for random seed (based on clock time)
    
    // filename
    string filehead("Tests/parline0");
    string datafile = filehead+"_data.txt";
    string pointfile = filehead+"_pts.txt";
    string labelfile = filehead+"_labels.txt";
    
    
    /*  ------  Begin code  ------ */

    // Create Stratification object
    Stratification mystrat(whichcase);
    
    // Create a proposal object
    Proposals myprops(sig,sigbdy,sigtan,lamgain,lamlose,seed);
    
    cout << "seed = " << myprops.seed() << endl;
    cout << "lamlose = " << lamlose << endl;
    cout << "lamgain = " << lamgain << ", should = " << sigbdy*lamlose << endl;
    cout << "lamsame = " << 1-lamlose-lamgain << endl;
    cout << "sig = " << sig << ", sigbdy = " << sigbdy << ", sigtan = " << sigtan << endl;
    
    
    // --------   Run SampleStrat   --------
    // Initialize sampler
    SampleStrat mysampler(&mystrat,&myprops,mystrat.p0,mystrat.L0);
    mysampler.setIfDebug(SampleStrat::cNo);  // Set debug flag
    
    // Sample!
    mysampler.sample(npts, dsave);
    
    
    // --------   Print Statistics and Diagnostics  --------
    
    // Print out rejection statistics
    //mysampler.print_rejections_manifold();
    mysampler.print_rejections();
    
    // Print out statistics
    VectorXd stats = mysampler.stats();
    cout << "Statistics: \n           " << stats.transpose() << endl;
    if(whichcase == 0) cout << "Should be  0.1457, 0.5795, 0.2748" << endl;
    //cout << "Ratios: \n" << stats.transpose()/stats(mysampler.nstats()-1) << endl;
    cout << "Ratios: \n" << stats.transpose()/stats(0)  << endl;
    
    if(whichcase == 0) cout << "Should be  1, 3.9762, 1.8856" << endl;
    if(whichcase == 1) cout << "Should be " << 0.735992 << endl;
    if(whichcase == 2) cout << "Should be " << 1.3333333 << endl;
    if(whichcase == 3) cout << "Should be " << 2.5620071 << endl;
    if(whichcase == 4) cout << "Should be " << sqrt(2) << endl;
    
    
    
    // --------   Write Data to File  --------
    
    // Save points and flags to file
    cout << "writing pts, labels to file " << endl;
    auto t1 = Clock::now();
    mysampler.write_pts(pointfile);
    mysampler.write_labels(labelfile);
    auto t2 = Clock::now();
    cout << "Writing to file took "
    << chrono::duration_cast<chrono::nanoseconds>(t2 - t1).count() / 1e9
    << " seconds" << std::endl;

    
    
    return 0;
}





