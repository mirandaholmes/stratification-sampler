//
//  test1.cpp
//
//  Test Stratification member functions for dealing with neighbours
//  
//
//  Created by Miranda Holmes-Cerfon on July 16, 2019.
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
    int npts = 1000*100;   // how many points to generate
    int dsave = 1;         // how often to save points
    
    // Stratification parameters
    int whichcase = 0;
    
    // Proposal parameters
    double sig = 1.1;
    double sigbdy = 0.5;   //Point::inf
    double sigtan = 0.9;
    double lamlose = 0.4;
    double lamgain = sigbdy*lamlose;  // 2 manif: want sigbdy*lamlose;     int seed = 0;     // = 0 for random seed (based on clock time)
    int seed = 0;
    
    
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
    ArrayXd stats = mysampler.stats();
    cout << "Statistics: \n           " << stats.transpose() << endl;
    cout << "Std dev: \n              " << mysampler.statsstd().transpose() << endl;
  
    
    // --------   Write Data to File  --------
   /* cout << "writing pts, labels to file " << endl;
    auto t1 = Clock::now();
    mysampler.write_pts("Tests/plane_pts.txt");
    mysampler.write_labels_idx("Tests/plane_labels.txt");
    //mysampler.write_data("Tests/plane_data.txt");
    auto t2 = Clock::now();
    cout << "Writing to file took "
    << chrono::duration_cast<chrono::nanoseconds>(t2 - t1).count() / 1e9
    << " seconds" << std::endl;
    */
   
   return 0;
}





