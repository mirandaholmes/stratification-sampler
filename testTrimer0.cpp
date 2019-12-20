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
    int npts = 1e5;   // how many points to generate
    int dsave = 10;         // how often to save points
    
    // Stratification parameters
    int whichcase = 0;
    int n = 3;        // number of spheres
    int dim = 2;      // dimension of spheres
    double kappa = 1;  // sticky parameter for making a loop
    VectorXd params(3);
    params << n,dim,kappa;

    // filename
    string filehead("Tests/trimer0Anthony");
    string datafile = filehead+"_data.txt";
    string pointfile = filehead+"_pts.txt";
    string labelfile = filehead+"_labels.txt";
    
    // Proposal parameters
    double sig = 0.5;
    double sigbdy = 0.4;   //Point::inf
    double sigtan = 0.3;
    double lamlose = 0.7;
    double lamgain = sigbdy*lamlose;  // 2 manif: want sigbdy*lamlose; 
    int seed = 0;     // = 0 for random seed (based on clock time)
    
    
    
    /*  ------  Begin code  ------ */

    // Create Stratification object
    Stratification mystrat(whichcase,params);
    
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
    mysampler.print_rejections_manifold();
    mysampler.print_rejections();
    
    // Print out statistics
    ArrayXd stats = mysampler.stats();
    cout << "Statistics: \n           " << stats.transpose() << endl;
    cout << "Std dev: \n              " << mysampler.statsstd().transpose() << endl;
    //cout << "f avg / manifold: " << stats(2)/stats(0)<< ", " << stats(3)/stats(1) << endl;
    
    
    // --------   Write Data to File  --------
    cout << "writing pts, labels to file " << endl;
    auto t1 = Clock::now();
    //mysampler.write_pts(pointfile);
    mysampler.write_labels(labelfile);
    mysampler.write_data(datafile);
    auto t2 = Clock::now();
    cout << "Writing to file took "
    << chrono::duration_cast<chrono::nanoseconds>(t2 - t1).count() / 1e9
    << " seconds" << std::endl;

   
   return 0;
}





