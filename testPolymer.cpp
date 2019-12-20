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
    int npts = 1e4 ;   // how many points to generate
    int dsave = 1;         // how often to save points
    int ifwritepts = 1;  // set to 1 to write points, labels to file
    
    // Stratification parameters
    int n = 8;        // number of spheres
    int dim = 3;      // dimension of spheres
    double kappa = 5;  // sticky parameter for off-backbone spheres
    VectorXd params(3);
    params << n,dim,kappa;
    int whichcase = 0;

    
    // filename
    string filehead("Tests/PolymerStrat");
    string datafile = filehead+"_data.txt";
    string pointfile = filehead+"_pts.txt";
    string labelfile = filehead+"_labels.txt";
    
    
    
    // Proposal parameters
    int seed = 0;     // = 0 for random seed (based on clock time)
    double sig = 0.3;//0.4;
    double sigbdy = 0.3;   //Point::inf
    double sigtan = 0.2;
    double lamlose = 0.4;
    double lamgain = 4*1./kappa * sigbdy*lamlose;  // smaller for higher kappa



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
    mysampler.print_rejections();
    
    // Print out statistics
    ArrayXd stats = mysampler.stats();
    cout << "Statistics: \n           " << stats.transpose() << endl;
    cout << "Std dev: \n              " << mysampler.statsstd().transpose() << endl;
  
    
    // --------   Write Data to File  --------
    cout << "writing data to file " << endl;
    auto t1 = Clock::now();
    mysampler.write_data(datafile);
    if(ifwritepts == 1) {
        cout << "writing pts, labels to file " << endl;
        mysampler.write_pts(pointfile);
        mysampler.write_labels(labelfile);
    }
    auto t2 = Clock::now();
    cout << "Writing to file took "
    << chrono::duration_cast<chrono::nanoseconds>(t2 - t1).count() / 1e9
    << " seconds" << std::endl;

   
   return 0;
}





