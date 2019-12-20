//
//  
//
//  Polymer adsorbing weakly to a wall, in d dimensions. 
//  The wall is defined by {(d-coord) = 0}. 
//  Last particle is tethered to the wall. 
// 
//  Option to add bending stiffness, to make it semi-flexible. 
//  No excluded volume interactions in this code; could add them. 
//  
//
//  Created by Miranda Holmes-Cerfon on November 1, 2019.
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
    int dsave = 5;         // how often to save points
    
    int ifwritepts = 0;  // set to 1 to write points, labels to file
    
    // Stratification parameters
    int n = 10;        // dimension of ambient space
    VectorXd params(n+1);
    params << n,2,2,2,2,3,3,3,1,1,1;  // should have n parameters following the n, 
    				    // equal to the n axes of the ellipse

    int whichcase = 1;  // 0 = surface area to 0-dim, 1 = surface area to vol

    
    // filename
    string filehead("Tests/ellipsoid");
    string datafile = filehead+"_data.txt";
    string pointfile = filehead+"_pts.txt";
    string labelfile = filehead+"_labels.txt";
    
    
    
    // Proposal parameters
    int seed = 0;     // = 0 for random seed (based on clock time)
    double sig = 0.6; 
    double sigtan = 0.3; 
    double sigbdy = 0.4;   
    double lamlose = 0.4;
    double lamgain =  sigbdy*lamlose;  



    /*  ------  Begin code  ------ */

    // Create Stratification object
    Stratification mystrat(whichcase,params);

    
    // Create a proposal object
    Proposals myprops(sig,sigbdy,sigtan,lamgain,lamlose,seed);
    
    cout << "seed = " << myprops.seed() << endl;
    cout << "whichcase = " << whichcase << endl;
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
  

    // compute factors needed to estimate surface area
    double pfac = 1.;   // product a1*a2*...*an
    for (int i=0;i<n; i++) {
        pfac *= params(i+1);
    }
    double vfac = 0;  // factor needed to get volume of n-dimensional sphere
    if(n==3) {
        vfac = 4.1887902;
    }
    if(n==10) {
        vfac = 2.550164;
    }

    //cout << "tfac = " << tfac << ", pfac = " << pfac << endl;


    // Print out surface area estimate
    cout << "Surface Area Estimate: ";
    if(whichcase == 0) {
        cout << stats(1)/stats(2) * 2  * exp(n*0.94)/exp(1*0.94)<< endl;
    }
    if(whichcase == 1) {
        cout << stats(1)/stats(2) * vfac * pfac << endl;
    }
     

    if(n==3) {
        double p = 1.6075;
        double ap = pow(params(1),p);
        double bp = pow(params(2),p);
        double cp = pow(params(3),p);
        cout << "Actual Surface Area: " << 12.5663706 * pow((ap*bp+ap*cp+bp*cp)/3,1/p) << endl;
    }


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





