//
//  EquationsPolymer.cpp
//  
//
//  Created by Miranda Holmes-Cerfon on September 3, 2019.
//
//
//  Define system-dependent quantities needed for Stratification object.
//  Must provide:
//    -- setup: sets variables and manifold list
//    -- eqs:   list of functions representing equations / inequalities
//    -- jac:   jacobian of above functions
//    -- feval: scalar function to sample on manifold
//    -- computestats: statistics to compute along the way
//
//  If the user wants to add more member variables or functions, these must be defined
//  in Stratification.hpp.
//  Or, another option would be to include a separate namespace in this file.
//
//
//

#include "Stratification.hpp"
#include "SampleStrat.hpp"
#include "Point.hpp"
#include "Move.hpp"
#include "Cluster.hpp"


// ====================================================
//           Stratification variables/functions
// ====================================================

// Scalar function to sample on manifold
double Stratification::feval(Point& x) {

    return 1.0;   
}


// (Optional) Function to compute statistics on manifolds
void Stratification::computestats(Point& x, VectorXd& stats){
    stats.setZero();
    stats(0) = x.neqns;   // number of equations
    if(x.neqns == 1) {
    	stats(1) = 1.;   // indicator function for surface of ellipsoid
    }
    if(whichcase == 1 && x.neqns == 0) {
    	stats(2) = 1.;   // indicator function for interior of ellipsoid
    }
    if(whichcase == 0 && x.neqns == params(0)) {
    	stats(2) = 1.;   // indicator function for 0-dimensional points
    }
}



// Initialize variables for Stratification
// Must set:
//   nvars
//   nfcns
//   nmanif
//   Llist (size, and values)
//   p0, L0 = initial condition
//   nstats
// Use whichcase or params or edges to set optional parameters
void Stratification::setup(void) {
    
    // tells sampler to find manifolds on the fly; no list of manifolds is provided
    nmanif = 0;
    
    // Number of statistics to compute
    nstats = 3;
    
    // other system-dependent parameters
    int n = (int) params[0];  // dimension of problem / number of variables
    nvars = n;   // number of variables to describe system
    if(whichcase == 0) {
    	nfcns = n;   // total number of functions
    }
    else {
    	nfcns = 1;
    }
     
    //debug
    cout << "n = " << n << endl;
    //end debug
    
    // abbreviations, for ease of reference
    int I = Stratification::cIn;
    int E = Stratification::cEq;
    
    
    // construct list of edges and fixed constraints
    FixFcns = VectorXi::Constant(nfcns, cVary);
    L0 = VectorXi::Constant(nfcns,E);
    if(whichcase == 0) {
        FixFcns(0) = cFix;
    }
    
    // Initial condition
    p0 = VectorXd::Zero(nvars);
    p0(0) = params[1];

    /*
    // debug
    cout << "p0= " << p0.transpose() << endl;
    cout << "L0 = " << L0.transpose() << endl;
    // end debug
    */
}



// System of equations
double Stratification::eqs(const int i, const VectorXd& x){

    int n = params(0);

    // ellipse constraint
    if(i == 0) {
        double val = 1.0;
        for(int j=0; j<n; j++){     // loop through dimensions
            val -= x(j)*x(j) / (params(j+1)*params(j+1));
        }
        return val;
    }

    // plane constraint; x(i) = 0 (i=1..n; skipping 0)
    else return x(i);

}


// Jacobian
void Stratification::jac(const int i, const VectorXd &x, VectorXd& row){
    
    int n = params(0);
    row.setZero();

    // ellipse constraint
    if(i==0) {
    	for(int j=0; j<n; j++) {
        	row(j) = -2.0*x(j) / (params(j+1)*params(j+1));
    	}
    }
    else {
        row(i) = 1;  // index of particle being constrained
    }
}



