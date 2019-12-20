//
//  EquationsWormChain.cpp
//  
//
//  Created by Miranda Holmes-Cerfon on July 16, 2019.
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


// ====================================================
//           Stratification variables/functions
// ====================================================

// Scalar function to sample on manifold
double Stratification::feval(Point& x) {
    double kfac = 1.;
    if(x.iL == 1) {
        kfac = 1.;   // put kappa here
    }
    return 1. ;
}


// (Optional) Function to compute statistics on manifolds
void Stratification::computestats(Point& x, VectorXd& stats){
    stats.setZero();
    stats(0) = x.iL;    // = 1 if bonded
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
    
    
    // Set number of statistics to compute
    nstats = 1;
    
    
    int I = Stratification::cIn;
    int E = Stratification::cEq;
    
    nvars = 10; // set to at least 2
    nfcns = 2;

    
    // construct list of edges and manifolds
    nmanif = 2;
    Llist.resize(nmanif,nfcns);
    Llist << I, I,   // interior
              E, I;   // plane
    

    
    // Now set initial condition
    L0 = Llist.row(1);
    p0 = VectorXd::Zero(nvars);
    
    // debug
    cout << "p0= " << p0.transpose() << endl;
    cout << "L0 = " << L0.transpose() << endl;
    // end debug

}



// System of equations
double Stratification::eqs(const int i, const VectorXd& x){
    switch(i) {
        case 0: {           // plane x(0)=0
            return 2*x(0);
            break;
        }
        case 1:            // inequality x(0) < 2
            return 2-x(0);
            break;
        default:
            return NAN;
    }
}


// Jacobian
void Stratification::jac(const int i, const VectorXd &x, VectorXd& row){
    row.setZero();
    switch(i) {
        case 0: {
            row(0) = 2.;
            break;
        }
        case 1:
            row(0) = -1.;
            break;
    }
}



