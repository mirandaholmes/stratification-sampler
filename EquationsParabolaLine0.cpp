//
//  EquationsParabolaLine.cpp
//  
//
//  Created by Miranda Holmes-Cerfon on 10/31/17.
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


// Initialize variables for Stratification
// Must set:
//   nvars
//   nfcns
//   nmanif
//   Llist (size, and values)
//   p0, L0 = initial condition
//   nstats
// There is an optional input parameter in case there are different cases to consider.
void Stratification::setup(void) {
    
    nvars = 3;
    nfcns = 3;
    
    int I = Stratification::cIn;
    int E = Stratification::cEq;

    nmanif = 0;
    FixFcns.resize(nfcns);
    FixFcns << cVary,cFix,cVary;
    
    L0.resize(nfcns);
    p0.resize(nvars);
    
    L0 << I,E,I;
    p0 << 0,0.5,0;

    // Set number of statistics to compute
    nstats = 3;
}


// System of equations
double Stratification::eqs(const int i, const VectorXd& x){
    switch(i) {
        case 0:
            return x(1)-x(0)*x(0);  // y - x^2
            break;
        case 1:
            return x(2);            // z = 0
        case 2:
            return -(x(1) - 2.0);   // const-y
            break;
        default:
            return NAN;
    }
}


// Jacobian
void Stratification::jac(const int i, const VectorXd &x, VectorXd& row){
    switch(i) {
        case 0:
            row(0) = -2*x(0);
            row(1) = 1.0;
            row(2) = 0;
            break;
        case 1:
            row(0) = 0;
            row(1) = 0;
            row(2) = 1.;
            break;
        case 2:
            row(0) = 0.0;
            row(1) = -1.0;
            row(2) = 0;
            break;
    }
}



// Scalar function to sample on manifold
double Stratification::feval(Point& x) {
    return 1.0;
}



// ====================================================
//           Statistics to sample
// ====================================================

// (Optional) Function to compute statistics on manifold
// Must return vector of size nstats
void Stratification::computestats(Point& x, VectorXd& stats){
    
    stats.setZero();
    stats(x.dim) += 1.0; // nstats = 3;
    
}










