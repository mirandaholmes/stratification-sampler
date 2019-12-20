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
    
    nvars = 2;
    nfcns = 2;
    
    int I = Stratification::cIn;
    int E = Stratification::cEq;

    
    // whole list of manifolds
    if(whichcase == 0) {
        nmanif = 4;
        // Set up Llist
        Llist.resize(nmanif,nfcns);
        Llist << I, I,   // interior
                  E, I,   // parabola
                  I, E,   // line
                  E, E;   // corner
    }
    
    // parabola + interior
    if(whichcase == 1) {
        nmanif = 2;
        Llist.resize(nmanif,nfcns);
        Llist << I, I,   // interior
                  E, I;   // parabola
    }
    
    // line + interior
    if(whichcase == 2) {
        nmanif = 2;
        Llist.resize(nmanif,nfcns);
        Llist << I, I,   // interior
                  I, E;   // line
    }
    
    // parabola + corners
    if(whichcase == 3) {
        nmanif = 2;
        Llist.resize(nmanif,nfcns);
        Llist << E, I,   // parabola
                  E, E;   // corner
    }
    
    // line + corners
    if(whichcase == 4) {
        nmanif = 2;
        Llist.resize(nmanif,nfcns);
        Llist << I, E,   // line
                  E, E;   // corner
    }
    
    // Now set initial conditions
    // possible initial conditions
    VectorXd pint(nvars), ppar(nvars), pcor(nvars), plin(nvars);
    VectorXi Lint(nfcns), Lpar(nfcns), Lcor(nfcns), Llin(nfcns);
    double sqrt2 = 1.414213562373095;
    pint << 0,0.5;
    ppar << 0,0;
    plin << 1.2,2;
    pcor << sqrt2, 2;
    Lint << I, I;
    Lpar << E, I;
    Llin << I, E;
    Lcor << E, E;
    
    // initial point (must be on stratification)
    p0 = pcor;
    L0 = Lcor;
    if(whichcase >= 0 && whichcase <= 2){
        p0 = pint;
        L0 = Lint;
    }
    
    // Set number of statistics to compute
    if(whichcase >= 1 || whichcase <= 4) nstats = 2;
    if(whichcase == 0) nstats = 3;
}


// System of equations
double Stratification::eqs(const int i, const VectorXd& x){
    switch(i) {
        case 0:
            return x(1)-x(0)*x(0);  // y - x^2
            break;
        case 1:
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
            break;
        case 1:
            row(0) = 0.0;
            row(1) = -1.0;
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
    
    if(whichcase == 3 || whichcase == 4) {
        stats(x.dim) += 1.0;
        // nstats = 2;
    }
    if(whichcase == 1 || whichcase == 2) {
        stats(x.dim-1) += 1.0;
        // nstats = 2;
    }
    if(whichcase == 0) {
        stats(x.dim) += 1.0;
        // nstats = 3;
    }
}










