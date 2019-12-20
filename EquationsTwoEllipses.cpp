//
//  EquationsTwoEllipses.cpp
//  
//
//  Created by Miranda Holmes-Cerfon on 10/31/17.
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

#include "Stratification.hpp"
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
//   params
//   nstats
// There is an optional input parameter in case there are different cases to consider.
void Stratification::setup(void) {
    
    int I = Stratification::cIn;
    int E = Stratification::cEq;
    int N = Stratification::cNone;
    
    
    nvars = 4;
    nfcns = 3;
    nmanif = 2;
    
    // Set up Llist
    Llist.resize(nmanif,nfcns);
    
    if(whichcase == 0) {
        Llist << E, E, N,   // two ellipses
                  E, I, N;   // one ellipse + interior of ellipse
    }
    if(whichcase == 1) {
        Llist << I, N, N,  // interior of ellipse
                  I, N, E; // interior + line
    }
    
    // number of statistics to compute
    nstats = 2;
    
    // Initial condition
    double a = params[0];
    double b = params[1];
    double c = params[2];
    double d = params[3];
    p0.resize(nvars);
    L0.resize(nfcns);
    if(whichcase == 0) {
        p0 << a,0.,c,0.;  // x=a,y=0,z=c,w=0
        L0 << E, E, N;
    }
    if(whichcase == 1) {
        p0 << 0.,0.,0.,0.;
        L0 << I, N, E;
    }
}


// System of equations
double Stratification::eqs(const int i, const VectorXd& x){
    double a = params[0];
    double b = params[1];
    double c = params[2];
    double d = params[3];
    
    switch(i) {
        case 0: {           // ellipse # 1
            return -(x(0)*x(0)/(a*a) + x(1)*x(1)/(b*b) - 1);  // (x/a)^2 + (y/b)^2 = 1
            break;
        }
        case 1:            // ellipse # 2
            return -(x(2)*x(2)/(c*c) + x(3)*x(3)/(d*d) - 1);  // (z/c)^2 + (w/d)^2 = 1
            break;
        case 2:            // line x=0
            return x(0);
            break;
        default:
            return NAN;
    }
}


// Jacobian
void Stratification::jac(const int i, const VectorXd &x, VectorXd& row){
    double a = params[0];
    double b = params[1];
    double c = params[2];
    double d = params[3];
    
    switch(i) {
        case 0: {
            row(0) = -2*x(0)/(a*a);
            row(1) = -2*x(1)/(b*b);
            row(2) = 0.0;
            row(3) = 0.0;
            break;
        }
        case 1:
            row(0) = 0.0;
            row(1) = 0.0;
            row(2) = -2*x(2)/(c*c);
            row(3) = -2*x(3)/(d*d);
            break;
        case 2:
            row(0) = 1.0;
            row(1) = 0.0;
            row(2) = 0.0;
            row(3) = 0.0;
    }
}

// Scalar function to sample on manifold
double Stratification::feval(Point& x) {
    return 1.0;
}


// (Optional) Function to compute statistics on manifold
void Stratification::computestats(const Point& x, VectorXd& stats){
    stats.setZero();
    if(whichcase == 0)  stats(x.dim-2) += 1.0;
     else if(whichcase == 1) stats(x.dim-3) += 1.0;
}











