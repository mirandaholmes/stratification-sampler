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
    double kfac = 1.;
    for(int i=0; i<x.neqns; i++) {
        kfac *= params(2);
    }
    return kfac * Cluster::getH(x);   // same sticky parameter for every bond
}


// (Optional) Function to compute statistics on manifolds
void Stratification::computestats(Point& x, VectorXd& stats){
    stats.setZero();
    int dim = x.pstrat->params(1);  // dimension of particles
    int n = x.pstrat->params(0);
    stats(0) = x.neqns;   // number of bonds
    //stats(1) = Cluster::dist(x.p,0,n-1,dim);  // distance 1-last
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
    nstats = 1;
    
    // other system-dependent parameters
    int n = (int) params[0];  // number of spheres
    int dim = (int) params[1];  // dimension of spheres
    nvars = n*dim;   // number of variables to describe system
    nedges = n*(n-1)/2;  // number of edges
    nfcns = nedges;      // total number of functions; no other constraints for now
    
    //debug
    cout << "n = " << n << endl;
    cout << "dim = " << dim << endl;
    cout << "k = " << params[2] << endl;
    //end debug
    
    // abbreviations, for ease of reference
    int I = Stratification::cIn;
    int E = Stratification::cEq;
    
    
    // construct list of edges and fixed edges
    FixFcns = VectorXi::Constant(nfcns, cVary);
    L0 = VectorXi::Constant(nfcns,I);
    edges.resize(nfcns,2);
    int k = 0;
    for(int ir=0; ir < n; ir++) {
        for(int ic=ir+1; ic<n; ic++) {
            edges(k,0) = ir;
            edges(k,1) = ic;
            if(ic == ir+1) {      // worm constraints
                FixFcns(k) = cFix;
                L0(k) = E;
            }
            k++;
        }
    }
    
    /*// debug
    cout << "Llist = \n" << Llist << endl;
    cout << "edges = \n" << edges.transpose() << endl;
    // end debug */

    
    // Now set initial condition
    p0 = VectorXd::Zero(nvars);
    // straight line, bent a little at the end
    
    for(int i=0; i<n; i++) {     
        p0(i*dim) = i-(n-2.)/2.;
    }
    p0((n-1)*dim) = n-2.-(n-2.)/2.; // now bend it a little
    p0((n-1)*dim+1) = 1;
    
/*
    // zig-zaggy line
    for(int i=0; i<n/2; i++) {
    	p0(2*i*dim) = i-(n/4.);
    	p0(2*i*dim + 1) = i-(n/4.);
    	p0((2*i+1)*dim) = i-(n/4.) + 1;
    	p0((2*i+1)*dim + 1) = i-(n/4.);
    }
*/
    
    // debug
    //cout << "p0= " << p0.transpose() << endl;
    //cout << "L0 = " << L0.transpose() << endl;
    // end debug

}



// System of equations
double Stratification::eqs(const int i, const VectorXd& x){
    int ir = edges(i,0);
    int ic = edges(i,1);
    int dim = params(1);
    
    double val = -1.0;
    for(int j=0; j<dim; j++){     // loop through dimensions
        val += (x(ir*dim+j) - x(ic*dim+j))*(x(ir*dim+j) - x(ic*dim+j));
    }
    return val;
}


// Jacobian
void Stratification::jac(const int i, const VectorXd &x, VectorXd& row){
    int ir = edges(i,0);
    int ic = edges(i,1);
    int dim = params(1);
    row.setZero();
    for(int j=0; j<dim; j++){     // loop through dimensions
        row(ir*dim+j) =  2*(x(ir*dim+j)-x(ic*dim+j));
        row(ic*dim+j) =  -2*(x(ir*dim+j)-x(ic*dim+j));
    }
}



