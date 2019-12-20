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
    double efac = 1.;

    int n = (int) params(0);  // number of spheres
    int dim = (int) params(1);  // dimension of spheres
    double kappa = params(2);
    double kbend = params(3);  

    for(int i=0; i < (x.neqns-(n-1)); i++) {
        kfac *= kappa;
    }

    // Add in bending energy

    if(whichcase == 1) {
        for(int i=0; i < n-2; i++) {
            double c=0;  // cos theta
            for(int k=0;k<dim; k++) {
                 c += (x.p(dim*(i+1)+k) - x.p(dim*(i)+k)) * (x.p(dim*(i+2)+k) - x.p(dim*(i+1)+k)) ;
            }
            efac *= exp(- kbend * (1-c));
        }
    }

    return kfac * efac * Cluster::getH(x);   
}


// (Optional) Function to compute statistics on manifolds
void Stratification::computestats(Point& x, VectorXd& stats){
    stats.setZero();
    int dim = params(1);  // dimension of particles
    int n = params(0);
    stats(0) = x.neqns - (n-1);   // number of wall bonds
    stats(1) = Cluster::dist(x.p,0,n-1,dim);  // distance 1-last
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
    nstats = 2;
    
    // other system-dependent parameters
    int n = (int) params[0];  // number of spheres
    int dim = (int) params[1];  // dimension of spheres
    nvars = n*dim;   // number of variables to describe system
    nedges = n-1;  // number of edges
    nfcns = nedges + n;   // total number of functions
    
    //debug
    cout << "n = " << n << endl;
    cout << "dim = " << dim << endl;
    cout << "k = " << params[2] << endl;
    cout << "kbend = " << params[3] << endl;
    //end debug
    
    // abbreviations, for ease of reference
    int I = Stratification::cIn;
    int E = Stratification::cEq;
    
    
    // construct list of edges and fixed constraints
    FixFcns = VectorXi::Constant(nfcns, cVary);
    L0 = VectorXi::Constant(nfcns,I);
    edges.resize(nedges,2);
    for(int k=0; k < n-1; k++) {
        FixFcns(k) = cFix;
        L0(k) = E;
        edges(k,0) = k;
        edges(k,1) = k+1;
    }
    FixFcns(n-1) = cFix;        // tie 1st particle to wall
    L0(n-1) = E;
    

/*
    // Initial condition: zig-zaggy line
    p0 = VectorXd::Zero(nvars);
    for(int i=0; i<n/2; i++) {
        p0(2*i*dim) = i;
        p0(2*i*dim + dim-1) = i;
        p0((2*i+1)*dim) = i;
        p0((2*i+1)*dim + dim-1) = i + 1;
    }
*/

    // Initial condition: zig-zaggy line in wall
    L0 = VectorXi::Constant(nfcns,E);
    p0 = VectorXd::Zero(nvars);
    for(int i=0; i<n/2; i++) {
        p0(2*i*dim) = i;
        p0(2*i*dim + 1) = i;
        p0((2*i+1)*dim) = i;
        p0((2*i+1)*dim + 1) = i + 1;
    }
    
    /*
    // debug
    cout << "p0= " << p0.transpose() << endl;
    cout << "L0 = " << L0.transpose() << endl;
    // end debug
    */
}



// System of equations
double Stratification::eqs(const int i, const VectorXd& x){

    int dim = params(1);

    // function is an edge constraint
    if(i < nedges) {
        int ir = edges(i,0);
        int ic = edges(i,1);
        double val = -1.0;
        for(int j=0; j<dim; j++){     // loop through dimensions
            val += (x(ir*dim+j) - x(ic*dim+j))*(x(ir*dim+j) - x(ic*dim+j));
        }
        return val;
    }

    // function is a wall constraint
    int ip = i-nedges;   // index of particle being constrained
    return x(ip*dim + dim-1);  // last dimension = 0. 

}


// Jacobian
void Stratification::jac(const int i, const VectorXd &x, VectorXd& row){
    
    int dim = params(1);
    row.setZero();

    // function is an edge constraint
    if(i < nedges) {
        int ir = edges(i,0);
        int ic = edges(i,1);
        for(int j=0; j<dim; j++){     // loop through dimensions
          row(ir*dim+j) =  2*(x(ir*dim+j)-x(ic*dim+j));
          row(ic*dim+j) =  -2*(x(ir*dim+j)-x(ic*dim+j));
      }
    }
    else {
        int ip = i-nedges;   // index of particle being constrained
        row(ip*dim + dim-1) = 1;
    }
}



