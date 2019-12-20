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
#include "Cluster.hpp"


// ====================================================
//           Stratification variables/functions
// ====================================================

// Scalar function to sample on manifold
double Stratification::feval(Point& x) {
    double kfac = 1.;
    if(x.iL == 1) {
        kfac = params(2);   // put kappa here
    }
    return kfac * Cluster::getH(x);
    return Cluster::getH(x);
}


// (Optional) Function to compute statistics on manifolds
void Stratification::computestats(Point& x, VectorXd& stats){
    stats.setZero();
    int dim = x.pstrat->params(1);  // dimension of particles
    int n = x.pstrat->params(0);
    stats(0) = x.iL;    // = 1 if bonded
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
    
    
    // Set number of statistics to compute
    nstats = 2;
    
    
    int I = Stratification::cIn;
    int E = Stratification::cEq;
    
    int n = (int) params[0];
    int dim = (int) params[1];
    nvars = n*dim;
    nfcns = n*(n-1)/2;
    nedges = nfcns;      // no other constraints for now
    
    //debug
    cout << "n = " << n << endl;
    cout << "dim = " << dim << endl;
    cout << "k = " << params[2] << endl;
    //end debug
    
    // construct list of edges and manifolds
    nmanif = 2;
    Llist = MatrixXi::Constant(nmanif,nfcns,I);  // all inequalities to start
    edges.resize(nfcns,2);
    int k = 0;
    for(int ir=0; ir < n; ir++) {
        for(int ic=ir+1; ic<n; ic++) {
            edges(k,0) = ir;
            edges(k,1) = ic;
            if(ic == ir+1) {      // worm constraints
                Llist(0,k) = E;
                Llist(1,k) = E;
            }
            k++;
        }
    }
    Llist(1,n-2) = E;    // loop constraint
    
    /*// debug
    cout << "Llist = \n" << Llist << endl;
    cout << "edges = \n" << edges.transpose() << endl;
    // end debug */

    
    // Now set initial condition
    L0 = Llist.row(0);
    p0 = VectorXd::Zero(nvars);
    for(int i=0; i<n; i++) {     // straight line (doesn't work; bend a little)
        p0(i*dim) = i-(n-2.)/2.;
    }
    p0((n-1)*dim) = n-2.-(n-2.)/2.;
    p0((n-1)*dim+1) = 1;
    
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





//namespace Cluster {

    /* UPDATE -- arbitrary dimension
    // Get center of mass of cluster
    void getCenterOfMass(Point& x, double &xc, double &yc, double &zc){
        int n = getN(x);
        xc=0; yc=0; zc=0;
        for (int i=0; i<n; i++) {
            xc += x(3*i+0);
            yc += x(3*i+1);
            zc += x(3*i+2);
        }
        xc /= n;
        yc /= n;
        zc /= n;
    }
     */
    
    /* UPDATE
    // Return the square root of the moment of inertia tensor (3d cluster)
    double getI3d(const VectorXd &x0) {
        int n = getN(x0);     // number of spheres
        double xc, yc, zc;   // coordinates of center of mass
        Matrix3d matI;       //  moment of inertia tensor
        getCenterOfMass(x0,xc,yc,zc);  // get center of mass
        double x,y,z;        // hold values temporarily
        
        // fill up moment of inertia tensor
        matI.fill(0);        // initialize to 0
        for (int i=0; i<n; i++) {
            x = x0(3*i + 0) - xc;
            y = x0(3*i + 1) - yc;
            z = x0(3*i + 2) - zc;
            matI(0,0) += y*y + z*z;
            matI(1,1) += x*x + z*z;
            matI(2,2) += x*x + y*y;
            matI(0,1) += -x*y;
            matI(0,2) += -x*z;
            matI(1,2) += -y*z;
        }
        matI(1,0) = matI(0,1);
        matI(2,0) = matI(0,2);
        matI(2,1) = matI(1,2);
        
        return sqrt(matI.determinant());
    }
     */
    
    /* UPDATE
    // Return vector of moments of inertias, given a matrix of cluster positions
    VectorXd getI(const MatrixXd &data) {
        int n = getN(data.row(0));
        int numpts = data.rows();
        VectorXd vecI(numpts);
        for (int i=0; i<numpts; i++) {
            vecI(i) = getI3d(data.row(i));
        }
        return vecI;
    }
     */
    
    
    /* UPDATE
     
    // Construct the vibrational factor for a list of clusters
    // T holds an equation class that can evaluate df
    // ** this function depends on specific form of the jacobian, for clusters **
    template <class T>
    VectorXd getH(const MatrixXd &data, const T& eqobj) {
        
        const double m_nulltol = 1e-8;   // tolerance level (eventually make class parameter)
        int numpts = data.rows();
        
        VectorXd  h(numpts);    // holds vibrational factors for each cluster
        MatrixXd R(eqobj.neqns(), eqobj.nvars());  // holds rigidity matrix
        VectorXd s(eqobj.neqns());     // holds singular values
        
        for(int k=0; k<numpts; k++) {
            eqobj.df(data.row(k),R);                      // set R to be jacobian
            R /= 2;                          // now it's the rigidity matrix
            s = R.jacobiSvd().singularValues();
            h(k) = 1;
            for (int i=0; i<s.size(); i++){
                if(abs(s(i)) > m_nulltol) {
                    h(k) /= s(i);
                }
            }
        }
        return h;
    }
     */
//}










