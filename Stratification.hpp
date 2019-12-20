//
//  Stratification.hpp
//  
//
//  Created by Miranda Holmes-Cerfon on 3/2/17.
//
//

#ifndef Stratification_hpp
#define Stratification_hpp

#include <stdio.h>
#include <iostream>
#include <Eigen/Dense>


using namespace std;
using namespace Eigen;


// Tell the compiler there is a class called "Point"
class Point;



class Stratification {
private:

    // Functions to be provided in a separate file
    void setup(void);   // initialize the variables above
    double eqs(const int, const VectorXd&);  // system of equations/inequalities
    void jac(const int, const VectorXd&, VectorXd&);  // Jacobian of above
    

    // hold information about neighbours, when we have the whole list (Llist)
    VectorXi mnbrN_gain;  // number of neighbours for each manifold, type gain dimension
    VectorXi mnbrN_lose;  // number of neighbours for each manifold, type lose dimension
    MatrixXi mnbrIL_gain;  // list of neighbours, in order. Row i is indices of nbrs of i.
    MatrixXi mnbrIL_lose;
    MatrixXi mnbrIchange_gain;  // each row holds Ichange for each neighbour
    MatrixXi mnbrIchange_lose;
    // note: storing neighbours and indices as matrices is not that efficient;
    // it would be more efficient to store them sparsely, as pairs or as pointers
    // to lists of different sizes. However, this would require more complex
    // memory management.
    
    
    // Neighbour-related functions
    void mconstructNeighbours();  // construct lists of neighbours
    
    // other private variables
    VectorXd mrow;   // construct when initialize, for evaluating jacobian

    
public:
    
    // Constructor #1:
    Stratification(int whichcase0 = 0, VectorXd params0 = VectorXd::Zero(1));

    // Labels
    static const int cEq    {1};   // Equation
    static const int cIn    {2};   // Inequality
    static const int cNone  {0};   // Neither; not used
    // Error code
    static const int cErr {-99};
    static const int cFix {1};     // function label is fixed (E/I/N)
    static const int cVary {0};    // function label can vary between E/I 
    

    // -----  All of these must be provided in "setup"  ----- //
    int nvars;   // number of variables
    int nfcns;   // number of function values total (equations + inequalities)
    int nmanif;  // number of manifolds in the stratification
                 // =0 means total number is unknown, so we calculate neighbour
                 // information on the fly
    MatrixXi Llist;   // holds lists of labels (one per row)
    VectorXi FixFcns;  // holds flag saying if function label is fixed or not (nmanif=0)
    
    
    // Initial condition: coordinates, and labels
    VectorXd p0;
    VectorXi L0;
    
    // scalar function to sample on manifold
    double feval(Point&);
    
    // Statistics to save on this stratification
    int nstats;  // number of statistics to compute (0 is don't compute)
    void computestats(Point&, VectorXd&);  // compute statistics of Points
    
    // Optional parameters, helpful for easily changing equations
    int whichcase;   // which case to consider (default = 0)
    VectorXd params;  // Parameters for equations
    MatrixXi edges;  // Edges, if we have a framework
    int nedges;      // total number of edges (not always = no. of functions)
    
    
    
    // -----  These are defined in Stratification.cpp  ----- //
    
    // get number of equations and inequalities, given flags
    int neqns(const VectorXi&);
    int neqns(const VectorXi&, int);  // # of eqns in L <= fcn #k (used in Move::getQinfo)
    int nineq(const VectorXi&);
    
    // evaluate equations, jacobian, and inequalities for given labels
    void evalEqns (const VectorXd&, const VectorXi&, VectorXd&);  // equations
    void evalJacob(const VectorXd&, const VectorXi&, MatrixXd&);  // jacobian
    void evalIneqs(const VectorXd&, const VectorXi&, VectorXd&);  // inequalities
    
    // evaluate equations, jacobian, and inequalities given index i
    double   evalEqns (const VectorXd&, int);   // equations, given index i
    VectorXd evalJacob(const VectorXd&, int);   // jacobian, given index i
    
    // Neighbour functions
    int nbrNIgain(Point&,VectorXi&);  // return number of neighbours, and list of ichange
    int nbrNIlose(Point&,VectorXi&);  // return number of neighbours, and list of ichange
    void nbrLgain(Point&, int, int, VectorXi&);  // return L of ith Lose neighbour
    void nbrLlose(Point&, int, int, VectorXi&);  // return L of ith Gain neighbour
    int getLidx(const VectorXi&);      // converts L to iL
    VectorXi getL(const int iL);       // converts iL to L

};

#endif /* Stratification_hpp */
