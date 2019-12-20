//
//  Point.hpp
//  
//
//  Created by Miranda Holmes-Cerfon on 3/2/17.
//
// This class holds data associated with a point on a stratification.
// The main data is
//
//       pstrat = pointer to a Stratification object
//       p = list of coordinates, a vector in R^d;
//       L = list of "labels", identifying which manifold it is on.
//
// The labels are used by the Stratification class so Point does not need to know
// how these are implemented; the implementation can be changed in the Stratification
// class. 
//
// Point has functions to calculate values of the equations, Jacobian, and
// inequalities, which all depend on the values of p:
//
//       getEqns(), getJacobian(), getIneqs() .
//
// These functions are implemented by calling the relevant functions in Manifold.
// The values are stored in public member vectors/matrices
//
//     EqVals, JacVals, IneqVals
//
//
// The function
//
//       getTangentSpace()
//
// calculates orthonormal bases for the tangent space and normal space at the point.
// These depend on the Jacobian. If the Jacobian doesn't exist, then getTangentSpace
// first calculates the Jacobian. The values are stored in public member matrices
//
//        Q, T, N.
//
//
// All variables are public for convenience and efficiency (e.g. when sending a vector
// as an argument to a function; we don't want to copy the vector.)
//
// HOWEVER, you SHOULD NOT SET the value of p,L directly -- rather, you should
// call the access functions:
//
//        setp(newp), setL(newL) .
//
// This ensures that the dimensions of the data vectors/matrices, and internal flags,
// are set correctly.
//
// There are four internal flags which say whether certain data has been calculated:
//
//        iftouchEqs, iftouchJac, iftouchIneq, iftouchTanSpace, iftouchNbrs .
//
// Each time the given data is calculated these flags are set to "yes".
// If you ask to calculate data and it has already been calculated, nothing happens.
// The function getTangentSpace() uses these flags because it depends on the Jacobian
// but if it has been calculated, then it can skip this (expensive) step.
//
// The flags are reset to "no" if you change p (by calling setp(..).)
//
//
// If you have a point already defined and you want to efficiently "erase" its data without
// erasing the point itself and without setting any new data, you should
// call the function:
//
//        resetFlags().
//
// This sets the "touch" variables to the default values indicating the point contains
// no data. Then, even if the values of the equations, Jacobian, etc have been set
// with previous data, any functions which call Point will treat it as if it has no data.
//
//
//
//    ** COMMENT on neighbour functions (added Oct 16 2017)
//
// -- don't set them directly; rather use getNbrs(bdytol)
//
//
//


#ifndef Point_hpp
#define Point_hpp

#include <stdio.h>
#include <Eigen/Dense>
#include <limits>              // for representing infinity
#include "Stratification.hpp"



using namespace std;
using namespace Eigen;



class Point {
public:
    
    // Constructor #1; Input = Stratification S
    Point(Stratification*);
    
    // Constructor #2; Input = Stratification S, p, L
    Point(Stratification*, VectorXd&, VectorXi&);
    
    // Constructor #3; Input = Stratification S, L
    Point(Stratification*, VectorXi&);
    
    
    // Variables defining exactly where it is on stratification
    Stratification *pstrat;  // which stratification it belongs to
    const int nfcns;    // total number of functions in the list (length of L)
    const int nvars;    // total number of variables (length of x)
    VectorXd p;         // coordinates of point
    VectorXi L;         // "Labels" identifying manifold
    int iL;             // index of manifold
    int neqns;          // number of equations
    int nineq;          // number of inequalities
    int dim;            // dimension of manifold (currently = nvars - neqns)
    
    // Neighbour variables
    int ngain;          // number of neighbours, gain
    int nlose;          // number of neighbours, lose
    VectorXi ichange_lose;  // ichange for "lose" neighbours, all possible
    VectorXi ichange_gain;  // ichange for "gain" neighbours, all possible
    VectorXi nbrIlose;  // ichange for "lose" neighbours within sigbdy
    VectorXi nbrOlose;  // order of "lose" nbrs within sigbdy (out of all strat nbrs)
    VectorXd nbrDlose;  // distances to boundary for neighbours within sigbdy
    // aug 19 VectorXi nbrLose_ngain;  // number of "gain" for each lose neighbour
    
    
    // Values we calculate from the stratification
    VectorXd EqVals;      // values of equations
    VectorXd IneqVals;    // values of inequalities
    MatrixXd JacVals;     // jacobian
    MatrixXd Q;           // jacobian^T
    MatrixXd T;           // tangent space, orthonormal basis
    MatrixXd N;           // normal space, orthonormal basis
    double f;           // value of feval at x
    
    // Evaluate equations, inequalities, jacobian, feval
    void getEqns(void);
    void getJacobian(void);
    void getIneqs(void);
    void getF(void);
    
    // Calculate tangent space
    void getTangentSpace(void);
    
    // calculate pseudo-determinant of a matrix A with given tolerance 
    // Used for clusters, to evaluate f(x), the function to sample
    double pseudodet(const MatrixXd& A, double tol=1e-6);
    
    // Neighbour functions
    void getNbrs(double tolbdy = inf);  // compute neighour information
    double distXk(int);  // est. distance from x to q_k=0, in direction Proj(grad(q_k))
    double distXk(int,const VectorXd&);  // est. distance from x to q_k=0, in direction v
    
    
    // Use these functions to manipulate internal values
    void setp(const VectorXd&);  // change x-data
    void setL(const VectorXi&);  // change L-data
    void copy(const Point&);     // copy y's data into x
    void copyStratNbrs(const Point&);  // copy just strat neighbours into x
    
    // change internal values
    void resetFlags(void);        // reset "touch" variables
    void resetDimensions(void);   // reset vectors/matrices/neqns/nineq, depending on L
    
    // binary variables; =Yes if variables have been set, =No otherwise
    int iftouchJac;
    int iftouchIneq;
    int iftouchEq;
    int iftouchTanSpace;
    int iftouchF;
    int iftouchNbrs;
    int iftouchStratNbrs;
    
    
    // Flags and constants
    static const int cYes {1};
    static const int cNo  {0};
    constexpr static const double inf = std::numeric_limits<double>::infinity();
    
    
    // For debugging
    void here(int);
};


#endif /* Point_hpp */
