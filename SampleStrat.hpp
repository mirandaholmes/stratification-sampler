//
//  SampleStrat.hpp
//  
//
//  Created by Miranda Holmes-Cerfon on 2017-03-14.
//
//

#ifndef SampleStrat_hpp
#define SampleStrat_hpp

#include <stdio.h>
#include <iostream>
#include <iomanip>           // for setw
#include <fstream>           // for writing to files
#include <cmath>
#include <ctime>             // for timing operations
//#include <time.h>            // time
//#include <stdlib.h>          // srand, rand
//#include <random>
#include <Eigen/Dense>
#include "Stratification.hpp"
#include "Point.hpp"
#include "Move.hpp"
#include "Proposals.hpp"


using namespace std;
using namespace Eigen;



// *************************************************** //
// *****   Sample points on a stratification    ****** //
// *************************************************** //
//
// Must provide it with a Stratification object, and a point
// on the stratification.
//
class SampleStrat {
private:
    Stratification *mstrat;   // object containing equations,jacobian, inequalities
    Proposals *mproposals;    // object containing proposal functions and densities
    
    // Parameters for sampling algorithm
    int      mnpts;          // how many points to generate total
    int      mnsave;         // how many points are saved
    int      mprec_pts = 6;     // precision for writing output
    int      mhasinit = SampleStrat::cNoInit; // whether has an initial point on strat
    double   mtol     = 1e-10;    // tolerance for Newton
    int      mmaxIter = 12;       // max. no. iterations for Newton
    
    
    // Save statistics on manifold (optional; must be written by the user)
    ArrayXXd mData;  // statistics to save
    
    
    // Data and statistics from sampling algorithm
    Point    mx0;              // initial point
    MatrixXd mPts;             // list of points sampled
    MatrixXi mLabels;          // list of flags sampled
    double   mtime;            // time it took to sample (in seconds)
    Vector3i mmoves;           // number of each kind of move (gain, lose, same)
    Vector3i mrej;             // total number of rejections
    Vector3i mrej_newton;      // no. of rejections due to Newton's method
    Vector3i mrej_bdy;         // no. of rejections due to stepping outside boundary
    Vector3i mrej_metropolis;  // no. of rejections due to Metropolis-Hastings step
    Vector3i mrej_reverse;     // no. of rejections due to reverse Newton step
    Vector3i mrej_alpha;       // no. of rejections due to alpha < 0
    
    // Statistics for each manifold
    MatrixXi mmovesL;           // number of each kind of move (gain, lose, same)
    MatrixXi mrejL;             // total number of rejections
    MatrixXi mrejL_newton;      // no. of rejections due to Newton's method
    MatrixXi mrejL_bdy;         // no. of rejections due to stepping outside boundary
    MatrixXi mrejL_metropolis;  // no. of rejections due to Metropolis-Hastings step
    MatrixXi mrejL_reverse;     // no. of rejections due to reverse Newton step
    MatrixXi mrejL_alpha;       // no. of rejections due to alpha < 0
    
    
    // Initialize all these matrices
    void initializeStatistics(void);

    
    // Internal Constants
    static const int cDimSame = Move::cDimSame;
    static const int cDimLose = Move::cDimLose;
    static const int cDimGain = Move::cDimGain;
    
     
public:
    // Constructor
    SampleStrat(Stratification *strat,
                Proposals *prop,
                VectorXd &p0,
                VectorXi &L0);
    
    // check that initial point is on the manifold
    int check_onstrat(Point &pt);
    
    // sampling function
    int   sample(int npts=1, int dsave=1);  // sampling algorithm; the main code

    // functions to let you see the internal variables
    int       npts()   const  { return mnpts; }
    int       nsave()  const  { return mnsave; }
    const MatrixXd* pts()   const { return &mPts; }
    const MatrixXi* labels() const { return &mLabels; }
    const ArrayXXd* data() const { return &mData; }
    int nstats()     const { return mstrat->nstats; }
    ArrayXd stats() const { return mData.colwise().sum()/mnsave;  }
   // VectorXd statsstd() const { return square(mData).colwise().sum()/mnsave - mData.colwise().sum().square()/(mnsave*mnsave);}
    static ArrayXd std(ArrayXXd& a) { return (a.square().colwise().sum()/a.rows() - a.colwise().sum().square()/(a.rows()*a.rows())).sqrt(); }
    ArrayXd statsstd(void) { return std(mData); }
    
    // functions to set internal variables
    void setTol    (double tol)   { mtol = tol; }
    void setPrec   (int prec)     { mprec_pts = prec; }
    void setNiter  (int niter)    { mmaxIter = niter; }
    
    // write points and/or flags to file
    int write_pts(string ptsfile, int dwrite=1);
    int write_labels(string flagfile, int dwrite=1);
    int write_labels_idx(string flagfile, int dwrite=1);  // writes indices, not labels
    int write_data(string datafile, int dwrite=1);
    //int write_data(const char* datafile, int dwrite=1);
    
    // Print out statistics after sampling
    void print_rejections(void);
    void print_rejections_manifold(void);
    
    
    // For debugging
    int ifdebug = SampleStrat::cNo;   // 0 for no, 1 for yes
    void setIfDebug(int i) { ifdebug = i; }
    
    
    // Flags
    // success/fail flags
    static const int cSuccess { 0 };
    static const int cFail    { 1 };
    // whether or not have initial condition
    static const int cHasInit   { 0 };
    static const int cNoInit    { 1 };
    // yes/no flags
    static const int cYes       { 0 };
    static const int cNo        { 1 };

    
    // For debugging
    void here(int);
    
    // TO DO:
    // -- save total number of rejections (so can call program multiple times and save
    //         previous statistics)
    //

};






    
#endif /* SampleStrat_hpp */
