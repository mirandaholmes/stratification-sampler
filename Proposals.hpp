//
//  Proposals.hpp
//  
//
//  Created by Miranda Holmes-Cerfon on 10/12/17.
//
//

#ifndef Proposals_hpp
#define Proposals_hpp

#include <stdio.h>
#include <Eigen/Dense>
#include <cmath>
#include <random>
#include "Point.hpp"
#include "Move.hpp"
#include "Stratification.hpp"


using namespace std;
using namespace Eigen;


class Proposals {
private:
    
    // Parameters
    double msig;         // step size in tangent space
    double msigbdy;      // boundary distance parameter
    double msigtan;     // step size in tangent space for gain/lose proposals
    double mlamgain;     // total prob of gaining a dimension
    double mlamlose;  // minimum prob of proposing move "same"
    
    
    // Random number generators and variables
    mt19937                       mrng;        // defines the generator
    int                           mseed;       // seed for rng
    normal_distribution<double>   mZdist;      // standard random numbers
    uniform_real_distribution<double> mUdist;  // uniform on [0,1]
    void         initializeRandom(void);       // initialize random number generators
    
    
public:
    // Constructor #1: all parameters
    Proposals(double sig,
              double sigbdy,
              double sigtan,
              double lamgain,
              double lamlose,
              int seed=0);
    
    
    // Proposal step v, and calculate density
    void v_propose(Move&);
    void v_density(Move&);
    
    // Proposal labels & movetype, and calculate density
    void lam_propose(Move&);
    void lam_density(Move&);
    
    
    // Generate random numbers
    double unifRand(void);   // output a uniform random number on [0,1]
    double normRand(void);   // output a standard normal random number

    
    // Set internal variables
    void setSeed  (int    seed)  { mrng.seed(seed); mseed = seed; }
    int  seed(void) const  { return mseed;  }
    

    // mathematical constants
    constexpr static const double pi = 3.141592653589793238462643383279502884;
    constexpr static const double sqrtpi = 1.77245385090551602729816748334115;
    constexpr static const double sqrt2pi = 2.50662827463100050241576528481105;
    
    // calculate surface area of d-ball
    static double area_dball(int dim) {
        return 2*pow(sqrtpi,dim) / tgamma((double)dim/2.0);
    }
    
    // For debugging
    void here(int);
};

#endif /* Proposals_hpp */



