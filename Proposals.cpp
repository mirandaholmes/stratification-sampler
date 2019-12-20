//
//  Proposals.cpp
//  
//
//  Created by Miranda Holmes-Cerfon on 10/12/17.
//
//
//  Simple2: test anisotropic v-proposals. 
//
//
//

#include "Proposals.hpp"





int whichcaseGain = 0;  // type of Gain proposal; all work
int whichcaseLose = 0;  // type of Lose proposal; all work

double param = 0;  // extra parameter, if needed


// Constructor #1: all parameters
Proposals::Proposals(double sig,
                     double sigbdy,
                     double sigtan,
                     double lamgain,
                     double lamlose,
                     int seed)
: msig(sig), msigbdy(sigbdy), msigtan(sigtan),
  mlamgain(lamgain), mlamlose(lamlose),
  mseed(seed), mZdist(0.,1.), mUdist(0.0,1.0)
{
    initializeRandom();     // initialize random number generators
}




// for debugging
void Proposals::here(int i) {
    cout << "here(" << i << ")" << endl;
}




// ====================================================
//           Random number generators
// ====================================================


// Initialize random number generator
void  Proposals::initializeRandom(void) {
    if(mseed==0) {
        mseed = random_device()(); // choose seed randomly based on device configurations
    }
    mrng.seed(mseed);    // initialize generator for random number stream
}


// Output a uniform random number on [0,1]
double Proposals::unifRand(void) {
    return mUdist(mrng);
}


// Output a standard normal random number 
double Proposals::normRand(void) {
    return mZdist(mrng);
}




// ====================================================
//           v-proposal
// ====================================================


// Propose a vector v and save as mymove.mv
void Proposals::v_propose(Move& mymove) {
    
    // extract variables for ease of writing
    int nvars = mymove.px->nvars;
    int xdim = mymove.px->dim;
    
    // compute tangent space and gradq info
    mymove.px->getTangentSpace();
    mymove.getQinfo();
    VectorXd v(nvars);
    
    /*
    //debug
    cout << "x.L = " << mymove.px->L.transpose() << endl;
    cout << "x.p = " << mymove.px->p.transpose() << endl;
    cout << "x->JacVals = " << mymove.px->JacVals << endl;
    cout << " dx->T = " << mymove.px->T << endl;
    cout << "nvars = " << nvars << ", xdim = " << xdim << endl;
    cout << "tmove = " << mymove.Tmove() << endl;
    cout << "gradq = " << mymove.gradq() << endl;
    cout << "gradqt = " << mymove.gradqt() << endl;
     */
    
    // Movetype == Same
    if(mymove.movetype == Move::cDimSame) {
        if(xdim == 0) {
            v.setZero();;
        }
        else {
            VectorXd r(xdim);  // holds steps in each tangent space direction
            for (int i=0; i<xdim; i++) r(i) = msig * normRand();
            const MatrixXd* Tmove = mymove.Tmove();
            v = (*Tmove)*r;   // random step in tangent space
        }
    }

    // Movetype == Gain
    if(mymove.movetype == Move::cDimGain) {
        
        // -----  for debugging (simple)  -----
        if(whichcaseGain == 1) {
            VectorXd r(xdim+1);  // holds steps in each tangent space direction
            for (int i=0; i<xdim+1; i++) r(i) = msig * normRand();
            const MatrixXd* Tmove = mymove.Tmove();;
            v = (*Tmove)*r;   // random step in tangent space
        }
        
        
        // ----  actual (desired) proposal  -----
        if(whichcaseGain == 0) {
            // chose step in direction of constraint: uniform
            double yorth; // amount move in direction gradq
            if(mymove.sides == 1) yorth = msigbdy * unifRand();  // uniform on [0,sigbdy]
            else {
                yorth = msigbdy * (2*unifRand() - 1);  // uniform on [-sigbdy,sigbdy]
            }
            v = yorth * mymove.ugradqt();
            
            // chose steps in direction orthogonal to constraint
            // distrn is normal, variance sigtan*yorth
            if(xdim >= 1) {
                VectorXd r(xdim);
                for (int i=0; i<xdim; i++) r(i) = msigtan * yorth * normRand(); //
                const MatrixXd* TperpQ = mymove.TperpQ();
                v = v + (*TperpQ)*r;
            }
        }
        
    }
    
    // Movetype == Lose
    if(mymove.movetype == Move::cDimLose) {
        
        // simplest case
        // it works
        if(whichcaseLose == 1) {
            VectorXd r(xdim);  // holds steps in each tangent space direction
            for (int i=0; i<xdim; i++) r(i) = msig * normRand();
            const MatrixXd* Tmove = mymove.Tmove();
            v = (*Tmove)*r;   // random step in tangent space
        }
        
        // same as simplest case, calculated differently
        // it works
        if(whichcaseLose == 2) {
            v = msig * normRand() * (-mymove.ugradqt());
            if(xdim >= 2) {
                VectorXd r(xdim-1);
                for (int i=0; i<xdim-1; i++) r(i) = msig * normRand();
                const MatrixXd* TperpQ = mymove.TperpQ();
                v = v + (*TperpQ)*r;
            }
        }
        
        
        // ----  actual (desired) proposal  -----
        if(whichcaseLose == 0) {
            
            // choose step in direction of -gradq
            v = -mymove.ugradqt();
            if(mymove.sides == 2) {
                // value of function #ichange
                double fq = mymove.px->pstrat->evalEqns(mymove.px->p,mymove.ichange);
                if(fq < 0) v = -v;
            }
            
            
            // choose steps in orthogonal directions
            if(xdim >= 2) {
                VectorXd r(xdim-1);
                for (int i=0; i<xdim-1; i++) r(i) = msigtan * normRand();
                const MatrixXd* TperpQ = mymove.TperpQ();
                v = v + (*TperpQ)*r;
                
                
                //debug
                if(0) {
                    cout << "=============  proposing: ===========" << endl;
                    cout <<" x.p = " << mymove.px->p.transpose() << endl;
                    cout << "-ugradt = " << -mymove.ugradqt().transpose() << endl;
                    cout << "TperpQ = \n" << (*TperpQ) << endl;
                    cout <<" TperpQ*ugradt = " <<(*TperpQ).transpose()*-mymove.ugradqt() << endl;
                    cout << "r = " << r.transpose() << endl;
                    cout << "v = " << v.transpose() << endl;
                }
                //end debug
            }
        }

        // normalize v
        v = v / v.norm();
    }
    
    mymove.setv(v);
}


// ====================================================
//           v-density
// ====================================================

// Calculate density for a given move
// x,v must be set
void Proposals::v_density(Move& mymove) {
    
    double dens;
    
    // extract variables for ease of writing
    const int nvars = mymove.px->nvars;
    const int xdim = mymove.px->dim;
    const int xeqns = mymove.px->neqns;
    
    VectorXd v = mymove.v();
    
    // Movetype == Same
    if(mymove.movetype == Move::cDimSame) {
        if(xdim == 0)  dens = 1.0;   // on a point
        else dens = exp(-0.5*v.squaredNorm()/(msig*msig)) / pow(sqrt2pi*msig,xdim);
    }
    
    // Movetype == Gain
    if(mymove.movetype == Move::cDimGain) {
        
        mymove.px->getTangentSpace();
        mymove.getQinfo();
        
        
        if(whichcaseGain == 1 ) {
            dens = exp(-0.5*v.squaredNorm()/(msig*msig)) / pow(sqrt2pi*msig,xdim+1);
        }
        
        // ----  actual (desired) proposal  -----
        if(whichcaseGain == 0) {
            // get step in direction gradq
            double yorth = v.dot(mymove.ugradqt()); // amount in direction gradq
            double len;  // length of interval for uniform distribution
            
            // reject (set dens=0) if it's too far, not within indicator function
            if(mymove.sides == 1) {
                len = msigbdy;
                if(yorth > msigbdy || yorth < 0) {
                    dens = 0.0;
                    mymove.setVxy(dens);
                    return;
                }
            }
            if(mymove.sides == 2) {
                len = 2*msigbdy;
                if(yorth > msigbdy || yorth < -msigbdy) {
                    dens = 0.0;
                    mymove.setVxy(dens);
                    return;
                }
            }
            
            // Density in direction gradqt (when move is accessible)
            dens = 1.0/len;
            
            // Density in orthogonal directions
            if(xdim >= 1) {
                // find amount in each direction orthogonal to gradq (ie in x.T)
                VectorXd r(xdim);
                const MatrixXd *TperpQ = mymove.TperpQ();
                r = (*TperpQ).transpose()*v; // ** prev px->T
                dens *= exp(-0.5*r.squaredNorm() / (msigtan*msigtan*yorth*yorth))
                / pow(sqrt2pi*msigtan*abs(yorth),xdim);
            }
        }
    }
    
    // Movetype == Lose
    if(mymove.movetype == Move::cDimLose) {
        
        if(whichcaseLose == 1) {
            dens = 1/area_dball(xdim);
        }
        
        
        if(whichcaseLose == 2) {
            dens = 1/area_dball(xdim);
        }
        
        
        // -----  actual (desired) proposal  -----
        if(whichcaseLose == 0) {
           if(xdim <= 1) {
               dens = 1.0;
            }
            else {     // xdim >= 2
                mymove.px->getTangentSpace();
                mymove.getQinfo();
                
                double yq = v.dot(-mymove.ugradqt()); // amount in direction -ugradqt
                if(mymove.sides == 1 && yq < 0) {  // reject if it's not in same dir as -ugradqt
                    dens = 0.0;
                    mymove.setVxy(dens);
                    return;
                }
                if(mymove.sides==2) {
                    double fq = mymove.px->pstrat->evalEqns(mymove.px->p,mymove.ichange);
                    if((yq < 0 && fq > 0) || (yq > 0 && fq < 0)) {
                        dens = 0.0;
                        mymove.setVxy(dens);
                        return;
                    }
                }
                
                // vector actually proposed (unnormalized version)
                VectorXd vorig(nvars);
                vorig = v / abs(yq); // scale so abs(component) in direction ugradqt = 1
                
                // extract amount in directions orthogonal to gradq
                VectorXd r(xdim-1);
                const MatrixXd *TperpQ = mymove.TperpQ();
                r = (*TperpQ).transpose()*vorig;
                
                dens = exp(-0.5*r.squaredNorm() / (msigtan*msigtan)) / pow(sqrt2pi*msigtan,xdim-1);

                // get jacobian factor
                double jacv = 1.0/pow(abs(yq),mymove.py->dim+1);
                dens *= jacv;
                
                
                if(mymove.sides == 2) dens = dens;  // can go in either direction: +-ugradqt
                
                //debug
                if(0) {
                    cout << "......." << endl;
                    cout << "density:" << endl;
                    cout << "r = " << r.transpose() << endl;
                    cout << "vorig = " << vorig.transpose() << endl;
                    cout << "vorig* -ugradqt = " << vorig.dot(-mymove.ugradqt()) << endl;
                    //end debug
                }
            }
        }
        
    }
    
    mymove.setVxy(dens);
}

    



// ====================================================
//           lam-proposal
// ====================================================


// Propose a movetype & label change
// sets mymove.movetype, mymove.py->L, mymove.ichange
void Proposals::lam_propose(Move& mymove) {
    
    // ------  Get necessary neighbour info  ----- //
    // Get neighbour info for x
    mymove.px->getNbrs(msigbdy);

    int nlose = mymove.px->nlose;   // number of lose neighbours
    int ngain = mymove.px->ngain;   // number of gain neighbours

    
    // Total probs of gain, lose, same (must sum to 1), to be calculated
    double pgain, plose, psame;
    pgain = mlamgain;
    plose = mlamlose;
    if(ngain == 0)  pgain = 0;
    if(nlose == 0) plose = 0;
    psame = 1 - plose - pgain;

    // ------  Choose a move  ----- //
    // Decide which movetype
    double u = unifRand();
    
    // Move "Same"
    if(u < psame) {
        mymove.setmove(mymove.px->L,Move::cDimSame);
        return;
    }
    // Move "Gain"
    if(u < psame + pgain) {
        // choose a neighbour uniformly and get its data
        int whichnbr = floor(ngain * unifRand());
        int ichange = mymove.px->ichange_gain(whichnbr);
        VectorXi newL; //(mymove.px->nfcns);
        mymove.px->pstrat->nbrLgain(*(mymove.px), whichnbr, ichange, newL); // get newL
        
        // set new variables in movetype
        mymove.setmove(newL, Move::cDimGain, ichange);
        return;
    }
    // Move "Lose"
    else {
        double ulose = unifRand();
        int i = floor(nlose * ulose);
        
        // choose lose neighbour "i"
        int whichnbr = mymove.px->nbrOlose(i);
        int ichange = mymove.px->nbrIlose(i);
        VectorXi newL; //(mymove.px->nfcns);
        mymove.px->pstrat->nbrLlose(*(mymove.px), whichnbr, ichange, newL); // get newL
        
        // set new variables in movetype
        mymove.setmove(newL, Move::cDimLose, ichange);
        return;
    }
}



// ====================================================
//           lam-density
// ====================================================

// Compute density for movetype & label change
void Proposals::lam_density(Move& mymove) {

    // ------  Get necessary neighbour info  ----- //
    // Get neighbour info for x
    mymove.px->getNbrs(msigbdy);

    int nlose = mymove.px->nlose;   // number of lose neighbours
    int ngain = mymove.px->ngain;   // number of gain neighbours

    
    // Total probs of gain, lose, same (must sum to 1), to be calculated
    double pgain, plose, psame;
    pgain = mlamgain;
    plose = mlamlose;
    if(ngain == 0)  pgain = 0;
    if(nlose == 0) plose = 0;
    psame = 1 - plose - pgain;
    
    
    // ------  Calculate density  ----- //
    int movetype = mymove.movetype;
    double dens;
    if(movetype == Move::cDimSame) {
        dens = psame;
    }
    else if(movetype == Move::cDimGain) {
        dens = pgain / ngain;
    }
    else if(movetype == Move::cDimLose) {
        if(nlose == 0) {
            dens = 0.0;  // could happen in the reverse case, if jump too far away from boundary (because of curvature effects)
        }
        
        // make sure this particular neighbour is accessible (i.e. is among neighbours)
        double w = mymove.px->distXk(mymove.ichange);
        if(abs(w) > msigbdy) {  // *** CHECK for what kind of move (one-sided/two-sided)
            // ** w right now is always positive. So must check value of function.
            dens = 0.0;
        }
        else {
            dens = plose / nlose;
        }
    }
    
    // set density
    mymove.setLxy(dens);
}





