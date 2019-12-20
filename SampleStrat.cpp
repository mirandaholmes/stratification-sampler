//
//  SampleStrat.cpp
//  
//
//  Created by Miranda Holmes-Cerfon on 2017-03-14.
//
//

#include "SampleStrat.hpp"



// for debugging
void SampleStrat::here(int i) {
    cout << "here(" << i << ")" << endl;
}




// Constructor, with x0, L0, but no f
SampleStrat::SampleStrat(Stratification *strat,
                         Proposals *prop,
                         VectorXd &p0,
                         VectorXi &L0)
: mstrat(strat), mproposals(prop),mx0(mstrat, p0, L0) {
    initializeStatistics();  // set sizes and initialize statistics matrices
    mhasinit = check_onstrat(mx0);  // check initial point is on manifold
}


// Initialize statistics matrices
void SampleStrat::initializeStatistics(void) {
    mrej            << 0,0,0;
    mrej_newton     << 0,0,0;
    mrej_bdy        << 0,0,0;
    mrej_metropolis << 0,0,0;
    mrej_reverse    << 0,0,0;
    mrej_alpha      << 0,0,0;
    mmoves          << 0,0,0;
    
    const int nmanif = mstrat->nmanif;
    if(nmanif > 0) {  // if we have a list of manifolds
        mmovesL.setZero(nmanif,3);
        mrejL.setZero(nmanif,3);
        mrejL_newton.setZero(nmanif,3);
        mrejL_bdy.setZero(nmanif,3);
        mrejL_metropolis.setZero(nmanif,3);
        mrejL_reverse.setZero(nmanif,3);
        mrejL_alpha.setZero(nmanif,3);
    }
}


// check that initial point is on the manifold
int SampleStrat::check_onstrat(Point &pt){
    int neqns = pt.neqns;
    int nineq = pt.nineq;
    // Check that it satisfies the equations
    if(neqns > 0) {
        VectorXd feval(neqns);
        mstrat->evalEqns(pt.p, pt.L, feval);
        
        if(feval.cwiseAbs().maxCoeff() > mtol) {
            cout << "Oops! You picked a point that is not on the stratification: \n"
            << "   it doesn't satisfy the equations." << endl;
            cout << "Functions = " << feval.transpose() << endl;
            return SampleStrat::cNoInit;
        }
    }
    // Check that it satisfies the inequalities
    if(nineq > 0) {
        VectorXd h(nineq);
        mstrat->evalIneqs(pt.p,pt.L,h);
        if(h.minCoeff() < 0-mtol){
            cout << "Oops! You picked a point that is not on the stratification: \n"
            << "   it doesn't satisfy the inequalities." << endl;
            return SampleStrat::cNoInit;
        }
    }
    // Check that the manifold is included in the stratification
    if(pt.pstrat->nmanif > 0) {
        int iL = pt.pstrat->getLidx(pt.L);
        if(iL == Stratification::cErr) {
            cout << "Oops! You picked a point that is not on the stratification: \n"
                 << "   that manifold is not included." << endl;
            return SampleStrat::cNoInit;
        }
    }
    return SampleStrat::cHasInit;
}


// Write points to file, at intervals of dwrite points
int SampleStrat::write_pts(string ptsfile,
                              int dwrite) {
    if(mhasinit == SampleStrat::cNoInit) {
        cout << "write_flags: no data to write." << endl;
        return SampleStrat::cFail;
    }
    
    
    ofstream myfile1 (ptsfile); // SHOULD check for problems
    myfile1.precision(mprec_pts);   // Set precision
    
    // Write points to file
    int i=0;
    while(i<mnsave) {
        myfile1 << mPts.row(i) << endl;
        i += dwrite;
    }
    myfile1.close();
    
    // return
    return SampleStrat::cSuccess;
}
// Write flags to file, at intervals of dwrite points
int SampleStrat::write_labels(string flagfile,
                                int dwrite) {
    if(mhasinit == SampleStrat::cNoInit) {
        cout << "write_labels: no data to write." << endl;
        return SampleStrat::cFail;
    }
    
    ofstream myfile2 (flagfile); // SHOULD check for problems
    int i=0;
    while(i<mnsave) {
        myfile2 << mLabels.row(i) << endl;
        i += dwrite;
    }
    myfile2.close();
    return SampleStrat::cSuccess;
}

// Write flag indices to file, at intervals of dwrite points
int SampleStrat::write_labels_idx(string flagfile,
                             int dwrite) {
    if(mhasinit == SampleStrat::cNoInit) {
        cout << "write_labels_idx: no data to write." << endl;
        return SampleStrat::cFail;
    }
    
    ofstream myfile2 (flagfile); // SHOULD check for problems
    int i=0;
    VectorXd L;
    while(i<mnsave) {
        myfile2 << mstrat->getLidx(mLabels.row(i)) << endl;
        i += dwrite;
    }
    myfile2.close();
    return SampleStrat::cSuccess;
}

// Write data to file, at intervals of dwrite points
int SampleStrat::write_data(string datafile,
                          int dwrite) {
    
    ofstream myfile1 (datafile); // SHOULD check for problems
    myfile1.precision(mprec_pts);   // Set precision
    
    // Write data to file
    int i=0;
    while(i<mnsave) {
        myfile1 << mData.row(i) << endl;
        i += dwrite;
    }
    myfile1.close();
    
    // return
    return 0;
}



// Print rejection statistics to screen
void SampleStrat::print_rejections(void) {
    
    if(mhasinit == SampleStrat::cNoInit) {
        cout << "print_stats: no data to print." << endl;
        return;
    }
    
    mrej << 0,0,0;
    for(int i=0;i<3;i++) {
        mrej(i) = mrej_newton(i) + mrej_bdy(i) + mrej_metropolis(i)
        + mrej_reverse(i) + mrej_alpha(i);
    }
    
    int w  = 12;   // width of boxes in table, for absolute numbers
    int w2 = 12;   // width for rates
    
    cout << "\n=========  Rejection Statistics:  =========" << endl;
    /* cout << "Total number of points : " << mnpts           << endl;
    cout << "-----  Absolute numbers: " << setw(w) << "Total" << " |"
         << setw(w) << "Same"
         << setw(w) << "Lose"
         << setw(w) << "Gain"  << endl;
    cout << "Rejections (All)       : "
         << setw(w) << mrej.sum() << " |"
         << setw(w) << mrej(0)
         << setw(w) << mrej(1)
         << setw(w) << mrej(2)  << endl;
    cout << "Rejections (Newton)    : "
         << setw(w) << mrej_newton.sum() << " |"
         << setw(w) << mrej_newton(0)
         << setw(w) << mrej_newton(1)
         << setw(w) << mrej_newton(2)  << endl;
    cout << "Rejections (Alpha)     : "
         << setw(w) << mrej_alpha.sum() << " |"
         << setw(w) << mrej_alpha(0)
         << setw(w) << mrej_alpha(1)
         << setw(w) << mrej_alpha(2)  << endl;
    cout << "Rejections (Boundary)  : "
         << setw(w) << mrej_bdy.sum() << " |"
         << setw(w) << mrej_bdy(0)
         << setw(w) << mrej_bdy(1)
         << setw(w) << mrej_bdy(2)  << endl;
    cout << "Rejections (Metropolis): "
         << setw(w) << mrej_metropolis.sum() << " |"
         << setw(w) << mrej_metropolis(0)
         << setw(w) << mrej_metropolis(1)
         << setw(w) << mrej_metropolis(2)  << endl;
    cout << "Rejections (Reverse)   : "
         << setw(w) << mrej_reverse.sum() << " |"
         << setw(w) << mrej_reverse(0)
         << setw(w) << mrej_reverse(1)
         << setw(w) << mrej_reverse(2)  << endl;
     */
    
    // Set precision
    streamsize prec_init = std::cout.precision();  // current precision
    cout.precision(3);
    
    // Write rates
    cout << "-----  Rates  ---------:" << setw(w2) << "    Total" << "  |"
         << setw(w2) << "Same"
         << setw(w2) << "Lose"
         << setw(w2) << "Gain"  << endl;
    cout << "Moves                  : "
         << setw(w2) << " " << " |"
         << setw(w2) << (double)mmoves(0)/mnpts
         << setw(w2) << (double)mmoves(1)/mnpts
    << setw(w2) << (double)mmoves(2)/mnpts  << endl;
    cout << "---------------------------------------"
         << "--------------------------------------- " << endl;
    cout << "Rejections (All)       : "
         << setw(w2) << (double)mrej.sum()/mnpts << " |"
         << setw(w2) << (double)mrej(0)/mmoves(0)
         << setw(w2) << (double)mrej(1)/mmoves(1)
         << setw(w2) << (double)mrej(2)/mmoves(2)  << endl;
    cout << "Rejections (Newton)    : "
         << setw(w2) << (double)mrej_newton.sum()/mnpts << " |"
         << setw(w2) << (double)mrej_newton(0)/mmoves(0)
         << setw(w2) << (double)mrej_newton(1)/mmoves(1)
         << setw(w2) << (double)mrej_newton(2)/mmoves(2)  << endl;
    cout << "Rejections (Alpha)     : "
         << setw(w2) << (double)mrej_alpha.sum()/mnpts << " |"
         << setw(w2) << (double)mrej_alpha(0)/mmoves(0)
         << setw(w2) << (double)mrej_alpha(1)/mmoves(1)
         << setw(w2) << (double)mrej_alpha(2)/mmoves(2)  << endl;
    cout << "Rejections (Boundary)  : "
         << setw(w2) << (double)mrej_bdy.sum()/mnpts << " |"
         << setw(w2) << (double)mrej_bdy(0)/mmoves(0)
         << setw(w2) << (double)mrej_bdy(1)/mmoves(1)
         << setw(w2) << (double)mrej_bdy(2)/mmoves(2)  << endl;
    cout << "Rejections (Metropolis): "
         << setw(w2) << (double)mrej_metropolis.sum()/mnpts << " |"
         << setw(w2) << (double)mrej_metropolis(0)/mmoves(0)
         << setw(w2) << (double)mrej_metropolis(1)/mmoves(1)
         << setw(w2) << (double)mrej_metropolis(2)/mmoves(2)  << endl;
    cout << "Rejections (Reverse)   : "
         << setw(w2) << (double)mrej_reverse.sum()/mnpts << " |"
         << setw(w2) << (double)mrej_reverse(0)/mmoves(0)
         << setw(w2) << (double)mrej_reverse(1)/mmoves(1)
         << setw(w2) << (double)mrej_reverse(2)/mmoves(2)  << endl;
    cout << "---------------------------------------"
        << "--------------------------------------- " << endl;
    
    // Return to default precision
    cout.precision(prec_init);
    
    cout << "-----  Timing  -----" << endl;
    cout << "Total time    : " << mtime << " seconds" << endl;
    //cout << "Time per point: " << mtime/mnpts << " seconds" << endl;
    double nn1 = 1e6;
    cout << "You can generate " << nn1 << " points in " << mtime*nn1/(double)mnpts
    << " seconds";
    cout << " or " << mtime*nn1/(double)mnpts/60.0
    << " minutes." << endl;
}


// Print statistics to screen
void SampleStrat::print_rejections_manifold(void) {
    
    if(mhasinit == SampleStrat::cNoInit) {
        cout << "print_stats_manifold: no data to print." << endl;
        return;
    }
    
    // Set precision and table formatting
    streamsize prec_init = std::cout.precision();  // current precision
    cout.precision(3);
    int w2 = 12;   // width of boxes in table
    
    const int nmanif = mstrat->nmanif;
    const int nfcns = mstrat->nfcns;
    VectorXi L(nfcns);
    
    cout << "\n=========  Rejection Statistics by manifold  =========" << endl;
    for (int i=0;i < nmanif; i++) {
        L = mstrat->getL(i);
        
        cout << "---------------------------------------" << endl;
        cout << "   Manifold " << i << ",  neqns = " << mstrat->neqns(L)
             << ", dim = " << mstrat->nvars - mstrat->neqns(L) << endl;
        cout << "---------------------------------------" << endl;
        
        int np = mmovesL.row(i).sum();
        for(int j=0;j<3;j++) {
            mrejL(i,j) = mrejL_newton(i,j) + mrejL_bdy(i,j) + mrejL_metropolis(i,j)
            + mrejL_reverse(i,j) + mrejL_alpha(i,j);
        }
        
        
        // Write rates
        cout << "-----  Rates  -----    :" << setw(w2) << "  Total" << " |"
        << setw(w2) << "Same"
        << setw(w2) << "Lose"
        << setw(w2) << "Gain"  << endl;
        cout << "Moves                  : "
        << setw(w2) << (double)np/mnpts << " |"
        << setw(w2) << (double)mmovesL(i,0)/np
        << setw(w2) << (double)mmovesL(i,1)/np
        << setw(w2) << (double)mmovesL(i,2)/np  << endl;
        //cout << "---------------------------------------"
        //<< "--------------------------------------- " << endl;
        cout << "Rejections (All)       : "
        << setw(w2) << (double)mrejL.row(i).sum()/np << " |"
        << setw(w2) << (double)mrejL(i,0)/mmovesL(i,0)
        << setw(w2) << (double)mrejL(i,1)/mmovesL(i,1)
        << setw(w2) << (double)mrejL(i,2)/mmovesL(i,2)  << endl;
        cout << "Rejections (Newton)    : "
        << setw(w2) << (double)mrejL_newton.row(i).sum()/np << " |"
        << setw(w2) << (double)mrejL_newton(i,0)/mmovesL(i,0)
        << setw(w2) << (double)mrejL_newton(i,1)/mmovesL(i,1)
        << setw(w2) << (double)mrejL_newton(i,2)/mmovesL(i,2)  << endl;
        cout << "Rejections (Alpha)     : "
        << setw(w2) << (double)mrejL_alpha.row(i).sum()/np << " |"
        << setw(w2) << (double)mrejL_alpha(i,0)/mmovesL(i,0)
        << setw(w2) << (double)mrejL_alpha(i,1)/mmovesL(i,1)
        << setw(w2) << (double)mrejL_alpha(i,2)/mmovesL(i,2)  << endl;
        cout << "Rejections (Boundary)  : "
        << setw(w2) << (double)mrejL_bdy.row(i).sum()/np << " |"
        << setw(w2) << (double)mrejL_bdy(i,0)/mmovesL(i,0)
        << setw(w2) << (double)mrejL_bdy(i,1)/mmovesL(i,1)
        << setw(w2) << (double)mrejL_bdy(i,2)/mmovesL(i,2)  << endl;
        cout << "Rejections (Metropolis): "
        << setw(w2) << (double)mrejL_metropolis.row(i).sum()/np << " |"
        << setw(w2) << (double)mrejL_metropolis(i,0)/mmovesL(i,0)
        << setw(w2) << (double)mrejL_metropolis(i,1)/mmovesL(i,1)
        << setw(w2) << (double)mrejL_metropolis(i,2)/mmovesL(i,2)  << endl;
        cout << "Rejections (Reverse)   : "
        << setw(w2) << (double)mrejL_reverse.row(i).sum()/np << " |"
        << setw(w2) << (double)mrejL_reverse(i,0)/mmovesL(i,0)
        << setw(w2) << (double)mrejL_reverse(i,1)/mmovesL(i,1)
        << setw(w2) << (double)mrejL_reverse(i,2)/mmovesL(i,2)  << endl;
        
    }

}




// ====================================================
//     Sampling algorithm: the heart of the class
// ====================================================
int SampleStrat::sample(int npts, int dsave) {
    
    if(mhasinit == SampleStrat::cNoInit) {
        cout << "No initial condition on the stratification " << endl;
        cout << "   Returning to main " << endl;
        return SampleStrat::cFail;
    }
    
    cout << "SampleStrat::sample:\n  Sampling manifold, dsave = " << dsave
         << ", npts = " << npts << endl;
    
    
    // save internal values passed through function interface
    mnpts = npts;
    
    // extract variables for ease of writing
    const int nvars = mstrat->nvars;
    const int nfcns = mstrat->nfcns;
    const int nmanif = mstrat->nmanif;
    
    
    // Initialize vectors and points
    Point x(mstrat);     // start point
    Point y(mstrat);     // end point
    Point xrev(mstrat);  // reverse move // **UPDATE -- need now?
    int iL;               // which manifold we're on (for recording statistics)


    // Other variables needed in the code
    double pacc;                      // total metropolis factor
    int newtonflag, moveflag, newtonflag_rev;
    clock_t start_time, end_time;     // start,end times
    
    
    // Initialize matrix to hold data
    mnsave = floor(npts / dsave);
    mPts   = MatrixXd::Zero(mnsave,nvars);
    mLabels = MatrixXi::Zero(mnsave,nfcns);
    if(mstrat->nstats > 0) {
        mData = MatrixXd::Zero(mnsave,mstrat->nstats);
    }
    VectorXd tempdata(mstrat->nstats);
    int ipt = 0;

    // Initialize points
    x.copy(mx0);
    x.getTangentSpace();
    x.getF();

    
    // *************  Loop through MCMC steps  ***************
    start_time = clock();            // start timing now
    for(int i=0;i<npts;i++){
    
        //debug
        if(ifdebug == cYes) {
            cout << " -----------" << endl;
            cout << "i = " << i << endl;
        }
        // end debug
        

        // save x, flags (save here, since don't know when reject & restart loop)
        if(i != 0 && i%dsave==0) {
            mPts.row(ipt) = x.p;
            mLabels.row(ipt) = x.L;
            if(mstrat->nstats > 0) {   // save statistics
                mstrat->computestats(x,tempdata);
                mData.row(ipt) = tempdata;
                //mstrat->computestats(x,mStats);
            }
            ipt++;
        }
        
        // reset proposal point (more efficient than defining anew each iteration)
        y.resetFlags();
        y.setL(x.L); // **why do this? (in case move is Same?)
        
        // Create a move
        Move mymove(&x, &y);

        
        // =============  Choose manifold to move to (movetype, L)  ==============
        
        mproposals->lam_propose(mymove); // sets mymove.movetype, y.L, mymove.ichange, mymove.sides
        
        //debug
        if(ifdebug == cYes) {
            cout << "movetype = " << mymove.movetype << ", ichange = " << mymove.ichange << ", sides = " << mymove.sides << endl;
        }
        //end debug
    
        
        // save statistics of move
        mmoves(mymove.movetype) += 1;  // total number of moves of each kind
        if(nmanif > 0) {
            iL = mstrat->getLidx(x.L);
            mmovesL(iL,mymove.movetype)++;     // no. of movetype moves, from each manifold
        }
        
        
        // =============  Choose step v in tangent step and project  ==============

        // choose step v
        // sets mymove.mv
        mproposals->v_propose(mymove);

        
        // take step and project
        // sets y.p, alpha (if successful)
        newtonflag = mymove.takestep(mtol,mmaxIter);
        
        /*
        cout << "x = " << x.p.transpose() << endl;
        cout << "L = " << x.L.transpose() << endl;
        cout << "v = " << mymove.v().transpose() << endl;
        cout << "newtonflag = " << newtonflag << endl;
         */
        
        
        // Check for Newton failure
        if(newtonflag==Move::cNewtonFail) {
            mrej_newton(mymove.movetype)++;
            if(nmanif > 0) {
                mrejL_newton(iL,mymove.movetype)++;
            }
            continue;
        }
        // Check for alpha failure
        if(newtonflag == Move::cAlphaFail) {
            mrej_alpha(mymove.movetype)++;
            if(nmanif > 0) {
                mrejL_alpha(iL,mymove.movetype)++;
            }
            continue;
        }

        
        
        // debug
        if(ifdebug == cYes) {
            cout << "forward move: " << endl;
            cout << "  movetype = " << mymove.movetype;
            if(mymove.movetype == 0) cout << ": Same " << endl;
            else if(mymove.movetype == 1) cout << ": Lose "<< endl;
            else if(mymove.movetype == 2) cout << ": Gain "<< endl;
            cout << "  x.p = " << x.p.transpose() << endl;
            cout << "  y.p = " << y.p.transpose() << endl;
            cout << "  alpha = " << mymove.alpha << endl;
        }
        //end debug
        
        
        
        // =============  Inequality check  =============
        if(y.nineq > 0) {
            y.getIneqs();
            if((y.IneqVals).minCoeff() < 0-mtol){
                mrej_bdy(mymove.movetype)++;
                if(nmanif > 0) {
                    mrejL_bdy(iL,mymove.movetype)++;
                }
                continue;
            }
        }
        
        
        // =============  Metropolis step  =============
        
        // Calculate tangent space at y and feval at y
        y.getTangentSpace();
        y.getF();

        
        // Get densities for forward move
        mproposals->v_density(mymove);     // sets vxy
        mproposals->lam_density(mymove);   // sets lxy
        mymove.getJac();                   // sets jxy
        mymove.setFxy(x.f);                // sets fxy
        

        // Create reverse move and get densities
        Move mymove_rev(&y,&x);       // creates reverse move
        mymove_rev.reverse(mymove);   // fills in manifold data by reversing forward move
        mymove_rev.getv();            // get reverse move v
        mproposals->v_density(mymove_rev);     // sets vxy
        mproposals->lam_density(mymove_rev);   // sets lxy
        mymove_rev.getJac();                   // sets jxy
        mymove_rev.setFxy(y.f);                // sets fxy
        
        
        //debug
        if(ifdebug == cYes){
            cout << "----\nmovetype = " << mymove.movetype<< ": " ;
            if(mymove.movetype == 0) cout << "Same " << endl;
            else if(mymove.movetype == 1) cout << "Lose "<< endl;
            else if(mymove.movetype == 2) cout << "Gain "<< endl;
            cout << "Forward: ";
            mymove.printDensities();
            cout << "Reverse: ";
            mymove_rev.printDensities();
            cout << "Acceptance prob: " << mymove_rev.axy() / mymove.axy() << endl;
        }
        //end debug
        
        
        // Compute Metropolis-Hastings probability, and accept or reject
        pacc = mymove_rev.axy() / mymove.axy();  // acceptance probability
        if(pacc > 1) pacc = 1;
        
        if(mproposals->unifRand() > pacc) {
            mrej_metropolis(mymove.movetype)++;
            if(nmanif > 0) {
                mrejL_metropolis(iL,mymove.movetype)++;
            }
            continue;
        }
        
        
        
        // ===========  Reverse projection: recover x = y + v' + Qy*a  ================
        
        // Set up a move for testing reverse projection
        xrev.resetFlags();
        xrev.setL(x.L); // sets labels of point to be found by reverse solver
        Move moverevcheck(&y,&xrev);
        moverevcheck.reverse(mymove);  // copy reversed information from mymove
        moverevcheck.setv(mymove_rev.v());  // set proposal step
        
        
        // Take a step in reverse direction. Sets xrev.p, possibly alpha
        newtonflag_rev = moverevcheck.takestep(mtol,mmaxIter);
        
        
        // debug
        if(ifdebug == cYes) {
            cout << "reverse move: " << endl;
            cout << "  movetype = " << mymove_rev.movetype << ": ";
            if(mymove_rev.movetype == 0) cout << "Same " << endl;
            else if(mymove_rev.movetype == 1) cout << "Lose "<< endl;
            else if(mymove_rev.movetype == 2) cout << "Gain "<< endl;
            cout << "  alpha = " << moverevcheck.alpha << endl;
            cout << "  newtonflag_rev = " << newtonflag_rev  << endl;
            cout << "  vrev  = " << moverevcheck.v().transpose() << endl;
            cout << "  x.p-y.p = " << (x.p-y.p).transpose() << endl;
            cout << "  y.p    = " << y.p.transpose() << endl;
            cout << "  x.p    = " << x.p.transpose() << endl;
            cout << "  xrev.p = " << xrev.p.transpose() << endl;
        }
        //end debug
        
        
        // failed alpha-test (for lose move)
        if(newtonflag_rev == Move::cAlphaFail) {
            mrej_alpha(mymove.movetype)++;
            if(nmanif > 0) {
                mrejL_alpha(iL,mymove.movetype)++;
            }
            continue;
        }
        
        // Newton's method failed, or gave different point
        if(newtonflag_rev == Move::cNewtonFail || (x.p-xrev.p).norm() > 50*sqrt(nvars)*mtol) {
            
            //debug
            //if(newtonflag_rev == Move::cNewtonFail) cout << "NEWTON failed" << endl;
            //else cout << "DX failed" << endl;
            //end debug
            
            mrej_reverse(mymove.movetype)++;
            if(nmanif > 0) {
                mrejL_reverse(iL,mymove.movetype)++;
            }
            continue;
        }
        
        
        // ==============  Accept move, since it passed all tests  ==============
        // Copy y to x. This copies all tangent space and neighbour data too.
        x.copy(y);
        
        
    }    // end loop through MCMC steps
    
    
    // save last value ended at, for restarting
    mx0.copy(x);
    mhasinit = SampleStrat::cHasInit;  // now we can restart at the last point
    
    // save point, if the modular arithmetic tells us to
    if (npts % dsave == 0) {
        mPts.row(ipt) = x.p;
        mLabels.row(ipt) = x.L;
        if(mstrat->nstats > 0) {
            mstrat->computestats(x,tempdata);
            mData.row(ipt) = tempdata;
        }
    }    
    // save time
    mtime = (clock() - start_time ) / (double) CLOCKS_PER_SEC;
    
    return SampleStrat::cSuccess;
    
}    // end sample()





