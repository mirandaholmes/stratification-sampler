//
//  Move.cpp
//  
//
//  Created by Miranda Holmes-Cerfon on 10/11/17.
//
//

#include "Move.hpp"

// Constructor
Move::Move(Point* px0, Point* py0)
: px(px0), py(py0) {     }



// Copy movetype data from forward move to *this, but reverse it.
void Move::reverse(const Move& forward) {
    // Set movetype
    if(forward.movetype == cDimSame) {
        movetype = cDimSame;
    }
    else if(forward.movetype == cDimGain) {
        movetype = cDimLose;
        ichange = forward.ichange;
        sides = forward.sides;
    }
    else if(forward.movetype == cDimLose) {
        movetype = cDimGain;
        ichange = forward.ichange;
        sides = forward.sides;
    }
}


// Compute acceptance probability
double Move::axy(void) {
    if(iftouchJxy==cYes && iftouchVxy==cYes && iftouchLxy==cYes && iftouchFxy==cYes) {
        return mjxy*mvxy*mlxy*mfxy;
    }
    else {
        return NAN;
    }
}


void Move::printDensities(void) {
    cout << "Vxy = " << mvxy;
    cout << ", Lxy = " << mlxy;
    cout << ", Jxy = " << mjxy;
    cout << ", fxy = " << mfxy;
    cout << ", Axy = " << axy();
    cout << endl;
}



// =====================================================================
//         Return internal vectors/matrices
// =====================================================================

// Return gradq
VectorXd Move::gradq(void) {
    getQinfo();
    return mgradq;
}

// Return gradqt
VectorXd Move::gradqt(void) {
    getQinfo();
    return mgradqt;
}

// return ugradqt (& calculate it if it doesn't yet exist)
VectorXd Move::ugradqt(void) {
    getQinfo();
    return mugradqt;
}

// return v (hopefully it has been calculated)
VectorXd Move::v(void) {
    if(iftouchV == cYes) return mv;
    else {
        cout << "ERROR in Move::v" << endl;
        VectorXd err(1); err(0) = NAN;
        return err;
    }
}

// return const pointer to Qperp (& calculate it if it doesn't yet exist)
const MatrixXd* Move::TperpQ(void) {
    getQinfo();
    return &mTperpQ;
}

// return const pointer to Tmove (& calculate it if it doesn't yet exist)
const MatrixXd* Move::Tmove(void) {
    getQinfo();
    return &mTmove;
}

// return const pointer to Vperp (& calculate it if it doesn't yet exist)
const MatrixXd* Move::Vperp(void) {
    getVperp();
    return &mVperp;
}



// =====================================================================
//         Set a move (L, movetype, ichange, sides)
// =====================================================================

// Give it L, movetype, ichange.
// Calculates sides.
void Move::setmove(VectorXi& newL, int newmovetype, int newichange) {
    py->setL(newL);
    movetype = newmovetype;
    ichange = newichange;
    sides = 0;
    // set sides
    if(movetype == cDimGain) {
        if(newL(ichange) == Stratification::cIn) sides = 1;
        else sides = 2;
    }
    if(movetype == cDimLose) {
        if((px->L)(ichange) == Stratification::cIn) sides = 1;
        else sides = 2;
    }
}



// Set a move by hand by giving movetype, ichange.
// Calculates L and sets py->L
// Must have called px->getNbrs() first.
// Note this only works if there is only one other manifold with that particular ichange;
// it won't work if there are different possibilities with different inequalities.
// However, this function is mainly intended to use during debugging.
void Move::setmove_hand(int newmovetype,int newichange) {
    // set movetype, ichange
    movetype = newmovetype;
    ichange = newichange;
    
    // movetype == Same: easy to set L
    if(movetype == cDimSame) {
        py->setL(px->L);
        return;
    }
    
    if(px->iftouchNbrs == cNo) {
        cout << "WARNING: Move::setmove: you haven't calculated the neighbours yet." << endl;
        cout << "y.L is not set" << endl;
        return;
    }
    
    // movetype == Gain or == Lose: find approxpriate neighbour
    if(movetype == cDimGain) {
        int ngain = px->ngain;
        VectorXi nbrI = px->ichange_gain;  // list of ichange
        for(int i=0; i<ngain; i++) {
            if(nbrI(i) == ichange) {
                VectorXi newL(px->nfcns);
                px->pstrat->nbrLgain(*px,i,ichange,newL);
                py->setL(newL);
                return;
            }
        }
    }
    if(movetype == cDimLose) {
        int nlose = px->nlose;
        for(int i=0; i<nlose; i++) {
            if(px->nbrIlose(i) == ichange) {
                VectorXi newL(px->nfcns);
                int whichnbr = px->nbrOlose(i);
                px->pstrat->nbrLlose(*px,whichnbr,ichange,newL);
                py->setL(newL);
                return;
            }
        }
    }
    cout << "WARNING: Move::setmove: No neighbour was found for "
    << "movetype==" << movetype << ", ichange == " << ichange << endl;
    cout << "y.L is not set" << endl;
    return;
}



// =====================================================================
//         Compute internal vectors/matrices
// =====================================================================


// Compute tangent space information depending on changed constraint:
// -- Lose: mgradq, mgradqt, mugradqt, mTperpQ, mQmove, mTmove
// -- Gain: mgradq, mgradqt, mugradqt, mTperpQ, mQmove, mTmove
// -- Same: mQmove, mTmove
void Move::getQinfo(void){
    
    if(iftouchQinfo == cNo) {
        px->getTangentSpace();
        
        // Same
        if(movetype == cDimSame) {
            mQmove = px->Q;
            mTmove = px->T;
            mgradq.resize(1); mgradq << NAN;
            mgradqt.resize(1); mgradqt << NAN;
            mugradqt.resize(1); mugradqt << NAN;
            mTperpQ.resize(1,1); mTperpQ << NAN;
        }
        
        // Gain
        if(movetype == cDimGain) {
            // construct Qmove first
            int neqns = px->neqns - 1;
            int nvars = px->nvars;
            if(neqns > 0) {
                mQmove.resize(nvars,neqns);
                int k = px->pstrat->neqns(px->L,ichange);  // which col to drop (start at 1)
                if(k == 1) {       // remove first equation
                    mQmove = px->Q.rightCols(neqns);
                }
                else if(k == neqns+1) {
                    mQmove = px->Q.leftCols(neqns);
                }
                else {
                    mQmove.leftCols(k-1) = px->Q.leftCols(k-1);
                    mQmove.rightCols(neqns-k+1) = px->Q.rightCols(neqns-k + 1);
                }
            }
            else { // don't need Qmove if neqns == 0, since it doesn't even exist
                mQmove.resize(1,1); mQmove << NAN;
            }
            
            // construct Tmove, via QR decomposition
            if(neqns > 0) {
                ColPivHouseholderQR<MatrixXd> qr(mQmove);
                MatrixXd wholeQ = qr.householderQ();
                mTmove = wholeQ.block(0,neqns,nvars,nvars-neqns);
            }
            else {  // neqns == 0
                mTmove = MatrixXd::Identity(nvars, nvars);
            }
            
            // construct mgradq, mgradqt, mugradqt
            mgradq = px->pstrat->evalJacob(px->p, ichange);
            mgradqt = mTmove*mTmove.transpose()*mgradq;  // project gradq onto tan space for moving
            mugradqt = mgradqt / mgradqt.norm();
            
            // *** NOTE: could try using x.N as the projection: Pperp = x.N*x.N.transpose().
            // *** projection matrix is then (eye - Pperp), not mTmove*mTmove.transpose().
            // Tmove is x.T, ugradqt.
            // *** NO! need orth decomposition of Qmove, since N includes gradq.
            //     so still need to do QR decomposition.
            // ******* should check x.t*ugradqt = 0 *******

            // construct mTperpQ
            //getUperp(mugradqt,mTperpQ);
            // ****** CHECK -- I think TperpQ = x.T   *********
            //  *** if so, then no need for QR decomposition here.
            
            mTperpQ = px->T;
        }
        
        // Lose
        if(movetype == cDimLose) {
            // construct mQmove, mTmove
            mQmove = px->Q;
            mTmove = px->T;
            
            // construct mgradq, mgraduqt, mugradqt
            mgradq = px->pstrat->evalJacob(px->p, ichange);
            mgradqt = px->T*px->T.transpose()*mgradq;  // project gradq onto tan space for moving
            mugradqt = mgradqt / mgradqt.norm();
            
            // construct mTperpQ
            getUperp(mgradqt,mTperpQ);
        }
        
        iftouchQinfo = cYes;
    }
}

// get matrix orthogonal to V
void Move::getVperp(void) {
    if(iftouchVperp == cNo) {
        if(px->dim >= 2 && movetype == cDimLose && iftouchV == cYes) {
            mVperp.resize(px->nvars,px->dim - 1);
            getUperp(mv,mVperp);     // tangent space orthogonal to v
        }
        if(iftouchV == cYes) iftouchVperp = cYes;
    }
}


// Compute orth basis of tangent space orthogonal to given vector u
// Assumes that u is in tangent space already. If not, algorithm below is incorrect.
// Note: it doesn't check if u is in tangent space
void Move::getUperp(const VectorXd& u, MatrixXd& Uperp) {
    
    int nvars = px->nvars;
    int neqnsmove = px->neqns;
    if(movetype == cDimGain) neqnsmove = neqnsmove - 1;
    int dim = nvars - neqnsmove;
    
    // if we are on a point or a line return NAN
    if(dim <= 1) {
        Uperp.resize(1,1);
        Uperp(0,0) = NAN;
        return;
    }
    
    px->getTangentSpace();

    // project tangent space at x to space orthogonal to u
    MatrixXd Uproj = (MatrixXd::Identity(nvars,nvars) - u*u.transpose()/u.squaredNorm() );
    if(neqnsmove > 0) { // if no eqns, tangent space is just identity
        Uproj = Uproj * mTmove;
    }
    
    // compute QR decomposition of Uproj
    ColPivHouseholderQR<MatrixXd> qr(Uproj);
    MatrixXd wholeQ = qr.householderQ();
    
    // block is (start row, start col, #rows, #cols)
    Uperp = wholeQ.block(0,0,nvars,nvars-neqnsmove-1);  // orthog basis for Qxv
}



// Compute v given x, y, movetype, ichange, alpha
// Doesn't check for errors, so make sure these variables are already set.
void Move::getv(void) {
    px->getTangentSpace();
    getQinfo();
    VectorXd v(px->nvars);

    if(movetype == cDimSame) {
        if(px->dim == 0) v.setZero();
        else if(px->neqns == 0)  v = (py->p - px->p);
        else v = px->T*px->T.transpose() * (py->p - px->p);

        alpha = NAN;
    }
    if(movetype == cDimGain) {
        v = mTmove*mTmove.transpose() * (py->p - px->p); // compute v
        alpha = NAN;
    }
    if(movetype == cDimLose) {
        v = px->T*px->T.transpose() * (py->p - px->p);
        alpha = v.norm()*v.dot(mugradqt)/abs(v.dot(mugradqt));  // signed norm of v (sgn = dir of ugradqt)
        v = v/v.norm();
    }
    setv(v);
}

// Set v
void Move::setv(VectorXd v) {
    mv.resize(px->nvars);
    mv = v;
    iftouchV = cYes;
}




// =====================================================================
//     getJac: get Jacobian for projection from x -> y in takestep
// =====================================================================
// needs v, x.p, y.p to be set already
// sets jxy
// must be consistent with algorithm used in takestep()
void Move::getJac(void) {
    
    double jxy;   // value of jacobian
    
    if(movetype == cDimSame) {
        jxy = 1.0;
        // it's not actually 1, but it doesn't matter since the ratio jyx/jxy = 1,
        // so the actual value is never used
    }
    else if(movetype == cDimGain) {
        // get tangent spaces, if necessary
        px->getTangentSpace();
        py->getTangentSpace();
        getQinfo();
        
        // compute determinant
        jxy = abs((mTmove.transpose()*py->T).determinant());
    }
    else if(movetype == cDimLose) {
        if(py->dim == 0) { // return right away if we're on a point
            jxy = 1.0;
            setJxy(jxy);
            return;
        }
        // get tangent spaces, if necessary
        px->getTangentSpace();
        py->getTangentSpace();
        getVperp();
        
        // compute determinant
        jxy = abs((mVperp.transpose()*py->T).determinant()) / pow(fabs(alpha),py->dim);
    }
    
    setJxy(jxy);
}




// =====================================================================
//     takestep: take a step on manifold and project back
// =====================================================================
// v must be set
// sets y.p, alpha
// returns flag indicating if newton succeeded or not
int Move::takestep(double tol, double maxIter){
    
    px->getTangentSpace();
    getQinfo();  // needed for mQmove
    int nvars = px->nvars;
    alpha = NAN;

    
    // ------------  takestep: Same or Gain  ------------ //
    if(movetype == cDimSame || movetype == cDimGain) {
        
        int neqns = px->neqns;
        if(movetype == cDimGain) neqns = neqns - 1;

        VectorXd z0;
        z0 = px->p + mv;  // initial point for projection
        
        if(movetype == cDimSame) {
            py->copyStratNbrs(*px);
        }
        
        // there's at least one equation, so we need to project
        if(neqns > 0) {
            // construct labels
            VectorXi L;    // labels of x (for moving)
            if(movetype == cDimSame) L = px->L;        // use x's labels
            else if(movetype == cDimGain) L = py->L;   // use y's labels
            
            // Run Newton's method to project back to manifold
            MatrixXd temp(neqns, nvars);   // needed to calculate jacobian
            VectorXd a(neqns),da(neqns);   // current best solution, and recent change
            VectorXd q(neqns);             // current value of function
            MatrixXd jq(neqns,neqns);      // current value of jacobian
            double ferr = NAN;             // holds error of function
            
            // loop until get close enough to a solution
            a.setZero();
            for(int niter=0; niter < maxIter; niter++) {
                fsame(z0,mQmove,L,a,q);      // sets q
                ferr = q.cwiseAbs().maxCoeff();  // maximum error, pointwise
                
                // check if converged
                if(ferr < tol) {
                    py->setp(z0 + mQmove*a);  // Set proposal point
                    return cNewtonSuccess;
                }
                
                // If not converged, calculate jacobian and update solution
                dfsame(z0,mQmove,L,a,jq,temp);    // sets jq
                da = jq.partialPivLu().solve(-q);  // faster, but less accurate
                //da = jq.fullPivLu().solve(-q);  // slower, more accurate
                a = a + da;
            }
            // maximum number of iterations exceeded; return NaN
            return cNewtonFail;
        }
        // no equations; we're moving in a Euclidean space
        else {
            py->setp(z0);
            return cNewtonSuccess;
        }
    }
    
    // ------------  takestep: Lose  ------------ //
    // Project to lower-dimensional manifold in direction of unit vector mv
    if(movetype == cDimLose) {
        
        int neqns_new = py->neqns;  // should be px->neqns+1
        px->getTangentSpace();
        
        // Construct initial condition for a
        VectorXd a(neqns_new),da(neqns_new);   // current best solution, and recent change
        VectorXd a0(neqns_new);  // initial condition for a
        if(neqns_new > 1) {  // at least 2 equations
            a0.head(neqns_new - 1) = VectorXd::Zero(neqns_new - 1);
        }
        a0(neqns_new - 1) = px->distXk(ichange,mv); // initial guess for distance in direction v
        // must be deterministic
        // originally used 1, but this gave bad results for higher dimensions

        
        // Set up matrices needed for Newton
        MatrixXd temp_jac(neqns_new,nvars);
        VectorXd temp_row(nvars,1);
        MatrixXd Qx_v(nvars,neqns_new);  // ( Q | v )
        
        
        if(neqns_new > 1) {
            Qx_v.leftCols(neqns_new-1) = px->Q;
            Qx_v.col(neqns_new-1) = mv;
        }

        VectorXd q(neqns_new);             // current value of function
        MatrixXd jq(neqns_new,neqns_new);      // current value of jacobian
        double ferr = NAN;             // holds error of function

        
        // Run Newton's method to project back to manifold
        a = a0;
        for(int niter=0; niter < maxIter; niter++) {
            flose(neqns_new,a,q);      // sets q
            ferr = q.cwiseAbs().maxCoeff();  // maximum error, pointwise
            
            // check if converged
            if(ferr < tol) {
                // only keep solution if alpha < 0, to ensure reproducibility for reverse move
                if(a(neqns_new-1) < 0) {
                    alpha = a(neqns_new-1);
                    return cAlphaFail;
                }

                // Newton success, so set new point
                alpha = a(neqns_new-1);
                if(neqns_new > 1) {
                    py->setp(px->p + alpha*mv + (px->Q)*a.head(neqns_new-1));
                } else {
                    py->setp(px->p + alpha*mv);
                }
                return cNewtonSuccess;
            }

            // If not converged, calculate jacobian and update solution
            dflose(neqns_new,a,Qx_v,jq,temp_jac, temp_row);    // sets jq
            da = jq.partialPivLu().solve(-q);  // faster, but less accurate
            //da = jq.fullPivLu().solve(-q);  // slower, more accurate
            a = a + da;
        }
        // maximum number of iterations exceeded; return newton error code
        return cNewtonFail;
    }
    
    return Move::cNewtonFail;
}






              
              
              
              

