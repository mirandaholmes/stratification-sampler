//
//  Point.cpp
//  
//
//  Created by Miranda Holmes-Cerfon on 3/2/17.
//
//

#include "Point.hpp"

// for debugging
void Point::here(int i) {
    cout << "here(" << i << ")" << endl;
}



// Constructor #1: Stratification only
Point::Point(Stratification *pS)
: pstrat(pS), nfcns(pstrat->nfcns), nvars(pstrat->nvars) {
    resetFlags();
    p = VectorXd::Zero(nvars);
    L = VectorXi::Zero(nfcns);
}

// Constructor #2; Input = Stratification S, x, L
Point::Point(Stratification *pS, VectorXd& p0, VectorXi& L0)
: pstrat(pS), p(p0), L(L0),
  nfcns(pstrat->nfcns), nvars(pstrat->nvars)  {
      resetFlags();
      resetDimensions();
}

// Constructor #3; Input = Stratification S, L
Point::Point(Stratification *pS, VectorXi& L0)
: pstrat(pS), L(L0),
  nfcns(pstrat->nfcns), nvars(pstrat->nvars) {
      resetFlags();
      resetDimensions();
      p = VectorXd::Zero(nvars);
}


// Reset touch variables
void Point::resetFlags(void) {
    iftouchJac      = cNo;
    iftouchIneq     = cNo;
    iftouchEq       = cNo;
    iftouchTanSpace = cNo;
    iftouchF        = cNo;
    iftouchNbrs     = cNo;
    iftouchTanSpace = cNo;
    iftouchStratNbrs= cNo;
}

// Reset data
void Point::resetDimensions(void) {
    neqns = pstrat->neqns(L);
    nineq = pstrat->nineq(L);
    dim = nvars - neqns;
    EqVals.resize(neqns);
    IneqVals.resize(nineq);
    JacVals.resize(neqns,nvars);
    resetFlags();
}

// Set x-data
void Point::setp(const VectorXd& pnew) {
    p = pnew;
    resetFlags();
}

// Set L-data
void Point::setL(const VectorXi& Lnew) {
    // Check if L changed
     if((L-Lnew).squaredNorm() > 0.25) {    
         L = Lnew;
         resetDimensions();
         // if nmanif > 0, then set iL
         if(pstrat->nmanif > 0) {
             iL = pstrat->getLidx(L);
         }
         else {
             iL = NAN;
         }
         return;
     }
    // resetFlags(); // commented out, since this is done in resetDimensions.
    // If L doesn't change, don't need to reset flags. 
}


// Get values of equations
void Point::getEqns(void){
    if(iftouchEq == cNo) {
        pstrat->evalEqns(p,L,EqVals);
        iftouchEq = cYes;
    }
}

// Get values of Jacobian
void Point::getJacobian(void){
    if(iftouchJac == cNo) {
        pstrat->evalJacob(p,L,JacVals);
        iftouchJac = cYes;
    }
}

// Get values of inequalities
void Point::getIneqs(void) {
    if(iftouchIneq == cNo) {
        pstrat->evalIneqs(p,L,IneqVals);
        iftouchIneq = cYes;
    }
}

void Point::getF(void) {
    if(iftouchF == cNo) {
        f = pstrat->feval(*this);
        iftouchF = cYes;
    }
}

// Get tangent space
void Point::getTangentSpace(void) {
    if(iftouchTanSpace == cNo) {
        if(neqns > 0 && (dim > 0)) {
            if(iftouchJac == cNo) {
                getJacobian();
            }
            Q = JacVals.transpose();
            ColPivHouseholderQR<MatrixXd> qr(Q);
            MatrixXd wholeQ = qr.householderQ();
            N = wholeQ.block(0,0,nvars,neqns); // block is (start row, start col, #rows, #cols)
            T = wholeQ.block(0,neqns,nvars,nvars-neqns);
        } else if (neqns == 0){
            T = MatrixXd::Identity(nvars,nvars);
        } else if (dim == 0) {
            Q = MatrixXd::Identity(nvars, neqns);
            N = MatrixXd::Identity(nvars, neqns);
        }
        iftouchTanSpace = cYes;
    }
}


// Copy y's data into x (assuming from the same stratification)
void Point::copy(const Point& y) {
    
    // copy coordinates & labels
    setp(y.p);
    setL(y.L);
    
    // Do this after setp, setL since they reset the flags
    iftouchJac      = y.iftouchJac;
    iftouchEq       = y.iftouchEq;
    iftouchIneq     = y.iftouchIneq;
    iftouchTanSpace = y.iftouchTanSpace;
    iftouchF        = y.iftouchF;
    iftouchNbrs     = y.iftouchNbrs;
    iftouchStratNbrs= y.iftouchStratNbrs;
    
    
    // Update values of equations, inequalities, jacobian, tangent space if these exist
    if(iftouchEq == cYes) {
        EqVals = y.EqVals;
    }
    if(iftouchIneq == cYes) {
        IneqVals = y.IneqVals;
    }
    if(iftouchJac == cYes) {
        JacVals = y.JacVals;
    }
    if(iftouchF == cYes) {
        f = y.f;
    }
    if(iftouchTanSpace == cYes) {
        Q = y.Q;
        T = y.T;
        N = y.N;
    }
    if(iftouchNbrs == cYes) {
        ngain = y.ngain;
        nlose = y.nlose;
        nbrIlose = y.nbrIlose;
        nbrOlose = y.nbrOlose;
        nbrDlose = y.nbrDlose;
    }
    if(iftouchStratNbrs == cYes) {
        ngain = y.ngain;
        nlose = y.nlose;
        ichange_lose = y.ichange_lose;
        ichange_gain = y.ichange_gain;
    }
}

// copy just strat neighbours into x
void Point::copyStratNbrs(const Point& y) {
    iftouchStratNbrs= y.iftouchStratNbrs;
    if(iftouchStratNbrs == cYes) {
        ngain = y.ngain;
        nlose = y.nlose;
        ichange_lose = y.ichange_lose;
        ichange_gain = y.ichange_gain;
    }
}



// psuedodeterminant of matrix A, with tolerance mtol for zero eigenvalues
double Point::pseudodet(const MatrixXd& A, double tol)  {
    VectorXd s = A.jacobiSvd().singularValues();
    double pdet = 1;
    for (int i=0; i<s.size(); i++){
        if(abs(s(i)) > tol) {  
            pdet *= abs(s(i));
        }
    }
    return pdet;
}


// ------------------------------
//      Neighbour functions
// ------------------------------


// Estimate signed distance to boundary, in direction of given (unit) vector v
// Estimate as
//             w = -q(x) / |gradq(x) . v|
//
// Best guess is y = x + w * v
//
double Point::distXk(int k, const VectorXd& v) {
    double fq = pstrat->evalEqns(this->p,k); // value of equation
    VectorXd gradqk  = pstrat->evalJacob(this->p, k);  // gradient of constraint
    double w = -fq / gradqk.dot(v);  // distance w
    return w;
}


// Estimate signed distance to boundary in direction of Proj(grad(function k)),
// where Proj(v) is the projection onto the tangent space at x.
// Estimate as in distXbdy(k,v), with v = gradqk.
// Code repeated here to avoid calculating gradqk again.
// Best estimate for point is
//           y = x + w * Proj(gradqk(x))/|Proj(gradqk(x))|
//
double Point::distXk(int k) {
    getTangentSpace();
    double w;
    double fq = pstrat->evalEqns(this->p,k); // value of equation
    VectorXd gradqk = pstrat->evalJacob(this->p, k);  // gradient of constraint
    if(dim > 0) {
        w = fabs(fq) / (T.transpose()*gradqk).norm();  // ** changed dec 17: used to be -fq **
        //VectorXd projqk = T*T.transpose()*gradqk;  // gradient, projected onto tangent space
        //w = -fq / ( gradqk.dot(projqk/projqk.norm()) );
    }
    else w = NAN;
    
    return w;
}


// Compute number of lose & gain neighbours, and the indices of neighbours
// For lose, considers those within distance tolbdy of boundary
// Distance has a sign, but it gets all neighbours where absolute value
// of distance is within threshhold tolbdy. 
// Default value (if none is provided) is numerical infinity;
// see http://en.cppreference.com/w/cpp/types/numeric_limits/infinity
//
void Point::getNbrs(double tolbdy) {
    if(iftouchNbrs == cYes) {
        return;
    }
    
    // compute # neighbours, lists of ichange
    if(iftouchStratNbrs == cNo) {
        ngain = pstrat->nbrNIgain(*this,ichange_gain);
        nlose = pstrat->nbrNIlose(*this,ichange_lose);
        iftouchStratNbrs = cYes;
    }

    // If there are no lose neighbours, return immediately
    if(nlose==0) {
        iftouchNbrs = cYes;
        return;
    }

    // For each possible Lose neighbour, estimate distance w to boundary
    // Estimate solution y as y = x + w Proj(gradq(x)) / |Proj(gradq(x))|
    // Linearize -->  w = -q(x) / |gradq(x)|
    nbrOlose.resize(nlose);    // "index" (order) of neighbour with correct distance
    nbrIlose.resize(nlose);   // "ichange" for each neighbour within correct dist
    nbrDlose.resize(nlose);    // "w" (dist) for each neighbour within correct dist
    int nw = 0;    // total number within correct distance
    for(int i=0; i < nlose; i++) {
        int k = ichange_lose(i);
        double w = distXk(k);
        if(fabs(w) < tolbdy) {      // Check if dist is below threshhold
            nbrOlose(nw) = i;      // order of "lose" neighbour which changed
            nbrIlose(nw) = ichange_lose(i);  // label which changed
            nbrDlose(nw) = w;
            nw++;
        }
    }
    
    // Store variables
    nlose = nw;
    nbrOlose = nbrOlose.head(nw).eval();
    nbrIlose = nbrIlose.head(nw).eval();
    nbrDlose = nbrDlose.head(nw).eval();
    
    
    // We're done. Change flag.
    iftouchNbrs = cYes;
    
    
    /* Removed aug 19
     
    // For each lose neighbour, calculate ngain = no. of gain nbrs
    // *** this is pretty inefficient; might think about indexing neighbours by iL in Stratification ***
    nbrLose_ngain = VectorXi::Zero(nw);   // reset this vector and create it again
    VectorXi nbrL(nfcns);
    for(int whichnbr=0; whichnbr < nw; whichnbr++) {
        pstrat->nbrLlose(*this, whichnbr, nbrL);    // changes nbrL  *** inefficient to require nbrL!
        nbrLose_ngain(whichnbr) = pstrat->nbrNgain(nbrL);
    }
    */
    
}














