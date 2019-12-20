//
//  Stratification.cpp
//  
//
//  Created by Miranda Holmes-Cerfon on 3/2/17.
//
//

#include "Stratification.hpp"
#include "Point.hpp"


// Constructor with functions, jacobian, and optional flag check
Stratification::Stratification(int whichcase0, VectorXd params0)
: whichcase(whichcase0), params(params0)
{
    setup();    // initialize variables: nvars, nfcns, nmanif, Llist
    mrow = VectorXd::Zero(nvars);  // for calculating jacobian

    if(nmanif > 0) {
        mconstructNeighbours();  // construct neighbour lists and nbr information
    }
}



// Return number of equations in L
int Stratification::neqns(const VectorXi& L) {
    int n = 0;
    for (int i=0; i < nfcns; i++) {
        if(L(i) == cEq) n++;
    }
    return n;
}


// Return number of equations in L up to (and including) function k
// Used in Move::getQinfo
int Stratification::neqns(const VectorXi& L, int k) {
    int n = 0;
    for (int i=0; i <= min(k,nfcns-1); i++) {
        if(L(i) == cEq) n++;
    }
    return n;
}



// Return number of inequalities in L
int Stratification::nineq(const VectorXi& L) {
    int n = 0;
    for (int i=0; i < nfcns; i++) {
        if(L(i) == cIn) n++;
    }
    return n;
}

// Evaluate Equations
void Stratification::evalEqns (const VectorXd& p, const VectorXi& L, VectorXd& vals) {
    int counter = 0;
    for(int i=0; i < nfcns; i++) {
        if(L(i) == cEq) {
            vals(counter) = eqs(i,p);
            counter++;
        }
    }
}

// Evaluate Inequalities
void Stratification::evalIneqs (const VectorXd& p, const VectorXi& L, VectorXd& vals) {
    int counter = 0;
    for(int i=0; i < nfcns; i++) {
        if(L(i) == cIn) {
            vals(counter) = eqs(i,p);
            counter++;
        }
    }
}

// Evaluate Jacobian
void Stratification::evalJacob (const VectorXd& p, const VectorXi& L, MatrixXd& jacmtx) {
    int counter_row = 0;
    for(int i=0; i < nfcns; i++) {
        if(L(i) == cEq) {
            jac(i, p, mrow);  // evaluate row of jacobian
            jacmtx.row(counter_row) = mrow;
            counter_row++;
        }
    }
}

// Evaluate Equations, given index i
double Stratification::evalEqns (const VectorXd& p, int ind) {
    return eqs(ind, p);
}

// Evaluate Jacobian, given index i
VectorXd Stratification::evalJacob(const VectorXd& p, int ind) {
    jac(ind, p, mrow);
    return mrow;
}




// ----------------------------- //
//      Neighbour functions      //
// ----------------------------- //


// return number of neighbours, and list of ichange
int Stratification::nbrNIgain(Point& x,VectorXi& igain){
    if(nmanif > 0) {
        int iL = x.iL;  // index of current labels
        int ngain = mnbrN_gain(iL);
        if(ngain > 0) {
            igain =  mnbrIchange_gain.row(iL).head(ngain);
        }
        return ngain;
    }
    // nmanif = 0
    else {
        int ngain = 0;
        igain = VectorXi::Zero(nfcns);  //igain.resize(nfcns);
        for(int i=0;i<nfcns;i++) {
            if(x.L(i) == cEq && FixFcns(i) == cVary) {
                igain(ngain) = i;
                ngain++;
            }
        }
        /* removed sept 5, since we never access values beyond ngain anyways
         if(ngain > 0) {
            igain = igain.head(ngain).eval();
        }*/
        return ngain;
    }
}


// return number of neighbours, and list of ichange
int Stratification::nbrNIlose(Point& x,VectorXi& ilose){
    if(nmanif > 0) {
        int iL = x.iL;  // index of current labels
        int nlose = mnbrN_lose(iL);
        if(nlose > 0) {
            ilose =  mnbrIchange_lose.row(iL).head(nlose);
        }
        return nlose;
    }
    // nmanif = 0
    else {
        int nlose = 0;
        ilose = VectorXi::Zero(nfcns); //ilose.resize(nfcns);
        for(int i=0;i<nfcns;i++) {
            if(x.L(i) == cIn && FixFcns(i) == cVary) {
                ilose(nlose) = i;
                nlose++;
            }
        }
        /* removed sept 5, since we never access values beyond ngain anyways
        if(nlose > 0) {
            ilose = ilose.head(nlose).eval();
        }
         */
        return nlose;
    }
}


// return labels of given neighbour, and index of changed equation, type gain
void Stratification::nbrLgain(Point& x, int whichnbr, int ichange, VectorXi& newL) {
    if(nmanif > 0) {
        int idx = mnbrIL_gain(x.iL,whichnbr);
        newL =  Llist.row(idx);
    }
    // nmanif = 0
    else {
        newL = x.L;
        newL(ichange) = cIn;
    }
}

// same as above, for lose
void Stratification::nbrLlose(Point& x, int whichnbr, int ichange, VectorXi& newL) {
    if(nmanif > 0) {
        int idx = mnbrIL_lose(x.iL,whichnbr);
        newL =  Llist.row(idx);
    }
    // nmanif = 0
    else {
        newL = x.L;
        newL(ichange) = cEq;
    }
}


 // converts L to iL
int Stratification::getLidx(const VectorXi& L) {
    if(nmanif > 0) {
        VectorXi Lt;
        for(int iL=0; iL < nmanif; iL++){
            Lt = Llist.row(iL);
            if((L-Lt).array().abs().any()==0) {  // L==Lt
                return iL;
            }
        }
    }
    return cErr;  // error; L not found in list
}

// converts iL to L
VectorXi Stratification::getL(const int iL) {
    if(nmanif > 0) {
        return Llist.row(iL);
    }
    else {
        VectorXi err(1);
        err(0) = cErr;
        return err;
    }
}



// -----------------------------------
//      Construct neighbour lists
// -----------------------------------
// Only for case where have entire list of manifolds
// Assumes that manifolds are connected if they differ by exactly
// 1 equation
void Stratification::mconstructNeighbours(){
        
    VectorXi L1;
    VectorXi L2;
    int n1, n2;
    
    int cErr = cErr;
    
    // containers to hold data about neighbours
    mnbrN_gain.resize(nmanif);
    mnbrN_lose.resize(nmanif);
    mnbrIchange_gain.resize(nmanif,nmanif);
    mnbrIchange_lose.resize(nmanif,nmanif);
    mnbrIL_gain.resize(nmanif,nmanif);
    mnbrIL_lose.resize(nmanif,nmanif);
    mnbrIchange_gain.fill(cErr);  // initialize to err so don't accidentally use
    mnbrIchange_lose.fill(cErr);
    mnbrIL_gain.fill(cErr);
    mnbrIL_lose.fill(cErr);
    mnbrN_gain.fill(0);
    mnbrN_lose.fill(0);
    
    // loop through manifolds and check for adjacency
    for(int i=0; i< nmanif; i++) {
        for(int j=i+1; j<nmanif; j++) {
            L1 = Llist.row(i);
            L2 = Llist.row(j);
            n1 = neqns(L1);
            n2 = neqns(L2);
            

            // L1-->L2 could be gain dim, L2-->L1 could be lose dim
            if(n1 == n2 + 1) {
                int ndiff = 0;
                int ichange = cErr;
                for(int k=0; k<nfcns; k++) {
                    if(L1(k)==cEq) {
                        if(L2(k)!=cEq) {
                            ndiff++;
                            ichange = k;
                        }
                    }
                    if(ndiff > 1) {
                        break;
                    }
                }
                // L1, L2 differ by exactly 1 equation
                if(ndiff == 1) {
                    mnbrIL_gain(i,mnbrN_gain(i)) = j;
                    mnbrIL_lose(j,mnbrN_lose(j)) = i;
                    mnbrIchange_gain(i,mnbrN_gain(i)) = ichange;
                    mnbrIchange_lose(j,mnbrN_lose(j)) = ichange;
                    mnbrN_gain(i)++;
                    mnbrN_lose(j)++;
                }
            }
            
            // L2-->L1 could be gain dim, L1-->L2 could be lose dim
            if(n2 == n1 + 1) {
                int ndiff = 0;
                int ichange = cErr;
                for(int k=0; k<nfcns; k++) {
                    if(L2(k)==cEq) {
                        if(L1(k)!=cEq) {
                            ndiff++;
                            ichange = k;
                        }
                    }
                    if(ndiff > 1) {
                        break;
                    }
                }
                
                // L1, L2 differ by exactly 1 equation
                if(ndiff == 1) {
                    mnbrIL_lose(i,mnbrN_lose(i)) = j;
                    mnbrIL_gain(j,mnbrN_gain(j)) = i;
                    mnbrIchange_lose(i,mnbrN_lose(i)) = ichange;
                    mnbrIchange_gain(j,mnbrN_gain(j)) = ichange;
                    mnbrN_lose(i)++;
                    mnbrN_gain(j)++;
                }
            }
            
        }  // end loop through j
    }  // end loop through i
}








