//
//  Move.hpp
//  
//
//  Created by Miranda Holmes-Cerfon on 10/11/17.
//
//
//
//
// **** DESCRIPTION!!!!
//
//
//
//
//
//
//
//
//
//
//
//
//
//

#ifndef Move_hpp
#define Move_hpp

#include <stdio.h>
#include <Eigen/Dense>
#include <cmath>
#include "Point.hpp"

using namespace std;
using namespace Eigen;


class Move {
private:
    
    // Member variables which are calculated (hence made private)
    VectorXd mgradq;     // gradient of constraint which changed
    VectorXd mgradqt;    // gradient of constraint *in tangent space*
    VectorXd mugradqt;   // unit vector in direction gradq
    MatrixXd mTperpQ;    // orthonormal basis of tangent space normal to gradq
    MatrixXd mTmove;    // orth basis of tangent space for moving, (Gain)
    MatrixXd mQmove;    // constraints for moving (Same or Gain)
    VectorXd mv;        // proposed move in tangent space
    MatrixXd mVperp;    // tangent space orthogonal to v

    
    // Densities
    double mjxy;      // jacobian factor, from projection of step
    double mvxy;      // v-factor
    double mlxy;      // lambda-factor, from changing labels/manifolds
    double mfxy;      // f-factor, from function want to evaluate
    
    
    // binary variables; =Yes if variables have been set, =No otherwise
    int iftouchQinfo = cNo;
    int iftouchJxy = cNo;
    int iftouchVxy = cNo;
    int iftouchLxy = cNo;
    int iftouchFxy = cNo;
    int iftouchV   = cNo;
    int iftouchVperp = cNo;

    
public:
    
    // Constructor
    Move(Point*, Point*);
    
    // Move information (public for convenient access to their data -- be careful!)
    Point *px;       // starting point
    Point *py;       // ending point
    double alpha;    // length of proposed move (movetype==lose)
    int movetype;    // type of move: same, gain, lose
    int ichange;     // which label changed (if any)
    int sides;       // number of sides for a label change (1 if E->I, 2 if E->N)
    
    // Function to copy manifold data from forward move, but reverse it
    void reverse(const Move&);
    
    // Functions to compute derived quantities
    void getQinfo(void);   // compute mgradq, mgradqt, mugradqt, mTperpQ, mQmove, mTmove
    void getUperp(const VectorXd&, MatrixXd& Uperp);  // basis of getQperp, getVperp
    void getVperp(void); // get tangent space orthogonal to V
    void getJac(void);     // compute Jxy given x,y,movetype,ichange,v,alpha
    void getv(void);       // compute v (&alpha) given x,y,movetype,ichange
    void setv(VectorXd);  // set v, given vector
    
    // Function to take a step along manifold and project back
    int takestep(double, double);
    
     
    // Set a move (L, movetype, ichange)
    void setmove(VectorXi&, int, int newichange = -1);
    void setmove_hand(int,int);  // set it by hand; give movetype, ichange
    
    // Set internal density variables
    void setJxy(double d) { mjxy=d; iftouchJxy = cYes; }
    void setVxy(double d) { mvxy=d; iftouchVxy = cYes; }
    void setLxy(double d) { mlxy=d; iftouchLxy = cYes; }
    void setFxy(double d) { mfxy=d; iftouchFxy = cYes; }
    
    // Access internal variables
    double jxy(void) { return mjxy; }
    double vxy(void) { return mvxy; }
    double lxy(void) { return mlxy; }
    double fxy(void) { return mfxy; }
    double axy(void);                   // compute & return product of densities
    void printDensities(void);    // Print out values of all densities
    VectorXd gradqt(void);
    VectorXd ugradqt(void);
    VectorXd gradq(void);
    VectorXd v(void);
    const MatrixXd* TperpQ(void);
    const MatrixXd* Tmove(void);
    const MatrixXd* Vperp(void);
    
    
    // Function and Jacobian for Same / Gain moves
    void fsame(const VectorXd& z0, const MatrixXd& Q, const VectorXi &L, const VectorXd& a, VectorXd &fun)   {
        px->pstrat->evalEqns(z0 + Q*a , L, fun);
    }
    void dfsame(const VectorXd& z0, const MatrixXd& Q, const VectorXi &L, const VectorXd& a, MatrixXd &jac, MatrixXd &temp)  {
        px->pstrat->evalJacob(z0 + Q*a, L, temp);
        jac = temp * Q;
    }
    
    // Function and Jacobian for Lose moves
    void flose(int mneqns, const VectorXd& a, VectorXd &fun)   {
        if(mneqns > 1) {
            px->pstrat->evalEqns(px->p + a(mneqns-1)*mv + (px->Q)*a.head(mneqns-1),
                             px->L, fun);
            fun(mneqns-1) = px->pstrat->evalEqns(px->p + a(mneqns-1)*mv +
                                             (px->Q)*a.head(mneqns-1), ichange);
        } else {
            fun(0) = px->pstrat->evalEqns(px->p + a(mneqns-1)*mv, ichange);
        }
    }
    void dflose(int mneqns, const VectorXd& a, const MatrixXd& Qx_v, MatrixXd& jac, MatrixXd& temp_jac, VectorXd temp_row)  {
        if(mneqns > 1) {
            px->pstrat->evalJacob(px->p + a(mneqns-1)*mv + (px->Q)*a.head(mneqns-1),
                              px->L, temp_jac);
            temp_row = px->pstrat->evalJacob(px->p + a(mneqns-1)*mv + px->Q*a.head(mneqns-1),
                                         ichange);
            temp_jac.row(mneqns-1) = temp_row;
            jac = temp_jac*Qx_v; //
        } else {
            temp_row = px->pstrat->evalJacob(px->p + a(mneqns-1)*mv, ichange);
            jac = temp_row.transpose()*mv;
        }
    }
    
    // Flags
    static const int cYes {1};
    static const int cNo  {0};
    static const int cErr {-99};
    // for movetype
    static const int cDimSame   { 0 };   // move on same-dimensional manifold
    static const int cDimLose   { 1 };   // move to lower-dimensional manifold
    static const int cDimGain   { 2 };   // move to higher-dimensional manifold
    // for takestep
    static const int cNewtonSuccess { 0 };
    static const int cNewtonFail    { 1 };
    static const int cAlphaFail { 2 };
};



#endif /* Move_hpp */
