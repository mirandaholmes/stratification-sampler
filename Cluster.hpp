//
//  Cluster.cpp
//  
//
//  Created by Miranda Holmes-Cerfon on 2019-07-19.
//

#ifndef Cluster_hpp
#define Cluster_hpp

#include <stdio.h>
#include "Point.hpp"


// ====================================================
//           Cluster functions
// ====================================================

namespace Cluster {
    
    // Vibrational factor
    double getH(Point& x) {
        return 1./x.pseudodet(x.JacVals / 2.);
    }
    // **NOTE: doesn't account for extra constraints right now.
    // Assumes all functions are bond constraints.
    
    // cosine of angle between particles i,j,k
    double costh(VectorXd& x, int i, int j, int k, int dim) {
        VectorXd p1(dim);
        VectorXd p2(dim);
        for (int c=0; c<dim; c++) {
            p1(c) = x(dim*i+c) - x(dim*j+c);
            p2(c) = x(dim*k+c) - x(dim*j+c);
        }
        return p1.dot(p2)/p1.norm()/p2.norm();
    }
    
    // distance between particles i,j
    double dist(VectorXd& x, int i, int j, int dim){
        double d = 0;
        for(int k=0;k<dim; k++) {
            d += (x(dim*i+k) - x(dim*j+k)) * (x(dim*i+k) - x(dim*j+k)) ;
        }
        return sqrt(d);
    }
}


#endif /* Cluster_hpp */

