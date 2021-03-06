#include "utils.h"

#include <math.h>
#include <iostream>
using std::cout; using std::endl;


double bisection(double a, double b, std::function<double (double)> f, double tol, int max_iter){

    double c, fa, fb, fc;
    int n_iter;

    fa = f(a); fb = f(b);

    n_iter = 0;
    while(fabs(b-a) > tol && fa != fb && n_iter < max_iter){
        c = (a + b) / 2;
        fc = f(c);


        if(fc == 0){
            return c; 
        }
        else if(fa * fc > 0){
            a = c;
        }
        else{
            b = c;
        }
        n_iter++;
    }

    return (a + b) / 2;
}

double dist_aabb(const double * a_min, const double * a_max, const double * b_min, const double * b_max){
    double dist = 0;
    double delta;
    for(int c = 0; c < 2; c++){
    
        if(a_min[c] > b_max[c]){
               delta = a_min[c] - b_max[c];
               dist += delta * delta;
        }
        else if(b_min[c] > a_max[c]){
               delta = b_min[c] - a_max[c];
               dist += delta * delta;
        }
    }
    return dist;
}

double dist_rm(const double * cell_rm, const double * body_pos){
    double dist = 0;
    double delta;
    for(int c = 0; c < 2; c++){
        delta = cell_rm[c] - body_pos[c];
        dist += delta * delta;
    }
    dist = pow(dist, 0.5);
    return dist;
}
