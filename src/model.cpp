
#include "model.h"
#include "utils.h"
#include <cmath>

#include <iostream>

using std::array;

array<double, 2> eval_force_simple(const double * cell_rm, double cell_m, const double * body_pos, double body_m){
    
    double RLIMIT = 0.03;
    //double M = cell_m * body_m;
    array<double, 2> f;
    
    double dist = dist_rm(cell_rm, body_pos);
    double dx = cell_rm[0] - body_pos[0];
    double dy = cell_rm[1] - body_pos[1];

    if (dist < RLIMIT) {
        double scale = RLIMIT / fabs(dist);
        dx = dx * scale;
        dy = dy * scale;
        dist = RLIMIT;
    }


    f[0] = 0.0001 * cell_m * body_m * dx / pow(dist, 3.0);
    f[1] = 0.0001 * cell_m * body_m * dy / pow(dist, 3.0);
    return f;
}
