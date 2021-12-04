
#ifndef MODEL_H_N_BODY
#define MODEL_H_N_BODY

#include <array>

using std::array;

array<double, 2> eval_force(const double * r1, double m1, const double * r2, double m2);

array<double, 2> eval_force_simple(const double * cell_rm, double cell_m, const double * body_pos, double body_m);

#endif // MODEL_H_N_BODY
