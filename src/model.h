
#ifndef MODEL_H_N_BODY
#define MODEL_H_N_BODY

#include <array>

using std::array;

array<double, 2> eval_force_simple(const double * cell_rm, double cell_m, const double * body_pos, double body_m);

#endif // MODEL_H_N_BODY
