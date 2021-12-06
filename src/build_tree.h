

#ifndef _BUILD_TREE_H_N_BODY_
#define _BUILD_TREE_H_N_BODY_

#include <vector>
#include <utility>
#include <array>

#include "body.h"
#include "tree.h"
#include "mpi_types.h"

using std::vector; using std::array;
using std::pair;


typedef vector<pair<array<double, 2>, array<double, 2> > > bound_vec;

void build_simple_tree(const vector<Body > & bodies, Tree & tree, int rank);

void build_tree(const vector<Body > & bodies, const bound_vec & bounds, const bound_vec & other_bounds,
                const vector<pair<int, bool> > & partners, Tree & tree, int rank);

vector<Cell *> construct_received_trees(const vector<MPICell> & recv_cells);




#endif // _BUILD_TREE_H_N_BODY_
