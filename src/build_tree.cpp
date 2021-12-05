#include "build_tree.h"

#include<math.h>

#include <iostream>
#include <algorithm>
#include <mpi.h>
#include "mpi_types.h"

using std::cout; using std::endl;
using std::copy; using std::fill;
using std::begin; using std::end;



void build_simple_tree(const vector<Body > & bodies, Tree & tree, int rank){

    /* insert bodies residing on this process */
    for(const Body & b : bodies){
        if (b.m != -1){
            tree.insert_body(&b);
        }
    }               
}
