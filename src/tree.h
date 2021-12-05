#ifndef _TREE_H_N_BODY_
#define _TREE_H_N_BODY_

#include <string>
#include <vector>
#include <array>

#include "body.h"

using std::vector; using std::array;

struct Cell {
    double min_bounds[2];
    double max_bounds[2];
    double m = 0;   // mass 
    double rm[2]; // center of mass
    Cell * subcells[4]; // eight subcells maximum
    bool inserted = false;
    const Body * b = nullptr;
    bool send_as_leaf = false;
    int list_idx = -1; 
    int parent_idx = -1;
};

class Tree {
public:

    Tree(const double * min_bounds, const double * max_bounds, double theta, int rank);


    void insert_body(const double * pos, double m);
    
    void insert_body(const Body * b);

    void insert_emptycell(const double * min_bounds, const double * max_bounds);

    void insert_cell(const double * min_bounds, const double * max_bounds, double m, const double * rm);

    void insert_cell(Cell * cell);

    void cells_to_send(const double * min_bounds, const double * max_bounds,
                    int min_depth, vector<Cell*> & cells);


    void prune_tree(const double * min_bounds, const double * max_bounds);


    array<double, 2> compute_force(const Body * b);


    std::string to_string(bool fulltree=true) const;

    double mass();

    int size(bool complete_tree=true);
    

private:
    /*
    insert body into subtree
    
    @param c : root cell of the tree
    @param pos : position of body
    @param m : mass of body
    */
    void insert_body(Cell * cell, const double * pos, double m);
    
    /*
    inserts a body into the tree (expands it)
    */
    void insert_body(Cell * c, const Body * b);
    

    bool mac(const Cell * cell, const double * pos);

    array<double, 2> compute_force(const Cell * cell, const Body * b);


    int size(Cell * c, bool complete_tree);

    Cell * root;

    double m_theta;

};

bool coord_in_cell(Cell * cell, const double * pos);

bool same_cell(Cell * cell, const double * min_bounds, const double * max_bounds);

bool bounds_in_cell(Cell * cell, const double * min_bounds, const double * max_bounds);

#endif // _TREE_H_N_BODY_
