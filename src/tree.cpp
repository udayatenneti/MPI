#include "tree.h"

#include <iostream>
#include <sstream>
#include <string>
#include <algorithm>
#include <math.h>

#include "utils.h"
#include "model.h"

using std::cout; using std::endl;
using std::stringstream; using std::string;
using std::copy; using std::fill;
using std::array;



Tree::Tree(const double * min_bounds, const double * max_bounds, double theta, int rank){
    /* parameter for opening criterion */
    m_theta = theta;

    /* create root cell */
    root = new Cell();
    copy(min_bounds, &min_bounds[2], root->min_bounds);
    copy(max_bounds, &max_bounds[2], root->max_bounds);
    root->m = 0;
    fill(root->subcells, &root->subcells[4], nullptr);
}
    
void Tree::insert_body(const Body * b){
    insert_body(root, b);
}

void Tree::insert_body(Cell * cell, const Body * b){
    
    /* if cell does not contain a body and is a leaf */
    if(cell->b == nullptr and cell->subcells[0] == nullptr){ 
        cell->m = b->m; 
        copy(std::begin(b->pos), std::end(b->pos), cell->rm);
        cell->b = b;
    }
    else{
        /* if cell already contains a body, but is not split */
        if(cell->subcells[0] == nullptr){   
            /* split the cell  */
            for(int i = 0; i < 4; i++){
                /* create subcell */
                Cell * subcell = new Cell();

                subcell->m = 0;

                /* find min and max bounds of cell */
                for(int c = 0; c < 2; c++){
                    double half_side = (cell->max_bounds[c] - cell->min_bounds[c]) / 2;
                    int shift = ((i >> c) & 1);
                    subcell->min_bounds[c] = cell->min_bounds[c] + shift * half_side;
                    subcell->max_bounds[c] = cell->max_bounds[c] - !shift * half_side;
                }

                /* set as subcell */
                cell->subcells[i] = subcell;
                
                // insert the old body into new cell  
                if(cell->b != nullptr and coord_in_cell(cell->subcells[i], cell->rm)){
                    insert_body(cell->subcells[i], cell->b);
                    cell->b = nullptr;
                }
            }
        }

        /* insert new body into subcell */
        for(int i = 0; i < 4; i++){
            if(cell->subcells[i] != nullptr and coord_in_cell(cell->subcells[i], b->pos)){
                insert_body(cell->subcells[i], b);
                break;
            }
        }

        /* update mass and center of mass for cell */
        for(int c = 0; c < 2; c++){
            cell->rm[c] = (cell->m * cell->rm[c] + b->m * b->pos[c]) / (cell->m + b->m);
        }
        cell->m += b->m;
    }
}


bool Tree::mac(const Cell * cell, const double * pos){
    /* minimum squared distance between cell and domain */
    double d = dist_rm(cell->rm, pos);

    /* squared volume of cell */
    //double v2 = pow(cell_volume(cell),0.333*2);
    double l = (double) (cell->max_bounds[0] - cell->min_bounds[0]);
    /* opening criterion rule */
    if( (l / d) < m_theta){
        return true; //divide into subcells
    }
    else{
        return false; //calculate force from center of mass
    }
}

array<double, 2> Tree::compute_force(const Body * b){
    return compute_force(root, b);
}

array<double, 2> Tree::compute_force(const Cell * cell, const Body * b){
    if (b->m == -1){
        array<double, 2> ftot = {{0, 0}};
        return ftot;
    }
    /* if cell is non empty and we dont want to open it or it is leaf cell */
    if((cell->m > 0) and (cell->subcells[0] == nullptr) and (cell->b != nullptr)){
        if(cell->b->idx != b->idx){
            /* evaluate force */
            return eval_force_simple(cell->rm, cell->m, b->pos, b->m); //rm should be same as body pos!
        }
        /* is this really needed? */
        else{
            return {{0, 0}}; //same as
        }
    } else if((cell->m > 0) and (cell->b == nullptr) and (cell->subcells[0] != nullptr) and mac(cell, b->pos)){
        //force is calculated using
        return eval_force_simple(cell->rm, cell->m, b->pos, b->m);
    }

    /* accumulate force evaluation from subcells */
    array<double, 2> ftot = {{0, 0}};
    for(int i = 0; i < 4; i++){
        if(cell->subcells[i]!=nullptr){
            /* call force computation recursively */
            array<double, 2> f = compute_force(cell->subcells[i], b);

            /* accumulate */
            for(int c = 0; c < 2; c++){
                ftot[c] += f[c];
            }
        }
        else{
            break;
        }
    }
    return ftot;
}
    
double Tree::mass(){
    return root->m;
}

int Tree::size(bool complete_tree){
    return size(root, complete_tree); 
}

int Tree::size(Cell * c, bool complete_tree){
    int s = 0;
    int i;
    for(i = 0; i < 4; i++){
        if(c->subcells[i] != nullptr){
            s += size(c->subcells[i], complete_tree);
        }
        else{
            break;
        }
    }
    if(i == 0){
        return 1;
    }
    /*if(i == 0 and c->m != 0){
        return 1;
    }*/
    else{
        std::cout << "\rHit incorrectly " << "/" << std::endl;;
        if(complete_tree){
            s++;
        }
        return s; 
    }
}




bool coord_in_cell(Cell * cell, const double * pos){
    for(int c = 0; c < 2; c++){
        /* if coordinate does not lie in cell boundary */
        if(pos[c] < cell->min_bounds[c] or pos[c] > cell->max_bounds[c]){
            return false;
        }
    }
    return true;
}

bool same_cell(Cell * cell, const double * min_bounds, const double * max_bounds){
    for(int c = 0; c < 2; c++){
        /* if coordinate boundaries or not same */
        if(min_bounds[c] != cell->min_bounds[c] or max_bounds[c] != cell->max_bounds[c]){
            return false;
        }
    }
    return true;
}

double cell_volume(const Cell * cell){
    double v = 0;
    double tmp;
    for(int c = 0; c < 2; c++){
       /* coordinate side length */
       tmp = (cell->max_bounds[c] - cell->min_bounds[c]);
       if(v == 0){
            v = tmp;
       }
       else{
            v *= tmp;
       }
    }
    return tmp;
}


bool bounds_in_cell(Cell * cell, const double * min_bounds, const double * max_bounds){
    for(int c = 0; c < 2; c++){
        /* if coordinate boundary does not fit in cell*/
        if(cell->min_bounds[c] > min_bounds[c] or cell->max_bounds[c] < max_bounds[c]){
            return false;
        }
    }
    return true;
}
