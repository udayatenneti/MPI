/* input/output */
#include <iostream>

/* math */
#include <cmath>

/* srand */
#include <stdlib.h>

/* STL types */
#include <vector>
#include <array>

/* MPI library */
#include <mpi.h>

/* Custom MPI structs */
#include "mpi_types.h"

/* Body class */
#include "body.h"

/* Reading and writing */
#include "readwrite.h"

/* Tree building */
#include "tree.h"
#include "build_tree.h"

/* Input parsing */
#include "inputparser.h"

/* Body generator */


int main(int argc, char * argv[]){
    
    int size, rank, tmax, N;
    std::vector<Body> bodies;
    std::vector<Body> allBodies;
    std::vector<pair<double, int> > splits;
    vector<double> comp_time;
    double dt, min[2], max[2];
    double start_time, stop_time;
    InputParser ip;
    bool overwrite;
    
    /* Initialize MPI */
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if(!ip.parse(argc, argv)){
        if(rank == 0){
            std::cout << "Error: Invalid command line input" << std::endl;
            print_usage(argc, argv);
        }
        return -1; 
    }

    /* Initialize custom MPI structures */
    init_mpi_types();

    /* Read bodies from file */
    auto p = read_bodies(ip.in_file().c_str(), MPI_COMM_WORLD);
    bodies = p.first;
    allBodies = p.second.first;
    N = p.second.second;
    std::cout << "\rRank: " << rank << " my bodies: "<< bodies.size() << " all bodies: "<< allBodies.size() <<  "/" << std::endl;

    /* Write initial positions to file */
    overwrite = true;

    
    tmax = ip.n_steps(); // number of time steps
    dt = ip.time_step(); // time step
    
    if(ip.clock_run()){
        MPI_Barrier(MPI_COMM_WORLD);
        start_time = MPI_Wtime();
    }

    for(int t = 0; t < tmax; t++){
        
        /* Domain composition and transfer of bodies */
        for(int c = 0; c < 2; c++){
            min[c] = 0;
            max[c] = 4;
        }

        /* Build the local tree */
        //Tree tree(min, max, ip.bh_approx_constant());
        //build_tree(bodies, bounds, other_bounds, partners, tree, rank);
        Tree tree(min, max, ip.bh_approx_constant(), rank);
        build_simple_tree(allBodies, tree, rank);

        /* Compute forces */
        std::vector<array<double, 2> > forces;
        for(Body & b : bodies){
            /* time the computation */
            
            double start_time = MPI_Wtime();
            
            array<double, 2> f = tree.compute_force(&b);

            /* update the workload for the body */
            b.work = MPI_Wtime() - start_time;
            forces.push_back(f);
        }
        
        /* Update positions */
        for (unsigned int i = 0; i < bodies.size(); i++) {
            Body & b = bodies[i];
            if (b.m != -1){
                for(int c = 0; c < 2; c++){
                    double a = forces[i][c] / b.m;
                    b.pos[c] = b.pos[c] + b.vel[c] * dt + 0.5 * a * pow(dt, 2);
                    b.vel[c] = b.vel[c] + a * dt;
                    if (b.pos[c] < 0 || b.pos[c] > 4){
                        b.m = -1;
                    }
                }
            }
        }

        //share bodies with other processes...
        //updated allBodies = bodies + bodiesFromOthers
        int counts[size];
        int displacements[size];
        int per_process = N / size;
        int modulo = N % size;
        int displacement = 0;
        for (int i=0; i < size; i++){
            int offset = 0;
            if (modulo > 0){
                if ( i < modulo ){
                    offset = 1;
                }
            }
            int this_size = per_process + offset;
            counts[i] = this_size;
            displacements[i] = displacement;
            displacement += this_size;
        }
        /*if (rank==0){
            std::cout << "\rsize: " << N << std::endl;
            for (int i=0; i < size; i++){
                std::cout << "\rcounts: " << counts[i] << " displacements: " << displacements[i] << std::endl; 
            }
        }*/
        MPI_Barrier(MPI_COMM_WORLD);
        allBodies.clear();
        allBodies.resize(N);
        
        MPI_Allgatherv(&bodies.front(), bodies.size(), mpi_body_type, &allBodies.front(), counts, displacements, mpi_body_type, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        //MPI_Allgatherv(my_values, my_values_count, MPI_INT, buffer, counts, displacements, MPI_INT, MPI_COMM_WORLD);
        
        /*if (rank==0){
            for (int i = 0; i < allBodies.size(); i++){
                Body & b = allBodies[i];
                std::cout << "\r idx: " << b.idx << std::endl;
            }
        }*/

        MPI_Barrier(MPI_COMM_WORLD);
        /* Output */
        /* Print time step to stdout */
        if(rank == 0 and ip.verbose()){
            std::cout << "\rTime step: " << t + 1 << "/" << tmax;
            if(t == tmax - 1){
                std::cout << std::endl;
            }
        }



        if(ip.sampling_interval() == 1 or (t % ip.sampling_interval() == 0 and t != 0)){
            stop_time = 0;
            start_time = 0;
            /* Stop the time */
            if(ip.clock_run()){
                MPI_Barrier(MPI_COMM_WORLD);
                stop_time = MPI_Wtime(); 
            }

            /* Write tree */
            if(rank == 0 and ip.write_tree()){
                write_tree(ip.out_tree_file().c_str(), tree, true, overwrite);
            }
            
            /* Write tree size*/
            if(rank == 0 and ip.write_tree_size()){
                write_to_file(ip.out_tree_size_file().c_str(), tree.size(), overwrite);
            }

            /* Write positions */
            if(ip.write_positions()){
                //write_bodies(ip.out_file().c_str(), bodies, MPI_COMM_WORLD, false);
            }

            /* Write running time */
            if(ip.clock_run()){
                if(rank == 0){
                    write_to_file(ip.out_time_file().c_str(), stop_time - start_time, overwrite);
                }
                
            }

            if(overwrite){
                overwrite = false;
            }

            if(ip.clock_run()){
                // start the time
                MPI_Barrier(MPI_COMM_WORLD);
                start_time = MPI_Wtime();        
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    write_bodies(ip.out_file().c_str(), bodies, MPI_COMM_WORLD, true);
    /* Finalize */
    free_mpi_types();
    MPI_Finalize();
   
    if(ip.write_summary()){
        if(rank == 0){
            write_summary(ip, N, size);
        }
    }
}

