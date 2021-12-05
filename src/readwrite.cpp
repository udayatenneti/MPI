#include "readwrite.h"

#include <string>
#include <fstream>
#include <sstream>

#include <sys/stat.h>

bool file_exists(const std::string& filename)
{
    struct stat buf;
    if (stat(filename.c_str(), &buf) != -1)
    {
        return true;
    }
    return false;
}

std::pair<std::vector<Body>, std::pair<std::vector<Body>, int>> read_bodies(const char * filename, MPI_Comm comm){
    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    std::vector<Body> bodies;
    std::vector<Body> allbodies;
    /* Read bodies from file */
    if(rank != 0){
        MPI_Status status;
        int x;
        MPI_Recv(&x, 1, MPI_INT, rank - 1, 0,
                MPI_COMM_WORLD, &status);
    }

    std::ifstream infile;
    infile.open(filename);

    std::string line;
    int i = -1;
    int total;
    while (std::getline(infile, line)){
        std::istringstream iss(line);
        double x, y, vx, vy, m;
        
        int num;
        if(i==-1){
            if(!(iss >> total)){ break;}
            i = 0;
            continue;
        }
        if (!(iss >> num >> x >> y >> m >> vx >> vy)) { break; } // error
        if((i % size) == rank){
            bodies.push_back(Body{{x, y}, {vx, vy}, m, 1, num});
        }
        allbodies.push_back(Body{{x, y}, {vx, vy}, m, 1, num});
        i++;
    }
    infile.close();

    if(rank != size - 1){
        int x = 1;
        MPI_Send(&x, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
    }
    std::pair<std::vector<Body>, int> wholeLayer = std::make_pair(allbodies, total);
    return std::make_pair(bodies, wholeLayer);
}

void write_bodies(const char * filename, const std::vector<Body> & bodies, MPI_Comm comm, bool overwrite){
    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);


    std::ofstream myfile;
    if(rank != 0){
        MPI_Status status;
        int x;
        MPI_Recv(&x, 1, MPI_INT, rank - 1, 0,
                MPI_COMM_WORLD, &status);
        myfile.open(filename, std::ios::app);
    }else{
        if(!file_exists(filename) or overwrite){
            myfile.open(filename);
        }else{
            myfile.open(filename, std::ios::app);
        }
    }

    for(const Body & b : bodies){
        myfile << b.idx << " " << b.pos[0] << " " << b.pos[1] << " " << b.m << " " << b.vel[0] << " " << b.vel[1] << std::endl;
    }

    if(rank != size - 1){
        myfile.close();
        int x = 1;
        MPI_Send(&x, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
    }else{
        myfile << std::endl;
        myfile.close();
    }
}

void write_to_file(const char * filename, const double x, bool overwrite){
    std::ofstream myfile;
    if(!file_exists(filename) or overwrite){
        myfile.open(filename);
    }else{
        myfile.open(filename, std::ios::app);
    }
    myfile << x << std::endl;;
    myfile.close();
}
