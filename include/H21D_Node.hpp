#ifndef __H21D_Node__
#define __H21D_Node__

#include "Eigen/Dense"
#include "Matrix.hpp"
#include "LowRank.hpp"

class H21D_Node 
{
friend class H21D;
// All methods are declared as private since all usage happens from
// the friend class H21D:
private:
    H21D_Node(int node_number, int level_number, int n_start, int n_size);
    
    // Storing the information passed to constructor as attribute:
    int node_number, level_number;
    int n_start, n_size;
    // Storing the start locations and sizes for the children of the node:
    int c_start[2], c_size[2];

    // Method to print the parameters of the node(mainly used to debug)
    void printNodeDetails();

    // Different Operators
    Mat self_interaction;        // Needed only at the leaf level.
    Mat neighbor_interaction[2]; // Neighbor interaction only needed at the leaf level.
	Mat M2L[16];                 // M2L of inner interactions. This is done on the box [-L,L]^2.
};

#endif /*__H21D_Node__*/
