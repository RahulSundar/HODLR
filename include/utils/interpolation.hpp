// Gives the interpolation operator / L2L for 1D:
Mat getL2L1D(Vec x, Vec x_nodes)
{
    // Size is the number of points in x
    int N = x.size()
    // Rank is the number of points in x_nodes
    int rank = x_nodes.size()
    Mat L2L(N, rank);

    // Tiling x       to bring it to shape (size, rank)
    // Tiling x_nodes to bring it to shape (rank, rank)
    Mat x_temp(N, rank);
    Mat x_nodes_temp(rank, rank);

    #pragma omp parallel for
    for(int i = 0; i < rank; i++)
    {
        x_temp.col(i) = x;
        x_nodes_temp.col(i) = x_nodes;
    }

    // We will be using x_nodes to get the denominator for 
    // the lagrange polynomials i.e. Π(x_i - x_j) where i != j
    // Similarly, we will be using x to get the numerator for
    // the lagrange polynomials i.e. num(P_i) = Π(x - x_j) where i =! j 
    #pragma omp parallel for
    for(int i = 0; i < rank; i++)
    {
        x_temp.col(i)       = x - x_nodes(i);
        x_nodes_temp.col(i) = x_nodes - x_nodes(i);
        x_nodes_temp(i, i)  = 1
    }
    
    // Performing Π(x_i - x_j) where i != j
    x_nodes_temp = x_nodes_temp.rowwise().prod();
    // temp is used to get num(P_i) = Π(x - x_j) where i =! j
    Mat temp;

    for(int i = 0; i < rank; i++)
    {
        temp        = x_temp;
        temp.col(i) = 1;
        L2L.col(i)  = temp.rowwise().prod();
    }

    #pragma omp parallel for
    for(int i = 0; i < rank; i++)
    {
        L2L.row(i) = L2L.row(i) / x_nodes_temp;
    }

    return L2L;
}

Mat getM2L(Mat &nodes_1, array &nodes_2, Matrix* M)
{
    // Evaluating the Kernel function at the interpolation nodes:
    Mat M2L = M->getMatrix(0, 0, nodes_1.size(), nodes_2.size(), nodes_1, nodes_2);
    return M2L;
}

void getInterpolation()
{
    if(M.getDimensionality() == 1)
    {
        double c_targets, r_targets, c_sources, r_sources;
        // Determining the center and radius of targets:
        determineCenterAndRadius(target_coords, c_targets, r_targets);
        // Determining the center and radius of sources:
        determineCenterAndRadius(source_coords, c_sources, r_sources);

        // Obtain the scaled Chebyshev nodes for the targets:
        array nodes_targets, nodes_sources;

        scalePoints(0, 1, standard_nodes, c_targets, r_targets, nodes_targets);
        scalePoints(0, 1, standard_nodes, c_sources, r_sources, nodes_sources);

        // Standard Locations of the coordinates:
        array standard_targets, standard_sources;
        scalePoints(c_targets, r_targets, target_coords, 0, 1, standard_targets);
        scalePoints(c_sources, r_sources, source_coords, 0, 1, standard_sources);
    
        // Initializing U, S, V:
        U = af::constant(0, M.getNumRows(), n_nodes, f64);
        S = array(n_nodes, n_nodes, f64);
        V = af::constant(0, M.getNumCols(), n_nodes, f64);

        getL2L1DAF(standard_targets, standard_nodes, U);
        getM2LAF(nodes_targets, nodes_sources, M, S);
        getL2L1D(standard_sources, standard_nodes, V);
    }
}
