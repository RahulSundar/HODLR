// Gives the interpolation operator / L2L for 1D:
Mat getL2L1D(Vec x, Vec x_nodes)
{
    // Size is the number of points in x
    int N = x.size();
    // Rank is the number of points in x_nodes
    int rank = x_nodes.size();
    Mat L2L(N, rank);

    // Tiling x       to bring it to shape (size, rank)
    // Tiling x_nodes to bring it to shape (rank, rank)
    Mat x_temp(N, rank);
    Mat x_nodes_temp(rank, rank);

    #pragma omp parallel for
    for(int i = 0; i < rank; i++)
    {
        x_temp.col(i)       = x;
        x_nodes_temp.col(i) = x_nodes;
    }

    // We will be using x_nodes_temp to get the denominator for 
    // the lagrange polynomials i.e. Π(x_i - x_j) where i != j
    // Similarly, we will be using x_temp to get the numerator for
    // the lagrange polynomials i.e. num(P_i) = Π(x - x_j) where i =! j 
    #pragma omp parallel for
    for(int i = 0; i < rank; i++)
    {
        x_temp.col(i)       = x.array() - x_nodes(i);
        x_nodes_temp.col(i) = x_nodes.array() - x_nodes(i);
        x_nodes_temp(i, i)  = 1;
    }
    
    // Performing Π(x_i - x_j) where i != j
    Mat prod = x_nodes_temp.colwise().prod();
    // temp is used to get num(P_i) = Π(x - x_j) where i =! j
    Mat temp;

    for(int i = 0; i < rank; i++)
    {
        temp        = x_temp;
        temp.col(i) = Vec::Ones(N);
        L2L.col(i)  = temp.rowwise().prod();
    }


    #pragma omp parallel for
    for(int i = 0; i < N; i++)
    {
        L2L.row(i) = L2L.row(i).array() / prod.array();
    }

    return L2L;
}

// Gives the interpolation operator / L2L for 2D:
Mat getL2L2D(Vec x, Vec y, Vec nodes)
{
    int N_nodes = nodes.size();

    // Getting individual L2L operators:
    Mat L2L_x = getL2L1D(x, nodes);
    Mat L2L_y = getL2L1D(y, nodes);
    
    Mat L2L(x.size(), N_nodes * N_nodes);
    #pragma omp parallel for
    for(int i = 0; i < N_nodes * N_nodes; i++)
    {
        L2L.col(i) =  (L2L_x.col(i % N_nodes)).array() 
                    * (L2L_y.col(i / N_nodes)).array();
    }

    return L2L;
}

Mat getM2L(Mat &nodes_1, Mat &nodes_2, Matrix* M)
{
    // Evaluating the Kernel function at the interpolation nodes:
    Mat M2L = M->getMatrix(0, 0, nodes_1.size(), nodes_2.size(), nodes_1, nodes_2);
    return M2L;
}
