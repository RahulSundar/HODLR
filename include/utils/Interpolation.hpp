#include "Matrix.hpp"

// Returns the roots on the N-th Legendre polynomial
// These values are in the standard interval of [-1, 1]
Vec getLegendreNodes(const int N)
{
    // Holds the coefficients of the polynomial:
    Vec poly(N+1);
    poly << Vec::Zero(N), 1;

    // Getting the scaled companion matrix:
    // This has been obtained following 
    // https://docs.scipy.org/doc/numpy-1.12.0/reference/generated/numpy.polynomial.legendre.legcompanion.html
    Mat companion = Mat::Zero(N, N);
    Vec scale     = ((2 * Vec::LinSpaced(N, 0.5, N-0.5)).cwiseSqrt()).array().inverse();

    // Assigning the values to the sub-diagonal and super-diagonal:
    companion.diagonal(-1) =   (Vec::LinSpaced(N-1, 1, N-1)).array()
                             * scale.segment(0, N-1).array()
                             * scale.segment(1, N-1).array();

    companion.diagonal(1) = companion.diagonal(-1);
    
    Vec nodes = companion.eigenvalues().real();
    std::sort(nodes.data(), nodes.data() + nodes.size());
    return nodes;
}

// Returns the roots on the N-th Chebyshev polynomial of the first kind
// These values are in the standard interval of [-1, 1]
Vec getChebyshevNodes(const int N)
{   
    Vec nodes(N);
    for(int i = 0; i < N; i++)
        nodes(i) = -cos(PI * (2 * i + 1) / (2 * N));
    
    return nodes;
}

// Returns the equispaced nodes in the standard interval of [-1, 1]
Vec getEquispacedNodes(const int N)
{
    return Vec::LinSpaced(N, -1 + 1/double(N), 1 - 1/double(N)); 
}

Vec getStandardNodes(int N_nodes, std::string nodes_type)
{
    // Obtaining the standard Chebyshev nodes of the first kind:
    if(nodes_type == "CHEBYSHEV")
    {
        return getChebyshevNodes(N_nodes);
    }

    // Obtaining the standard Legendre nodes:
    else if(nodes_type == "LEGENDRE")
    {
        return getLegendreNodes(N_nodes);
    }

    // Obtaining the standard equispaced nodes:
    else if(nodes_type == "EQUISPACED")
    {
        return getEquispacedNodes(N_nodes);
    }

    else
    {
        std::cout << "Invalid choice for interpolation type" << std::endl;
        std::cout << "Please use either CHEBYSHEV, LEGENDRE or EQUISPACED" << std::endl;
        exit(1);
    }
}

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

    // Allowing broadcasting:
    #pragma omp parallel for
    for(int i = 0; i < rank; i++)
    {
        L2L.row(i) = L2L.row(i) / x_nodes_temp;
    }

    return L2L;
}

Mat getM2LAF(Mat &nodes_1, array &nodes_2, Matrix* M)
{
    // Evaluating the Kernel function at the Chebyshev nodes:
    M2L = M->getMatrix(0, 0, nodes_1, nodes_2);
    M2L.eval();
}
