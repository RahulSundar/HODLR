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
