#include "Matrix.hpp"
// #include "LowRank.hpp"
// #include "Interpolation.hpp"

// Derived class of Matrix which is ultimately
// passed to the HODLR_Tree class:
class Kernel : public Matrix 
{
public:
    Mat x, y;

    // Constructor:
    Kernel(int N) : Matrix(N) 
    {
        x = 2 * Mat::Ones(N, 1) + Mat::Random(N, 1);
        y = 7 * Mat::Ones(N, 1) + Mat::Random(N, 1);
    };

    dtype getMatrixEntry(Mat x, Mat y, int i, int j) 
    {
        return 1 / abs((x(i)-y(j)));
    }

    // Destructor:
    ~Kernel() {};
};

int main(int argc, char* argv[]) 
{
    // Declaration of Matrix object that abstracts data in Matrix:
    Kernel* K  = new Kernel(5);
    // LowRank* F = new LowRank(K, "rookPivoting");

    Mat B = K->getMatrix(0, 0, 5, 5);
    // Mat L, R, error;
    std::cout << B << std::endl;

    // F->getFactorization(L, R, 1e-12);
    // error = B - L * R.transpose();
    // std::cout << "Accuracy of Factorization using Rook Pivoting:" << error.cwiseAbs().maxCoeff() << std::endl;
}
