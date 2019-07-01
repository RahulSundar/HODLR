#include "Matrix.hpp"
#include "LowRank.hpp"

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

    Mat getX() override
    {
        return this->x;
    }

    Mat getY() override
    {
        return this->y;
    }

    dtype getMatrixEntryXY(Mat x, Mat y, int i, int j) override
    {
        return 1 / abs((x(i)-y(j)));
    }

    // Destructor:
    ~Kernel() {};
};

int main(int argc, char* argv[]) 
{
    // Size of the Matrix:
    int N = 10;

    // Declaration of Matrix object that abstracts data in Matrix:
    Kernel* K  = new Kernel(N);
    LowRank* F = new LowRank(K, "interpolation1d");

    Mat B = K->getMatrix(0, 0, N, N);
    Mat L, R, error;

    F->getFactorization(L, R, 3);
    error = B - L * R.transpose();
    std::cout << "Accuracy of Factorization using Rook Pivoting:" << error.cwiseAbs().maxCoeff() << std::endl;
}
