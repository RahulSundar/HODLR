#include "Matrix.hpp"

extern "C" 
{
    double get_matrix_entry(int* i, int* j);
}

class Kernel : public Matrix
{
private:
    Mat x;

public:
    // Constructor:
    Kernel(int N) : Matrix(N) 
    {
        x = (Mat::Random(N, 1)).real();
    };
    
    // In this example, we are illustrating usage using
    // the gaussian kernel:
    dtype getMatrixEntry(int i, int j) 
    {
        return get_matrix_entry(&i, &j);
    }

    // Destructor:
    ~Kernel(){}
};
