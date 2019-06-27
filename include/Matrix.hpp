#ifndef __Matrix__
#define __Matrix__

#ifdef MKL_ENABLED
    #define EIGEN_USE_MKL_ALL
#endif

#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
#include <set>
#include <vector>
// Used to dump data:
#include <fstream>
#include <complex>

#ifdef USE_FLOAT
    using dtype=float;
    using dtype_base=float;
    using Mat=Eigen::MatrixXf;
    using Vec=Eigen::VectorXf;
#endif

#ifdef USE_DOUBLE
    using dtype=double;
    using dtype_base=double;
    using Mat=Eigen::MatrixXd;
    using Vec=Eigen::VectorXd;
#endif

#ifdef USE_COMPLEX32
    using dtype=std::complex<float>;
    using dtype_base=float;
    using Mat=Eigen::MatrixXcf;
    using Vec=Eigen::VectorXcf;
    const std::complex<float> I(0.0, 1.0);
#endif

#ifdef USE_COMPLEX64
    using dtype=std::complex<double>;
    using dtype_base=double;
    using Mat=Eigen::MatrixXcd;
    using Vec=Eigen::VectorXcd;
    const std::complex<double> I(0.0, 1.0);
#endif

#define PI 3.141592653589793238462643383279502884

class Matrix 
{
public:

    // Size of the matrix:
    int N;
    // These would be declared when performing lowrank using interpolation:
    Mat x, y;

    // Modulo operator:
    // This is separately defined to make sure 
    // that positive values are always returned
    int mod(int a, int b)
    {
        return ((a % b + b) % b);
    }

    // Constructor:
    explicit Matrix(int N)
    {
        this->N         = N;
    }

    virtual Vec getX()
    {
        return x;
    }

    virtual Vec getY()
    {
        return y;
    }

    // Returns individual entries of the matrix:
    virtual dtype getMatrixEntry(int j, int k) 
    {
        return -1.0;
    }

    // Returns individual entries of the matrix, when entries are a function of x, y:
    virtual dtype getMatrixEntryXY(Mat x, Mat y, int i, int j) 
    {
        return getMatrixEntry(i, j);
    }

    Vec getRow(int j, int n_col_start, int n_cols, Mat x = Mat::Zero(1, 1), Mat y = Mat::Zero(1, 1));
    Vec getCol(int k, int n_row_start, int n_rows, Mat x = Mat::Zero(1, 1), Mat y = Mat::Zero(1, 1));
    Vec getDiag1(int j, int k, int n_rows, int n_cols, Mat x = Mat::Zero(1, 1), Mat y = Mat::Zero(1, 1));
    Vec getDiag2(int j, int k, int n_rows, int n_cols, Mat x = Mat::Zero(1, 1), Mat y = Mat::Zero(1, 1));
    Mat getMatrix(int j, int k, int n_rows, int n_cols, Mat x = Mat::Zero(1, 1), Mat y = Mat::Zero(1, 1));

    // Destructor:
    ~Matrix() {};
};

#endif /*__Matrix__*/
