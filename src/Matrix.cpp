#include "Matrix.hpp"

void Matrix::checkInterpolation()
{
    // Just used to set the flag is_interp:
    dtype temp = getMatrixEntry(Mat::Random(1, 1), Mat::Random(1, 1), 0, 0);
}

Vec Matrix::getRow(const int j, const int n_col_start, const int n_cols, Mat x, Mat y) 
{
    this->checkInterpolation();
    Vec row(n_cols);

    if(this->is_interp == true)
    {
        #pragma omp parallel for
        for(int k = 0; k < n_cols; k++) 
        {   
            row(k) = this->getMatrixEntry(j,  k + n_col_start);
        }
    }

    else
    {
        #pragma omp parallel for
        for(int k = 0; k < n_cols; k++) 
        {   
            row(k) = this->getMatrixEntry(j,  k + n_col_start);
        }
    }

    return row;
}

Vec Matrix::getCol(const int k, const int n_row_start, const int n_rows, Mat x, Mat y) 
{
    this->checkInterpolation();
    Vec col(n_rows);

    if(this->is_interp == true)
    {
        #pragma omp parallel for
        for (int j=0; j<n_rows; ++j) 
        {
            col(j) = this->getMatrixEntry(j + n_row_start, k);
        }
    }

    else
    {
        #pragma omp parallel for
        for (int j=0; j<n_rows; ++j) 
        {
            col(j) = this->getMatrixEntry(j + n_row_start, k);
        }
    }
    

    return col;
}

Vec Matrix::getDiag1(const int n_row_start, const int n_col_start, 
                     const int n_rows, const int n_cols, Mat x, Mat y) 
{
    this->checkInterpolation();
    int N = std::max(n_rows, n_cols);
    Vec diag(N);

    int row_ind, col_ind;
    if(this->is_interp == true)
    {
        #pragma omp parallel for
        for (int j = 0; j < N; ++j) 
        {   
            if(n_cols > n_rows)
            {
                row_ind = this->mod(n_row_start - n_col_start + j, n_rows);
                col_ind = j;
            }

            else
            {
                row_ind = j;
                col_ind = this->mod(n_col_start - n_row_start + j, n_cols);
            }
            
            diag(j) = this->getMatrixEntry(row_ind, col_ind);
        }
    }

    else
    {
        #pragma omp parallel for
        for (int j = 0; j < N; ++j) 
        {   
            if(n_cols > n_rows)
            {
                row_ind = this->mod(n_row_start - n_col_start + j, n_rows);
                col_ind = j;
            }

            else
            {
                row_ind = j;
                col_ind = this->mod(n_col_start - n_row_start + j, n_cols);
            }
            
            diag(j) = this->getMatrixEntry(row_ind, col_ind);
        }
    }
    
    return diag;
}

Vec Matrix::getDiag2(const int n_row_start, const int n_col_start, 
                     const int n_rows, const int n_cols, Mat x, Mat y) 
{
    this->checkInterpolation();
    int N = std::max(n_rows, n_cols);
    Vec diag(N);
    
    int row_ind, col_ind;
    if(this->is_interp == true)
    {
        #pragma omp parallel for
        for (int j = 0; j < N; ++j) 
        {   
            if(n_cols > n_rows)
            {
                row_ind = this->mod(n_row_start + n_col_start - j, n_rows);
                col_ind = j;
            }

            else
            {
                row_ind = j;
                col_ind = this->mod(n_col_start + n_row_start - j, n_cols);
            }

            diag(j) = this->getMatrixEntry(row_ind, col_ind);
        }
    }

    else
    {
        #pragma omp parallel for
        for (int j = 0; j < N; ++j) 
        {   
            if(n_cols > n_rows)
            {
                row_ind = this->mod(n_row_start + n_col_start - j, n_rows);
                col_ind = j;
            }

            else
            {
                row_ind = j;
                col_ind = this->mod(n_col_start + n_row_start - j, n_cols);
            }

            diag(j) = this->getMatrixEntry(row_ind, col_ind);
        }
    }

    return diag;
}

Mat Matrix::getMatrix(const int n_row_start, const int n_col_start, 
                      const int n_rows, const int n_cols, Mat x, Mat y) 
{
    Mat mat(n_rows, n_cols);
    
    if(this->is_interp == true)
    {
        #pragma omp parallel for
        for (int j=0; j < n_rows; ++j) 
        {
            #pragma omp parallel for
            for (int k=0; k < n_cols; ++k) 
            {
                mat(j,k) = this->getMatrixEntry(j + n_row_start, k + n_col_start);
            }
        }
    }

    else
    {
        #pragma omp parallel for
        for (int j=0; j < n_rows; ++j) 
        {
            #pragma omp parallel for
            for (int k=0; k < n_cols; ++k) 
            {
                mat(j,k) = this->getMatrixEntry(j + n_row_start, k + n_col_start);
            }
        }
    }

    return mat;
}
