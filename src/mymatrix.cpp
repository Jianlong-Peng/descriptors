/*=============================================================================
#         FileName: mymatrix.cpp
#             Desc: 
#           Author: jlpeng
#            Email: jlpeng1201@gmail.com
#         HomePage: 
#          Version: 0.0.1
#          Created: 2013-01-06 12:39:03
#   LastChange: 2013-01-24 11:56:13
#          History:
=============================================================================*/
#include <iostream>
#include <cstring>
#include <cstdlib>
#include <cassert>
#include "mymatrix.h"

using std::cerr;
using std::endl;

namespace Descriptors
{
    Matrix::Matrix(int row, int col, double v)
    {
        _row = row;
        _col = col;
        m = new double* [row];
        if(!m) {
                cerr << "matrix.cpp: out of memory!" << endl;
                exit(EXIT_FAILURE);
        }
        for(int i=0; i<row; ++i) {
                m[i] = new double [col];
                if(!m[i]) {
                        cerr << "matrix.cpp: out of memory!" << endl;
                        exit(EXIT_FAILURE);
                }
                for(int j=0; j<col; ++j)
                        m[i][j] = v;
        }
    }

    Matrix::Matrix(const Matrix &rm)
    {
        _row = rm.nrows();
        _col = rm.ncols();
        m = new double* [_row];
        if(!m) {
                cerr << "matrix.cpp: out of memory!" << endl;
                exit(EXIT_FAILURE);
        }
        for(int i=0; i<_row; ++i) {
                m[i] = new double [_col];
                if(!m[i]) {
                        cerr << "matrix.cpp: out of memory!" << endl;
                        exit(EXIT_FAILURE);
                }
                for(int j=0; j<_col; ++j)
                        m[i][j] = rm(i,j);
        }
    }

    Matrix::~Matrix()
    {
        for(int i=0; i<_row; ++i)
                delete[] m[i];
        delete[] m;
    }

    Matrix& Matrix::operator=(const Matrix &rm)
    {
        if(this == &rm)
                return *this;
        if(_row==0 && _col==0) {
                _row = rm.nrows();
                _col = rm.ncols();
                m = new double* [_row];
                if(!m) {
                        cerr << "matrix.cpp: out of memory!" << endl;
                        exit(EXIT_FAILURE);
                }
                for(int i=0; i<_row; ++i) {
                        m[i] = new double [_col];
                        if(!m[i]) {
                                cerr << "matrix.cpp: out of memory!" << endl;
                                exit(EXIT_FAILURE);
                        }
                }
        }
        assert(_row==rm.nrows() && _col==rm.ncols());
        for(int i=0; i<_row; ++i)
                for(int j=0; j<_col; ++j)
                        (*this)(i,j) = rm(i,j);
        return *this;
    }

    Matrix Matrix::operator+(const Matrix &rm)
    {
        assert(_row==rm.nrows() && _col==rm.ncols());
        Matrix result(_row,_col);
        for(int i=0; i<_row; ++i)
                for(int j=0; j<_col; ++j)
                        result(i,j) = (*this)(i,j) + rm(i,j);
        return result;
    }

    Matrix Matrix::operator-(const Matrix &rm)
    {
        assert(_row==rm.nrows() && _col==rm.ncols());
        Matrix result(_row,_col);
        for(int i=0; i<_row; ++i)
                for(int j=0; j<_col; ++j)
                        result(i,j) = (*this)(i,j) - rm(i,j);
        return result;
    }

    Matrix Matrix::operator*(const Matrix &rm)
    {
        assert(_col == rm.nrows());
        Matrix result(_row,rm.ncols());
        for(int i=0; i<_row; ++i) {
                double value(0.);
                for(int j=0; j<rm.ncols(); ++j) {
                        for(int k=0; k<_col; ++k)
                                value += ((*this)(i,k) * rm(k,j));
                        result(i,j) = value;
                }
        }
        return result;
    }


    // solve a set of linear equations
    // solution.nrows()+1 == solution.ncols()
    // length(answer) == solution.nrows()
    bool linearEquations(const Matrix &solution, double *answer)
    {
        assert(solution.nrows()+1 == solution.ncols());
        
        Matrix data(solution);
        int j, m, n;
        double tmp, *d;
        int nr_variable = data.nrows();
        
        //a temp array
        d = new double [nr_variable+1];
        if(d == NULL) {
                cerr << "linearAlgebra.c::linearEquations: not enough memory!" << endl;
                return false;
        }
        for (m = 0; m < nr_variable - 1; m ++)
        {
                //if data[m][m]==0, then exchange with line who got nonzero data[i][i]
                for (n = m + 1; n < nr_variable && data(m,m) == 0.0; n ++)
                {
                        if ( data(n,m) != 0.0)
                        {
                                memcpy(d, data[m], sizeof(double)*(nr_variable + 1));
                                memcpy(data[m], data[n], sizeof(double)*(nr_variable + 1));
                                memcpy(data[n], d, sizeof(double)*(nr_variable + 1));
                        }
                }
                //
                if (data(m,m) == 0.0)
                {
                        delete[] d;
                        return false;
                }
                
                //elimination
                for (n = m + 1; n < nr_variable; n ++)
                {
                        tmp = data(n,m) / data(m,m);
                        for (j = m; j <= nr_variable; j ++)
                                data(n,j) -= tmp * data(m,j);
                }
        }
        for (j = 0; j < nr_variable; j ++)
                d[j] = 0.0;
        //
        answer[nr_variable - 1] = data(nr_variable-1,nr_variable) / data(nr_variable-1,nr_variable-1);
        //
        for (m = nr_variable - 2; m >= 0; m --)
        {
                for (j = nr_variable - 1; j > m; j --)
                        d[m] += answer[j] * data(m,j);
                answer[m] = (data(m,nr_variable) - d[m]) / data(m,m);
        }
        
        delete[] d;
        return true;
    }
}


