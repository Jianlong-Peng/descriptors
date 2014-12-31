/*=============================================================================
#         FileName: mymatrix.h
#             Desc: 
#           Author: jlpeng
#            Email: jlpeng1201@gmail.com
#         HomePage: 
#          Version: 0.0.1
#          Created: 2013-01-06 12:35:58
#   LastChange: 2013-02-20 19:29:46
#          History:
=============================================================================*/

#ifndef  MATRIX_H
#define  MATRIX_H
#include "config.h"

namespace Descriptors {
    class DESCRIPTOR_API Matrix {
    public:
        Matrix(): m(0), _row(0), _col(0) {}
        Matrix(int row, int col, double v=0.);
        Matrix(const Matrix &rm);
        ~Matrix();

        int nrows() const {return _row;}
        int ncols() const {return _col;}

        double* operator[](int i) {return m[i];}
        const double* operator[](int i) const {return m[i];}
        double& operator()(int i, int j) {return m[i][j];}
        const double& operator()(int i, int j) const {return m[i][j];}
        Matrix& operator=(const Matrix &rm);
        Matrix operator+(const Matrix &rm);
        Matrix operator-(const Matrix &rm);
        Matrix operator*(const Matrix &rm);
    private:
        double **m;
        int _row; // number of rows
        int _col; // number of columns
    };

    // to solve linear equations:
    //         a11*x1 + a12*x2 + ... + a1n*xn = b1
    //         ...
    //         an1*x1 + an2*x2 + ... + ann*xn = bn
    // parameters:
    //   m: matrix of n*(n+1)
    //   answer: where solution will be saved. size = n
    // return variable:
    //   true if success, otherwise false
    bool DESCRIPTOR_API linearEquations(const Matrix &m, double *answer);
}

#endif   /* ----- #ifndef MATRIX_H  ----- */

