/*=============================================================================
#     FileName: dmath.h
#         Desc: 
#       Author: jlpeng
#        Email: jlpeng1201@gmail.com
#     HomePage: 
#      Created: 2014-03-13 19:45:44
#   LastChange: 2014-03-13 20:55:47
#      History:
=============================================================================*/

#ifndef  DMATH_H
#define  DMATH_H

#include "config.h"

/*******************************************************************************
Given a matrix a[1..m][1..n], this routine computes its singular value
decomposition, A = U.W.VT.  The matrix U replaces a on output.  The diagonal
matrix of singular values W is output as a vector w[1..n].  The matrix V (not
the transpose VT) is output as v[1..n][1..n].
*******************************************************************************/
void DESCRIPTOR_API svdcmp(double **a, int m, int n, double *w, double **v);

// leverage matrix
// H = M*(M'*M)^-1*M'
double **leverage(double **M, int row, int col);

#endif   /* ----- #ifndef DMATH_H  ----- */

