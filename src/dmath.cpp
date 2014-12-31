/*=============================================================================
#     FileName: dmath.cpp
#         Desc: 
#       Author: jlpeng
#        Email: jlpeng1201@gmail.com
#     HomePage: 
#      Created: 2014-03-13 19:46:03
#   LastChange: 2014-03-20 17:06:52
#      History:
=============================================================================*/
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "dmath.h"


/*******************************************************************************
Singular value decomposition program, svdcmp, from "Numerical Recipes in C"
(Cambridge Univ. Press) by W.H. Press, S.A. Teukolsky, W.T. Vetterling,
and B.P. Flannery
*******************************************************************************/

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define MAX(a,b) ((a) > (b) ? a : b)
#define MIN(a,b) ((a) < (b) ? a: b)

static double pythag(double a, double b)
/* compute (a^2 + b^2)^1/2 without destructive underflow or overflow */
{
    double absa,absb;
    absa=fabs(a);
    absb=fabs(b);
    if (absa > absb) return absa*sqrt(1.0+(absb/absa)*(absb/absa));
    else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+(absa/absb)*(absa/absb)));
}


void svdcmp(double **a, int m, int n, double *w, double **v)
/*******************************************************************************
Given a matrix a[1..m][1..n], this routine computes its singular value
decomposition, A = U.W.VT.  The matrix U replaces a on output.  The diagonal
matrix of singular values W is output as a vector w[1..n].  The matrix V (not
the transpose VT) is output as v[1..n][1..n].
*******************************************************************************/
{
    int flag,i,its,j,jj,k,l,nm;
    double anorm,c,f,g,h,s,scale,x,y,z,*rv1;

    rv1 = (double*)malloc(sizeof(double)*(n+1));
    g = scale = anorm = 0.0; /* Householder reduction to bidiagonal form */
    for(i=1; i<=n; i++) {
        l = i + 1;
        rv1[i] = scale * g;
        g = s = scale = 0.0;
        if(i <= m) {
            for(k=i; k<=m; k++)
                scale += fabs(a[k][i]);
            if(scale) {
                for(k=i; k<=m; k++) {
                    a[k][i] /= scale;
                    s += a[k][i]*a[k][i];
                }
                f = a[i][i];
                g = -SIGN(sqrt(s),f);
                h = f * g - s;
                a[i][i] = f - g;
                for(j=l; j<=n; j++){
                    for(s=0.0,k=i; k<=m; k++)
                        s += a[k][i]*a[k][j];
                    f = s / h;
                    for(k=i;k<=m;k++)
                        a[k][j] += f*a[k][i];
                }
                for(k=i;k<=m;k++)
                    a[k][i] *= scale;
            }
        }
        w[i] = scale * g;
        g = s = scale = 0.0;
        if(i<=m && i!=n) {
            for(k=l; k<=n; k++)
                scale += fabs(a[i][k]);
            if(scale) {
                for (k=l; k<=n; k++) {
                    a[i][k] /= scale;
                    s += a[i][k]*a[i][k];
                }
                f = a[i][l];
                g = -SIGN(sqrt(s),f);
                h = f * g - s;
                a[i][l] = f - g;
                for(k=l;k<=n;k++)
                    rv1[k] = a[i][k] / h;
                for(j=l;j<=m;j++) {
                    for(s=0.0,k=l; k<=n; k++)
                        s += a[j][k]*a[i][k];
                    for(k=l; k<=n; k++)
                        a[j][k] += s*rv1[k];
                }
                for(k=l; k<=n; k++)
                    a[i][k] *= scale;
            }
        }
        anorm = MAX(anorm,(fabs(w[i])+fabs(rv1[i])));
    }
    for (i=n; i>=1; i--) { /* Accumulation of right-hand transformations. */
        if(i < n) {
            if(g) {
                for(j=l; j<=n; j++) /* Double division to avoid possible underflow. */
                    v[j][i] = (a[i][j] / a[i][l]) / g;
                for(j=l; j<=n; j++) {
                    for(s=0.0,k=l; k<=n; k++)
                        s += a[i][k]*v[k][j];
                    for(k=l; k<=n; k++)
                        v[k][j] += s*v[k][i];
                }
            }
            for(j=l; j<=n; j++)
                v[i][j] = v[j][i] = 0.0;
        }
        v[i][i] = 1.0;
        g = rv1[i];
        l = i;
    }
    for(i=MIN(m,n); i>=1; i--) { /* Accumulation of left-hand transformations. */
        l = i+1;
        g = w[i];
        for(j=l; j<=n; j++)
            a[i][j] = 0.0;
        if(g) {
            g = 1.0/g;
            for(j=l; j<=n; j++) {
                for (s=0.0,k=l; k<=m; k++)
                    s += a[k][i]*a[k][j];
                f = (s/a[i][i])*g;
                for (k=i; k<=m; k++)
                    a[k][j] += f*a[k][i];
            }
            for (j=i; j<=m; j++)
                a[j][i] *= g;
        } else {
            for (j=i;j<=m;j++)
                a[j][i]=0.0;
        }
        ++a[i][i];
    }
    for (k=n; k>=1; k--) { /* Diagonalization of the bidiagonal form. */
        for (its=1; its<=30; its++) {
            flag = 1;
            for (l=k; l>=1; l--) { /* Test for splitting. */
                nm = l-1; /* Note that rv1[1] is always zero. */
                if ((double)(fabs(rv1[l])+anorm) == anorm) {
                    flag=0;
                    break;
                }
                if ((double)(fabs(w[nm])+anorm) == anorm)
                    break;
            }
            if (flag) {
                c=0.0; /* Cancellation of rv1[l], if l > 1. */
                s=1.0;
                for (i=l; i<=k; i++) {
                    f = s * rv1[i];
                    rv1[i] = c*rv1[i];
                    if ((double)(fabs(f)+anorm) == anorm)
                        break;
                    g = w[i];
                    h = pythag(f,g);
                    w[i] = h;
                    h = 1.0/h;
                    c = g*h;
                    s = -f*h;
                    for (j=1; j<=m; j++) {
                        y = a[j][nm];
                        z = a[j][i];
                        a[j][nm] = y*c+z*s;
                        a[j][i] = z*c-y*s;
                    }
                }
            }
            z = w[k];
            if (l == k) { /* Convergence. */
                if (z < 0.0) { /* Singular value is made nonnegative. */
                    w[k] = -z;
                    for (j=1; j<=n; j++)
                        v[j][k] = -v[j][k];
                }
                break;
            }
            if (its == 30) printf("no convergence in 30 svdcmp iterations");
            x = w[l]; /* Shift from bottom 2-by-2 minor. */
            nm = k-1;
            y = w[nm];
            g = rv1[nm];
            h = rv1[k];
            f = ((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
            g = pythag(f,1.0);
            f = ((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
            c = s = 1.0; /* Next QR transformation: */
            for (j=l;j<=nm;j++) {
                i = j+1;
                g = rv1[i];
                y = w[i];
                h = s*g;
                g = c*g;
                z = pythag(f,h);
                rv1[j] = z;
                c = f/z;
                s = h/z;
                f = x*c+g*s;
                g = g*c-x*s;
                h = y*s;
                y *= c;
                for (jj=1; jj<=n; jj++) {
                    x = v[jj][j];
                    z = v[jj][i];
                    v[jj][j] = x*c+z*s;
                    v[jj][i] = z*c-x*s;
                }
                z = pythag(f,h);
                w[j] = z; /* Rotation can be arbitrary if z = 0. */
                if (z) {
                    z = 1.0/z;
                    c = f*z;
                    s = h*z;
                }
                f = c*g+s*y;
                x = c*y-s*g;
                for (jj=1;jj<=m;jj++) {
                    y = a[jj][j];
                    z = a[jj][i];
                    a[jj][j] = y*c+z*s;
                    a[jj][i] = z*c-y*s;
                }
            }
            rv1[l] = 0.0;
            rv1[k] = f;
            w[k] = x;
        }
    }
    free(rv1);
}


// leverage matrix
// H = M*(M'*M)^-1*M'
double **leverage(double **M, int row, int col)
{
    int i,j,k,l;
    double **U = (double**)malloc(sizeof(double*)*col);
    double *W = (double*)malloc(sizeof(double)*col);
    double **V = (double**)malloc(sizeof(double*)*col);
    double **MI = (double**)malloc(sizeof(double*)*col); // (M'*M)^-1
    double **H = (double**)malloc(sizeof(double*)*row);
    for(i=0; i<col; ++i) {
        U[i] = (double*)malloc(sizeof(double)*col);
        V[i] = (double*)malloc(sizeof(double)*col);
        MI[i] = (double*)malloc(sizeof(double)*col);
    }
    for(i=0; i<row; ++i)
        H[i] = (double*)malloc(sizeof(double)*row);
    // 1. M'*M
    for(i=0; i<col; ++i) {
        for(j=0; j<col; ++j) {
            U[i][j] = 0.;
            for(k=0; k<row; ++k)
                U[i][j] += (M[k][i] * M[k][j]);
        }
    }
    // 2. MI = (M'*M)^-1
    // if A = U*W*V', then A+ = V*W+*U'
    svdcmp(U, col, col, W, V);
    for(i=0; i<col; ++i) {
        for(j=0; j<col; ++j) {
            MI[i][j] = 0.;
            for(k=0; k<col; ++k) {
                if(W[j] != 0.)
                    MI[i][j] += (V[i][k] / W[k] * U[j][k]);
            }
        }
    }
    // 3. H = M*MI*M'
    for(i=0; i<row; ++i) {
        for(j=0; j<row; ++j) {
            H[i][j] = 0.;
            for(k=0; k<col; ++k) {
                double tmp = 0.;
                for(l=0; l<col; ++l)
                    tmp += (M[i][l] * MI[l][k]);
                H[i][j] += (tmp * M[j][k]);
            }
        }
    }
    // 4. free
    for(i=0; i<col; ++i) {
        free(U[i]);
        free(V[i]);
        free(MI[i]);
    }
    free(U);
    free(V);
    free(W);
    free(MI);

    return H;
}


