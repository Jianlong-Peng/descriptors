/*=============================================================================
#     FileName: physicochemical.h
#         Desc: 
#       Author: jlpeng
#        Email: jlpeng1201@gmail.com
#     HomePage: 
#      Created: 2013-11-30 13:58:50
#   LastChange: 2013-12-27 14:23:49
#      History:
=============================================================================*/

#ifndef  PHYSICOCHEMICAL_H
#define  PHYSICOCHEMICAL_H

#include <vector>
#include <string>
#include "config.h"
#include <openbabel/mol.h>


// to calculate logP using OpenBabel
// if failed to do the calculation, return `1E+38`
// reference:
//   Wildman, S.A. and Crippen, G.M., J. Chem. Inf. Comput. Sci., 1999, 39, 868-873.
double DESCRIPTOR_API calcLogP(OpenBabel::OBMol &mol);

#endif   /* ----- #ifndef PHYSICOCHEMICAL_H  ----- */

