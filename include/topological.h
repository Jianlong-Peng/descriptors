/*=============================================================================
#     FileName: topological.h
#         Desc: topological descriptors
#       Author: jlpeng
#        Email: jlpeng1201@gmail.com
#     HomePage: 
#      Version: 0.0.1
#      Created: 2012-12-28 15:38:59
#   LastChange: 2014-03-13 20:16:56
#      History:
=============================================================================*/

#ifndef  TOPOLOGICAL_H
#define  TOPOLOGICAL_H
#include <openbabel/mol.h>
#include <vector>
#include <string>
#include "config.h"

// Wiener index
double DESCRIPTOR_API wiener(OpenBabel::OBMol &mol);
// calculate Randic index (order 0~3)
void DESCRIPTOR_API randic(OpenBabel::OBMol &mol, double value[4]);
// calcualte Kier&Hall index (order 0~3)
void DESCRIPTOR_API kier_hall(OpenBabel::OBMol &mol, double value[4]);
// Kier shape index (order 1~3)
// value[0]: Kier flexibility index
void DESCRIPTOR_API kier_shape(OpenBabel::OBMol &mol, double value[4]);
// Balaban index
// if codessa is true, the following bond weight will be applied to calculate shortest distance
// single: 1.;  double: 1./2;  triple: 1./3;  aromatic: 2./3
double DESCRIPTOR_API balaban(OpenBabel::OBMol &mol, bool codessa=true);
// information content index and its derivatives
// averIC, IC, averSIC, SIC, averCIC, CIC, averBIC, BIC (order 0~2)
// Attention:
// 1. for polycyclicaromatichydrocarbons (PHAs), the information content order 1 and 2 can't be
//    calculated as same as CODESSA!!!!!!
void DESCRIPTOR_API infoContent(OpenBabel::OBMol &mol, double value[3][8]);
// to calculate electrotopological state indices for each heavy atom
// for hydrogen atoms, a value of -1. will be given.
// if failed, return empty vector
// reference: Hall, L. J.; Kier, L. B. J. Chem. Inf. Comput. Sci. 1995, 35: 1039-1045.
std::vector<double> calcEState(OpenBabel::OBMol &mol);
// SPAN: \frac{largest shortest distance from atom i to other atoms}{largest shortest distance between any two atoms}
std::vector<double> calcSPAN(OpenBabel::OBMol &mol);

// GETAWAY(GEmotery, Topology, and Atom-Weights AssemblY)
// ref: Consonni V. et al. J. Chem. Inf. Comput. Sci. 2002, 42: 682-692
void DESCRIPTOR_API getaway(OpenBabel::OBMol &mol, double value[]);
#endif   /* ----- #ifndef TOPOLOGICAL_H  ----- */

