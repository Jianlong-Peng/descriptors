/*=============================================================================
#     FileName: geometrical.h
#         Desc: geometrical features
#       Author: jlpeng
#        Email: jlpeng1201@gmail.com
#     HomePage: 
#      Version: 0.0.1
#      Created: 2012-12-28 20:38:01
#   LastChange: 2013-12-27 16:52:23
#      History:
=============================================================================*/

#ifndef  GEOMETRICAL_H
#define  GEOMETRICAL_H
#include <vector>
#include <string>
#include <utility>
#include <openbabel/babelconfig.h>
#include <openbabel/mol.h>
#include "config.h"

// principle moment of inertia
std::vector<std::pair<int, double> > DESCRIPTOR_API calcMOI(OpenBabel::OBMol &mol);
// - int axes[3]: principle inertia axes (MOI: axes[0] < axes[1] < axes[2])
//   0: X-axis; 1: Y-axis; 2: Z-axis
//   if axes[0]==-1, then mol will be aligned to principle inertia,
//   otherwise, using XYZ-axes indicated by axes
// - if you called 'calcMOI' first, then you can do it as follows:
//     vector<pair<int,double> > moi = calcMOI(mol);
//     int axes[3] = {moi[0].first, moi[1].first, moi[2].first};
//     shadow(mol,axes,density,result);
// - return: [XYshadow,ZXShadow,YZShadow,XYShadow/XYR,ZXShadow/ZXR,YZShadow/YZR]
// - ref: Rohrbaugh R.H. etc. Analytica Chimica Acta. 1987, 199: 99-109
void DESCRIPTOR_API shadow(OpenBabel::OBMol &mol, int axes[3], int density, double result[6]);
// calculate molecular volume  using van der Waals radius
// 'axes': same as that of 'shadow'
// density: (recommend)>= 10
// return: [MV, MV/XYZBox]
void DESCRIPTOR_API volume(OpenBabel::OBMol &mol, int axes[3], int density, double result[2]);
// Molecular surface area
double DESCRIPTOR_API surfaceArea(OpenBabel::OBMol &mol);  // NOT implemented!!
// solvent accessible surface area
// ref: Shrake A. etc. J. Mol. Biol. 1973, 79(2): 351-371
std::vector<double> DESCRIPTOR_API calcSASA(OpenBabel::OBMol &mol);
// to calculate charged partial (solvent-accessible) surface area
// mol   : object of class OBMol
// pc    : partial charge, size == mol.NumAtoms()
// result: where results will be saved.
bool DESCRIPTOR_API calcCPSA(OpenBabel::OBMol &mol, std::vector<double> &pc, double result[38]);
//extern std::vector<std::string> DESCRIPTOR_API cpsa_names;
// topologiial polar surface area
// ref: Ertl P. etc., J. Med. Chem. 2000, 43(20): 3714-3717
double DESCRIPTOR_API calcTPSA(OpenBabel::OBMol &mol);
// Gravitation index   <-- need 3D coordinates !!!!!
// value[0]: summation over all pairs of atoms
// value[1]: summation over all bonded pairs of atoms
void DESCRIPTOR_API gravitation(OpenBabel::OBMol &mol, double value[]);
#endif   /* ----- #ifndef GEOMETRICAL_H  ----- */

