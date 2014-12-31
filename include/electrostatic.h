/*=============================================================================
#     FileName: electrostatic.h
#         Desc: 
#       Author: jlpeng
#        Email: jlpeng1201@gmail.com
#     HomePage: 
#      Version: 0.0.1
#      Created: 2013-01-08 11:04:34
#   LastChange: 2014-01-08 12:38:27
#      History:
=============================================================================*/

#ifndef  ELECTROSTATIC_H
#define  ELECTROSTATIC_H
#include <vector>
#include <string>
#include <openbabel/mol.h>
#include "config.h"

// to calculate atomic partial charge (Zefirov)
// parameters:
//   mol: molecule
//   scale_method: electronegativity
//     1 - Pauling scale
//     2 - Sanderson scale
// not same as those calculated by CODESSA !!!!!!!!!!
std::vector<double> DESCRIPTOR_API calcZefirovPC(OpenBabel::OBMol &mol, int scale_method=1);

// to calculate max/min partial charge for a * atom [Zefirov's PC]
// mol   : object of class OBMol
// result: final result
//         [0~9][0]: max partial charge for a specific atom
//         [0~9][1]: min partial charge for a specific atom
//         H,C,N,O,F,P,S,Cl,Br,I
//         [10][0]: max partial charge
//         [10][1]: min partial charge
//         [11][0]: polarity parameter (Qmax-Qmin)
//         [11][1]: polarity parameter / square distance
bool DESCRIPTOR_API calcMaxMinZefirovPC(OpenBabel::OBMol &mol, std::vector<double> &pc, 
        double result[12][2]);

// to calculate Topolographic electronic index [Zefirov's PC]
// result[0]: for all pairs of atoms
// result[1]: for all bonded pairs of atoms
bool DESCRIPTOR_API calcTopolElectronicIndex(OpenBabel::OBMol &mol, std::vector<double> &pc,
        double result[2]);

// - to calculate charged partial surface area [Zefirov's PC]
// for CPSA descriptors using Zefirov's PC, you can call 'calcCPSA' in 'geometrical.h'
// bool calcCPSA(OBMol &mol, vector<double> &pc, double result[38])
#endif   /* ----- #ifndef ELECTROSTATIC_H  ----- */

