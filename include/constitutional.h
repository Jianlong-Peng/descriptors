/*=============================================================================
#     FileName: constitutional.h
#         Desc: descriptors that reflect only molecular composition of the compound
#       Author: jlpeng
#        Email: jlpeng1201@gmail.com
#     HomePage: 
#      Version: 0.0.1
#      Created: 2012-12-28 14:31:25
#   LastChange: 2014-01-08 13:35:06
#      History:
=============================================================================*/

#ifndef  CONSTITUTIONAL_H
#define  CONSTITUTIONAL_H
#include <map>
#include <vector>
#include <string>
#include <openbabel/mol.h>
#include "config.h"

/* 23
 * natoms  : Number of atoms
 * nCatoms : Number of C  atoms
 * rnCatoms: Relative number of C  atoms
 * nHatoms : Number of H  atoms
 * rnHatoms: Relative number of H  atoms
 * nOatoms : Number of O  atoms
 * rnOatoms: Relative number of O  atoms
 * nNatoms : Number of N  atoms
 * rnNatoms: Relative number of N  atoms
 * nSatoms : Number of S  atoms
 * rnSatoms: Relative number of S  atoms
 * nFatoms : Number of F  atoms
 * rnFatoms: Relative number of F  atoms
 * nClatoms: Number of Cl atoms
 *rnClatoms: Relative number of Cl atoms
 * nBratoms: Number of Br atoms
 *rnBratoms: Relative number of Br atoms
 * nIatoms : Number of I  atoms
 * rnIatoms: Relative number of I  atoms
 * nPatoms : Number of P  atoms
 * rnPatoms: Relative number of P  atoms
 * MW      : molecular wieght
 * rMW     : average atomic weight
 */
std::vector<double> DESCRIPTOR_API atoms(OpenBabel::OBMol &mol);
/* 9
 * nbonds  : Number of bonds
 * nsbonds : Number of single bonds
 * rnsbonds: Relative number of single bonds
 * ndbonds : Number of double bonds
 * rndbonds: Relative number of double bonds
 * ntbonds : Number of triple bonds
 * rntbonds: Relative number of triple bonds
 * narbonds: Number of aromatic bonds
 *rnarbonds: Relative number of aromatic bonds
*/
std::vector<double> DESCRIPTOR_API bonds(OpenBabel::OBMol &mol);
/* 4
 * nrings : Number of rings
 * rnrings: Relative number of rings
 * nbrings: Number of benzene rings
 *rnbrings: Relative number of benzene rings
 */
std::vector<double> DESCRIPTOR_API rings(OpenBabel::OBMol &mol);
/*
 * 2
 * #HA: number of hydrogen bond acceptors
 * #HD: number of hydrogen bond donors
 */
std::vector<double> DESCRIPTOR_API hbonds(OpenBabel::OBMol &mol);
#endif   /* ----- #ifndef CONSTITUTIONAL_H  ----- */

