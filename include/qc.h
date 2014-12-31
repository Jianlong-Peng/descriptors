/*=============================================================================
#     FileName: qc.h
#         Desc: quantum-chemical features
#       Author: jlpeng
#        Email: jlpeng1201@gmail.com
#     HomePage: 
#      Version: 0.0.1
#      Created: 2012-12-30 18:20:19
#   LastChange: 2014-02-13 19:57:08
#      History:
=============================================================================*/
#ifndef  QC_H
#define  QC_H
#include <string>
#include <fstream>
#include <vector>
#include <openbabel/mol.h>
#include "config.h"
using std::string;
using std::ifstream;
using std::vector;

// "final heat of formation" and "No. of filled levels" can be read directly from a mopout file

// to read HOMO and LUMO coefficients.
// inf: opened mopout file. "EIGENVECTORS" shouldn't been read!
// homo_level: No. of filled levels
// fmo_energy: energy of HOMO-1, HOMO, LUMO, LUMO-1
// homo_coeff, lumo_coeff: where HOMO and LUMO coefficients will be saved.
//     .size() == no of atoms
//     [i].size() == no of valence atomic orbitals of atom i.
// return value
//   true - successfully read coefficients
//   false - otherwise
bool DESCRIPTOR_API read_fmo_coeff(ifstream &inf, int homo_level, double fmo_energy[4], 
        vector<vector<double> > &homo_coeff, vector<vector<double> > &lumo_coeff);
// to calculate nucleophilic reactivity indices for each atom
// homo_coeff: HOMO coefficients, size = number of atoms
// energy: HOMO energy
// return value:
//   nucelophilic reactivity indices of each atom
vector<double> DESCRIPTOR_API calcNA(vector<vector<double> > &homo_coeff, double energy);
// to calculate electrophilic reactivity indices for each atom
// lumo_coeff: LUMO coefficients, size = number of atoms
// energy: LUMO energy
// return value:
//   electrophilic reactivity indices of each atom
vector<double> DESCRIPTOR_API calcEA(vector<vector<double> > &lumo_coeff, double energy);
// to calculate 1-electron reactivity indices for each atom
// homo_coeff, lumo_coeff: HOMO/LUMO coefficients, size = number of atoms
// homo_energy, lumo_energy: HOMO/LUMO energies
// return value:
//   1-electron reactivity indices of each atom
vector<double> DESCRIPTOR_API calcRA(vector<vector<double> > &homo_coeff,
        vector<vector<double> > &lumo_coeff, double homo_energy, double lumo_energy);
// to calculate net atomic charge (Mulliken)
// inf: opened mopout file, "NET ATOMIC CHARGES" shouldn't been read before!!
// natoms: number of atoms;
// result: where net atomic charges saved.
// return value:
//     true - if sucess, otherwise false
bool DESCRIPTOR_API netAtomicChargeMulliken(ifstream &inf, int natoms, vector<double> &result);

// to extract atomic orbital electronic population
// inf: opened mopout file, "ATOMIC ORBITAL ELECTRON POPULATIONS" shouldn't been read before !!
// result: where atomic orbital electronic population saved
bool DESCRIPTOR_API getElectronPopulation(ifstream &inf, vector<double> &result);
// to calculate max SIGMA-SGIMA, SIGMA-PI, and PI-PI bond order
// inf: opened mopout file, "SIGMA-PI BOND-ORDER MATRIX" shouldn't been read before !!
// atom_symbol: atomic symbol
// result: where final result will be saved
//   [0,1,2]: each for max SIGMA-SIGMA, SIGMA-PI, and PI-PI bond order
bool DESCRIPTOR_API sigmaPiBond(ifstream &inf, vector<string> &atom_symbol, double result[3]);
// to calculate charged partial surface area [Quantum-chemical PC]
// Please use 'calcCPSA' decleared in 'geometrical.h'!

// to get energy partition terms in AM1 (one-center & two-center)
bool energyPartitionTerms(ifstream &inf, int natoms, double result[7][2]);
// to extract electrostatic potential charges
// "ELECTROSTATIC POTENTIAL CHARGES" shouldn't been read before!
bool netAtomicChargeESP(ifstream &inf, int natoms, vector<double> &result);
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */


// this structure was used by 'qcfeatures' to hold all QC features
struct DESCRIPTOR_API QCValues {
    double heat_formation[2]; // final heat of formation; */(# of atoms)
    double occ_levels[2]; // No. of occupied electronic levels; */(# of atoms)
    double fmo_energy[4]; // HOMO-1, HOMO, LUMO, LUMO+1
    // nucleophilic reactivity index
    // na[i][]: each for C, N, and O
    // na[][j]: each for min, max, and average
    double na[3][3];
    // electrophilic reactivity index
    double ea[3][3];
    // Fukui 1-electron reactivity index
    double ra[3][3];
    // net atomic charges (Mulliken)
    // charge[i][]: each for C, N, O, H, and all
    // charge[][j]: each for max and min
    double chargeMulliken[5][2];
    // Tot. point-charge comp. of the molecular dipole
    // Tot. hybridzation comp. of the molecular dipole
    // Tot. dipole of the molecule
    double dipole[3];
    // min/max atomic orbital electronic population
    double populationAO[2];
    // [i][]: SIGMA-SIGMA, SIGMA-PI, PI-PI bond order
    // [][j]: max, min
    double sigpiBond[3][2];
    // energyTerms[i][]:  E-N(1-center), E-E(1-center),
    //                    resonance, exchange, E-E(2-center), E-N(2-center), N-N(2-center)
    // energyTerms[][j]:  total, average
    double energyTerms[7][2];
    // electrostatic potential charges   (1e+38 if unavailable)
    // charge[i][]: each for C, N, O, H, and all
    // charge[][j]: each for max and min
    double chargeESP[5][2];
};

void DESCRIPTOR_API print_qc_result(std::ostream &os, const QCValues &result, std::string sep);

// a convenient function to calculate features, by calling the above functions
// mulliken, esp: to hold corresponding charges being extracted
// return:
//     true  - if success
//     false - otherwise
bool DESCRIPTOR_API qcfeatures(std::string mopout_file, QCValues &values, 
        vector<double> &mulliken, vector<double> &esp);
#endif   /* ----- #ifndef QC_H  ----- */

