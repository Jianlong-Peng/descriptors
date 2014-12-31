/*=============================================================================
#     FileName: electrostatic.cpp
#         Desc: 
#       Author: jlpeng
#        Email: jlpeng1201@gmail.com
#     HomePage: 
#      Version: 0.0.1
#      Created: 2013-01-08 11:05:08
#   LastChange: 2014-01-08 12:39:09
#      History:
=============================================================================*/
#include <vector>
#include <string>
#include <cassert>
#include <cmath>
#include <openbabel/mol.h>
#include "electrostatic.h"
#include "geometrical.h"
#include "tools.h"
#include "mymatrix.h"

using std::vector;
using std::string;
using namespace OpenBabel;

/*
static string temp_names[] = {"maxZefirovPCforH","minZefirovPCforH",
    "maxZefirovPCforC","minZefirovPCforC","maxZefirovPCforN","minZefirovPCforN",
    "maxZefirovPCforO","minZefirovPCforO","maxZefirovPCforF","minZefirovPCforF",
    "maxZefirovPCforP","minZefirovPCforP","maxZefirovPCforS","minZefirovPCforS",
    "maxZefirovPCforCl","minZefirovPCforCl","maxZefirovPCforBr","minZefirovPCforBr",
    "maxZefirovPCforI","minZefirovPCforI","maxZefirovPC","minZefirovPC",
    "polarity[Zefirov'sPC]","polarity/squareDistance[Zefirov'sPC]",
    "topolElectronicIndex(allPairs)[Zefirov'sPC]","topolElectronicIndex(allBonds)[Zefirov'sPC]"};
vector<string> DESCRIPTOR_API electrostatic_names(temp_names,temp_names+26);
*/

// reference
// 1. Oliferenko, etc. Doklady Chemistry, 2000, 375(5): 645-648
// 2. Oliferenko, etc. SAR and QSAR in Environmental Research, 2002, 13(2): 297-305
// 3. Alexander, etc. J. Phys. Org. Chem. 2001, 14: 355-369
// main idea:
// to solve the following equation: S*X = X0 ... (1)
// where, S is the special case of matrix of nodal conductivities, S=D-A+I
// (D: degree matrix of a graph, A: adjacency matrix of the graph, I: identity matrix)
// X is vector of unknown EN values, X0 is the vector of standard Pauling(Sanderson) EN values.
// after solving the equation (1), atomic partial charge can be calculated as
//                    q_i = (X_i-X0_i)/Xm
// where Xm is the geometrical mean of all atomic ENs (X0)
vector<double> calcZefirovPC(OBMol &mol, int scale_method)
{
    assert(scale_method==1 || scale_method==2);
    int natoms = static_cast<int>(mol.NumAtoms());
    double *en = new double [natoms];
    Descriptors::Matrix solution(natoms,natoms+1);
    FOR_ATOMS_OF_MOL(atom,mol) {
        int i = static_cast<int>(atom->GetIdx()) - 1;
        // adjacency matrix
        FOR_NBORS_OF_ATOM(nbor,*atom) {
            int j = static_cast<int>(nbor->GetIdx()) - 1;
            OBBond *b = mol.GetBond(&(*atom), &(*nbor));
            if(b->IsAromatic())
                solution(i,j) = -1.5;
            else
                solution(i,j) = -static_cast<double>(b->GetBondOrder());
        }
        // Sanderson electronegativity
        if(scale_method == 1)
            solution(i,natoms) = paulingEN(*atom);
        else
            solution(i,natoms) = sandersonEN(*atom);
        en[i] = solution(i,natoms);
    }
    // s = D-A+I
    for(int i=0; i<natoms; ++i) {
        double value(0.);
        for(int j=0; j<natoms; ++j)
            value += (-solution(i,j));
        solution(i,i) = value+1;
    }
    // solve linear equations
    double *answer = new double [natoms];
    if(!Descriptors::linearEquations(solution,answer)) {
        delete[] en;
        delete[] answer;
        return vector<double>();
    }
    // geometrical mean of original EN
    double gm = 1.;
    for(int i=0; i<natoms; ++i)
        gm *= en[i];
    gm = pow(gm,1.0/natoms);
    // partial charge
    vector<double> pc;
    for(int i=0; i<natoms; ++i)
        pc.push_back((answer[i]-en[i])/gm);

    delete[] en;
    delete[] answer;
    return pc;
}

bool calcMaxMinZefirovPC(OBMol &mol, vector<double> &pc, double result[12][2])
{
    for(int i=0; i<11; ++i) {
        result[i][0] = -1E+38; // max
        result[i][1] = +1E+38; // min
    }

    unsigned max_id=0, min_id=0;
    for(unsigned i=0; i<mol.NumAtoms(); ++i) {
        // max/min PC
        if(result[10][0] < pc[i]) {
            result[10][0] = pc[i];
            max_id = i+1;
        }
        if(result[10][1] > pc[i]) {
            result[10][1] = pc[i];
            min_id = i+1;
        }
        // max/min PC for a specific atom
        int id = -1;
        unsigned atomic_number = mol.GetAtom(i+1)->GetAtomicNum();
        if(atomic_number == 1)
            id = 0;
        if(atomic_number == 6)
            id = 1;
        else if(atomic_number == 7)
            id = 2;
        else if(atomic_number == 8)
            id = 3;
        else if(atomic_number == 9)
            id = 4;
        else if(atomic_number == 15)
            id = 5;
        else if(atomic_number == 16)
            id = 6;
        else if(atomic_number == 17)
            id = 7;
        else if(atomic_number == 35)
            id = 8;
        else if(atomic_number == 53)
            id = 9;
        else
            ;
        if(id == -1)
            continue;
        if(result[id][0] < pc[i])
            result[id][0] = pc[i];
        if(result[id][1] > pc[i])
            result[id][1] = pc[i];
    }
    result[11][0] = result[10][0] - result[10][1];
    OBAtom *max_atom = mol.GetAtom(max_id);
    OBAtom *min_atom = mol.GetAtom(min_id);
    double d = pos_distance_sq(max_atom->GetCoordinate(), min_atom->GetCoordinate());
    result[11][1] = result[11][0] / d;
    return true;
}

bool calcTopolElectronicIndex(OBMol &mol, vector<double> &pc, double result[2])
{
    result[0] = 0.;
    result[1] = 0.;
    // all pairs of atoms
    for(unsigned i=0; i<mol.NumAtoms()-1; ++i) {
        OBAtom *atomi = mol.GetAtom(i+1);
        for(unsigned j=i+1; j<mol.NumAtoms(); ++j) {
            OBAtom *atomj = mol.GetAtom(j+1);
            double d = pos_distance_sq(atomi->GetCoordinate(),atomj->GetCoordinate());
            result[0] += (abs(pc[i]-pc[j])/d);
        }
    }
    // all bonded pairs
    FOR_BONDS_OF_MOL(bond,mol) {
        OBAtom *batom = bond->GetBeginAtom();
        OBAtom *eatom = bond->GetEndAtom();
        double d = pos_distance_sq(batom->GetCoordinate(),eatom->GetCoordinate());
        result[1] += (abs(pc[batom->GetIdx()-1]-pc[eatom->GetIdx()-1])/d);
    }
    return true;
}


