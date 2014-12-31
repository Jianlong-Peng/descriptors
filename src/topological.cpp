/*=============================================================================
#     FileName: topological.cpp
#         Desc: topological descriptors
#       Author: jlpeng
#        Email: jlpeng1201@gmail.com
#     HomePage: 
#      Version: 0.0.1
#      Created: 2012-12-28 15:39:24
#   LastChange: 2014-03-20 17:08:01
#      History:
=============================================================================*/
#include <iostream>
#include <sstream>
#include <vector>
#include <numeric>
#include <string>
#include <map>
#include <bitset>
#include <utility>
#include <queue>
#include <iterator>
#include <cmath>
#include <cstdlib>
#include "topological.h"
#include "graph.h"
#include "tools.h"
#include "dmath.h"
#include <openbabel/mol.h>
#include <openbabel/graphsym.h>

using std::ostringstream;
using std::vector;
using std::accumulate;
using std::max_element;
using std::string;
using std::map;
using std::pair;
using std::make_pair;
using std::bitset;
using std::queue;
using namespace OpenBabel;

/*
static string str_names[] = {"Wiener","Randic0","Randic1","Randic2","Randic3",
    "Kier&Hall0","Kier&Hall1","Kier&Hall2","Kier&Hall3",
    "KierFlex","KierShape1","KierShape2","KierShape3","Balaban",
    "averIC(order0)","IC(order0)","averSIC(order0)","SIC(order0)","averCIC(order0)","CIC(order0)",
    "averBIC(order0)","BIC(order0)",
    "averIC(order1)","IC(order1)","averSIC(order1)","SIC(order1)","averCIC(order1)","CIC(order1)",
    "averBIC(order1)","BIC(order1)",
    "averIC(order2)","IC(order2)","averSIC(order2)","SIC(order2)","averCIC(order2)","CIC(order2)",
    "averBIC(order2)","BIC(order2)"};
vector<string> DESCRIPTOR_API topology_names(str_names,str_names+38);
*/

double wiener(OBMol &mol)
{
    double value = 0.;
    FOR_ATOMS_OF_MOL(atom,mol) {
        if(atom->IsHydrogen())
            continue;
        vector<double> path = shortestDist(mol,atom->GetIdx());
        value += accumulate(path.begin(),path.end(),0.);
    }
    return (0.5*value);
}

// calculate randic index (order 0~3)
void randic(OBMol &mol, double value[4])
{
    // 1. calculate degree of each vertex and Randic index order 0
    vector<double> degree(mol.NumAtoms(),0.);
    value[0] = 0.;
    FOR_ATOMS_OF_MOL(atom,mol) {
        if(atom->IsHydrogen())
            continue;
        FOR_NBORS_OF_ATOM(nbor,*atom)
            if(!nbor->IsHydrogen()) degree[atom->GetIdx()-1] += 1.;
        if(degree[atom->GetIdx()-1] == 0.)
            value[0] = 1E+38;
        else
            value[0] += (1.0/sqrt(degree[atom->GetIdx()-1]));
    }
    // 2. order 1~3
    for(int i=1; i<=3; ++i) {
        value[i] = 0.;
        vector<vector<unsigned> > path = pathOfLenN(mol,i);
        if(path.empty()) {
            value[i] = 1E+38;
            continue;
        }
        //std::cout << "path of length " << i << std::endl;
        for(vector<vector<unsigned> >::size_type j=0; j<path.size(); ++j) {
            double tmp = 1.;
            for(vector<unsigned>::size_type k=0; k<path[j].size(); ++k)
                tmp *= degree[path[j][k]-1];
            //copy(path[j].begin(),path[j].end(),std::ostream_iterator<unsigned>(std::cout," "));
            //std::cout << std::endl;
            value[i] += (1.0/sqrt(tmp));
        }
    }
    //std::cout << "~~~~~~~~~~~~~~END~~~~~~~~~~~~~~~~" << std::endl;
}

// calcualte Kier&Hall index (order 0~3)
void kier_hall(OBMol &mol, double value[4])
{
    // 1. calculate atomic connectivity and index of order 0
    vector<double> connectivity(mol.NumAtoms(),0);
    value[0] = 0.;
    FOR_ATOMS_OF_MOL(atom,mol) {
        if(atom->IsHydrogen())
            continue;
        int hydroNbors = 0;
        FOR_NBORS_OF_ATOM(nbor,*atom) {
            if(nbor->IsHydrogen())
                ++hydroNbors;
        }
        int val_electrons = getValElectrons(*atom);
        if(val_electrons == -1) {
            for(int i=0; i<4; ++i)
                value[i] = 1E+38;
            return ;
        }
        connectivity[atom->GetIdx()-1] = 
            ((double)(val_electrons-hydroNbors)/(atom->GetAtomicNum()-val_electrons-1));
        if(connectivity[atom->GetIdx()-1] == 0.)
            value[0] = 1E+38;
        else
            value[0] += (1.0/sqrt(connectivity[atom->GetIdx()-1]));
    }

    // 3. order 1~3
    for(int i=1; i<=3; ++i) {
        value[i] = 0.;
        vector<vector<unsigned> > path = pathOfLenN(mol,i);
        if(path.empty()) {
            value[i] = 1E+38;
            continue;
        }
        for(vector<vector<unsigned> >::size_type j=0; j<path.size(); ++j) {
            double tmp = 1.;
            for(vector<unsigned>::size_type k=0; k<path[j].size(); ++k)
                tmp *= connectivity[path[j][k]-1];
            value[i] += (1.0/sqrt(tmp));
        }
    }
}

// Kier shape index (order 1~3)
// value[0]: Kier flexibility index
void kier_shape(OBMol &mol, double value[4])
{
    // 1. number of heavy atoms
    int numHvyAtoms = 0;
    FOR_ATOMS_OF_MOL(atom,mol) {
        if(!atom->IsHydrogen())
            ++numHvyAtoms;
    }
    // 2. alpha
    double alpha = 0.;
    bool success(true);
    FOR_ATOMS_OF_MOL(atom,mol) {
        if(atom->IsHydrogen())
            continue;
        double radius = getCovalentRadius(*atom);
        if(radius == 0.) {
            success = false;
            break;
        }
        alpha += (radius / 0.77 - 1);
    }
    if(!success) {
        for(int i=0; i<4; ++i)
            value[i] = 1E+38;
        return ;
    }
    // 3. calculate Kier shape index
    // - order 1
    vector<vector<unsigned> > path = pathOfLenN(mol,1);
    if(path.empty()) {
        for(int i=0; i<4; ++i)
            value[i] = 1E+38;
        return ;
    }
    else
        value[1] = (numHvyAtoms+alpha)*pow(numHvyAtoms+alpha-1,2)/pow(path.size()+alpha,2);
    // - order 2
    path.clear();
    path = pathOfLenN(mol,2);
    if(path.empty()) {
        for(int i=2; i<4; ++i)
            value[i] = 1E+38;
        value[0] = 1E+38;
        return ;
    }
    else
        value[2] = (numHvyAtoms+alpha-1)*pow(numHvyAtoms+alpha-2,2)/pow(path.size()+alpha,2);
    // - order 3
    path.clear();
    path = pathOfLenN(mol,3);
    if(path.empty())
        value[3] = 1E+38;
    else {
        if(numHvyAtoms&1)
            value[3] = (numHvyAtoms+alpha-1)*pow(numHvyAtoms+alpha-3,2)/pow(path.size()+alpha,2);
        else
            value[3] = (numHvyAtoms+alpha-3)*pow(numHvyAtoms+alpha-2,2)/pow(path.size()+alpha,2);
    }
    // 4. Kier flexibility index
    value[0] = (value[1]*value[2])/numHvyAtoms;
}


// Balaban index
double balaban(OBMol &mol, bool codessa)
{
    // 1. calculate distance(shortest path) matrix
    // CODESSA use weighted shortest distance!!!!!!
    vector<vector<double> > distance;
    double numHvyAtoms = 0.;
    FOR_ATOMS_OF_MOL(atom,mol) {
        if(atom->IsHydrogen())
            distance.push_back(vector<double>());
        else {
            if(codessa)
                distance.push_back(weightedShortestDist(mol,atom->GetIdx()));
            else
                distance.push_back(shortestDist(mol,atom->GetIdx()));
            numHvyAtoms += 1.;
        }
    }
    // 2. calculate Balaban index
    double value=0., numHvyBonds=0.;
    FOR_BONDS_OF_MOL(bond, mol) {
        OBAtom *atom1 = bond->GetBeginAtom();
        OBAtom *atom2 = bond->GetEndAtom();
        if(atom1->IsHydrogen() || atom2->IsHydrogen())
            continue;
        unsigned i = atom1->GetIdx();
        unsigned j = atom2->GetIdx();
        double si = static_cast<double>(accumulate(distance[i-1].begin(),distance[i-1].end(),0.));
        double sj = static_cast<double>(accumulate(distance[j-1].begin(),distance[j-1].end(),0.));
        //value += pow(si*sj,-0.5);
        value += 1./sqrt(si*sj);
        numHvyBonds += 1.;
    }
    value *= (numHvyBonds/(numHvyBonds-numHvyAtoms+1+1));

    return value;
}


#if defined(_MSC_VER)
static inline double log2(double value)
{
    return log(value)/log(2.);
}
#endif


// do not include {n}th layer
// start with {0}th layer
static void refineNodeIdentifier(OBMol &mol, vector<int> &values, int n)
{
    if(values.empty()) {
        values.resize(mol.NumAtoms());
        for(unsigned i=0; i<mol.NumAtoms(); ++i)
            values[i] = static_cast<int>(mol.GetAtom(i+1)->GetAtomicNum());
    }
    vector<int> new_values(values.begin(),values.end());
    for(unsigned i=0; i<mol.NumAtoms(); ++i) {
        int temp = 0;
        queue<unsigned> toBeVisit;
        queue<unsigned> layer;
        bitset<200> visited;
        toBeVisit.push(i+1);
        layer.push(0);
        visited.set(i);
        while(!toBeVisit.empty()) {
            unsigned currentId = toBeVisit.front();
            unsigned currentLayer = layer.front();
            // add node identifier
            temp += values[currentId-1];
            FOR_NBORS_OF_ATOM(nbor,mol.GetAtom(currentId)) {
                unsigned nborId = nbor->GetIdx();
                if(visited.test(nborId-1))
                    continue;
                if(currentLayer < n) {
                    toBeVisit.push(nborId);
                    layer.push(currentLayer+1);
                    visited.set(nborId-1);
                }
                OBBond *bond = mol.GetBond(mol.GetAtom(currentId),mol.GetAtom(nborId));
                // add bond order
                temp += bond->GetBondOrder();
            }
            toBeVisit.pop();
            layer.pop();
        }
        new_values[i] = temp;
    }
    for(vector<int>::size_type i=0; i<values.size(); ++i)
        values[i] = new_values[i];
}

// strategy:
// 1. [atomicNumber, numTriple, numDouble, numSingle] ---hash--> integer
// 2. [atomicNumber, numTriple, numDouble, numSingle, numAromatic]  --hash--> integer
// currently, simply using 10^i*n to map the array!
static void assignNodeIdentifier(OBMol &mol, vector<int> &values)
{
    if(values.empty()) {
        values.resize(mol.NumAtoms());
        for(unsigned i=0; i<mol.NumAtoms(); ++i) {
            values[i] = static_cast<int>(mol.GetAtom(i+1)->GetAtomicNum()) * 1000;
            FOR_NBORS_OF_ATOM(nbor,mol.GetAtom(i+1)) {
                OBBond *bond = mol.GetBond(mol.GetAtom(i+1),&*nbor);
                switch(bond->GetBondOrder()) {
                    case 1: values[i] += 1; break;
                    case 2: values[i] += 10; break;
                    case 3: values[i] += 100; break;
                    default: std::cerr << "Error: invalid bond order " << bond->GetBondOrder() << std::endl; exit(EXIT_FAILURE);
                }
                /*
                if(bond->IsAromatic())
                    values[i] += 1;
                else if(bond->IsSingle())
                    values[i] += 10;
                else if(bond->IsDouble())
                    values[i] += 100;
                else if(bond->IsTriple())
                    values[i] += 1000;
                else {
                    std::cerr << "Error: invalid bond order " << bond->GetBondOrder() << std::endl;
                    exit(EXIT_FAILURE);
                }
                */
            }
        }
    } else {
        vector<int> new_values(values.begin(),values.end());
        for(unsigned i=0; i<values.size(); ++i) {
            FOR_NBORS_OF_ATOM(nbor,mol.GetAtom(i+1))
                values[i] += new_values[nbor->GetIdx()-1];
        }
    }
}
static double calcAverIC(OBMol &mol, vector<int> &node_values)
{
    map<int,int> type_count;  // key= "atomicNum_identify", value=frequency
    for(vector<int>::size_type i=0; i<node_values.size(); ++i) {
        //std::cout << " " << node_values[i];
        pair<map<int,int>::iterator,bool> iter
            = type_count.insert(make_pair(node_values[i],1));
        if(!iter.second)
            iter.first->second += 1;
    }
    //std::cout << std::endl;
    double value = 0.;
    for(map<int,int>::iterator iter=type_count.begin(); iter!=type_count.end(); ++iter) {
        double temp = static_cast<double>(iter->second) / node_values.size();
        value -= temp * log2(temp);
    }
    return value;
}
void infoContent(OBMol &mol, double value[3][8])
{
    unsigned natoms = mol.NumAtoms();
    double log2_natoms = log2(static_cast<double>(natoms));
    // 1. calculate average information content
    vector<int> node_values;
    for(int i=0; i<3; ++i) {
        assignNodeIdentifier(mol,node_values);
        //refineNodeIdentifier(mol,node_values,i);
        //for(vector<int>::iterator iter=node_values.begin(); iter!=node_values.end(); ++iter)
        //    std::cout << " " << *iter;
        //std::cout << std::endl;
        value[i][0] = calcAverIC(mol,node_values);
    }
    // 2. calculate other information content indices
    for(int i=0; i<3; ++i) {
        value[i][1] = value[i][0] * natoms;        // IC
        value[i][2] = value[i][0] / log2_natoms;   // averSIC
        value[i][3] = value[i][2] * natoms;
        //value[i][3] = value[i][1] / log2_natoms;   // SIC
        value[i][4] = log2_natoms - value[i][0];   // averCIC
        value[i][5] = value[i][4] * natoms;
        //value[i][5] = log2_natoms - value[i][1];   // CIC
        value[i][6] = value[i][0] / log2(static_cast<double>(mol.NumBonds())); // averBIC
        value[i][7] = value[i][6] * natoms;
        //value[i][7] = value[i][1] / log2(static_cast<double>(mol.NumBonds())); // BIC
    }
}

/*
// static method to calculate average information content
// called by infoContent
static double averIC(map<string,double> &type_count, unsigned natoms)
{
    double value = 0.;
    for(map<string,double>::iterator iter=type_count.begin(); iter!=type_count.end(); ++iter) {
        double temp = iter->second / natoms;
        value -= temp * log2(temp);
    }
    return value;
}
// information content index
// averIC, IC, averSIC, SIC, averCIC, CIC, averBIC, BIC (order 0~2)
// reference
//   Basak, S.C. et al. J. Pharm.Sci. 1984. 73: 429-437.
//   DOI: 10.1002/jps.2600730403
void infoContent(OBMol &mol, double value[3][8])
{
    map<string,double> type_count; // type_count.first: "atomicNum_bondorder_..."
    unsigned natoms = mol.NumAtoms();
    double log2_natoms = log2(static_cast<double>(natoms));
    // average information content (order 0)
    FOR_ATOMS_OF_MOL(atom,mol) {
        unsigned atomic_number = atom->GetAtomicNum();
        ostringstream os;
        os << atomic_number;
        // get bond orders
        vector<unsigned> bondorders;
        FOR_BONDS_OF_ATOM(bond,*atom)
            bondorders.push_back(bond->GetBondOrder());
        sort(bondorders.begin(),bondorders.end());
        for(vector<unsigned>::iterator iter=bondorders.begin(); iter!=bondorders.end(); ++iter)
            os << "_" << *iter;
        pair<map<string,double>::iterator,bool> iter
            = type_count.insert(make_pair(os.str(),1));
        if(!iter.second)
            iter.first->second += 1.;
    }
    value[0][0] = averIC(type_count, natoms);  // averIC
    // average information content (order 1)
    for(int j=0; j<8; ++j)
        value[1][j] = 1E+38;
    // average information content (order 2)
    type_count.clear();
    vector<unsigned> topoSymClass;
    OBGraphSym sym(&mol);
    sym.GetSymmetry(topoSymClass);
    for(vector<unsigned>::iterator i=topoSymClass.begin(); i!=topoSymClass.end(); ++i) {
        ostringstream os;
        os << *i;
        pair<map<string,double>::iterator,bool> iter
            = type_count.insert(make_pair(os.str(),1));
        if(!iter.second)
            iter.first->second += 1.;
    }
    value[2][0] = averIC(type_count,natoms);  // averIC
    // calculate other information content indices
    for(int i=0; i<3; ++i) {
        if(i == 1)
            continue;
        value[i][1] = value[i][0] * natoms;        // IC
        value[i][2] = value[i][0] / log2_natoms;   // averSIC
        value[i][3] = value[i][2] * natoms;
        //value[i][3] = value[i][1] / log2_natoms;   // SIC
        value[i][4] = log2_natoms - value[i][0];   // averCIC
        value[i][5] = value[i][4] * natoms;
        //value[i][5] = log2_natoms - value[i][1];   // CIC
        value[i][6] = value[i][0] / log2(static_cast<double>(mol.NumBonds())); // averBIC
        value[i][7] = value[i][6] * natoms;
        //value[i][7] = value[i][1] / log2(static_cast<double>(mol.NumBonds())); // BIC
    }
}
*/

vector<double> calcEState(OBMol &mol)
{
    // matrix of intrinsic state of heavy atoms
    vector<double> I(mol.NumAtoms());
    for(unsigned i=0; i<mol.NumAtoms(); ++i) {
        OBAtom *atom = mol.GetAtom(i+1);
        if(atom->IsHydrogen())
            I[i] = -1.;
        else {
            int n = getPrincipalQuantumNumber(*atom);
            if(n == -1)
                return vector<double>();
            unsigned delta = atom->GetHvyValence();
            int valElectrons = getValElectrons(*atom);
            if(valElectrons == -1)
                return vector<double>();
            int deltav = valElectrons - getNumOfHNbors(*atom);
            I[i] = (pow(2./n, 2) * deltav + 1.) / delta;
        }
    }
    vector<vector<int> > distance = floyd(mol);
    vector<double> estate(mol.NumAtoms(),-1.);
    for(unsigned i=0; i<mol.NumAtoms(); ++i) {
        if(I[i]==-1.)
            continue;
        double value = I[i];
        for(unsigned j=0; j<mol.NumAtoms(); ++j) {
            if(j==i || mol.GetAtom(j+1)->IsHydrogen())
                continue;
            value += ((I[i]-I[j])/pow(static_cast<double>(distance[i][j]+1),2));
        }
        estate[i] = value;
    }
    return estate;    
}

vector<double> calcSPAN(OBMol &mol)
{
    vector<vector<int> > distance = floyd(mol);
    int max_distance = 0;
    vector<double> span(mol.NumAtoms());
    for(unsigned i=0; i<mol.NumAtoms(); ++i) {
        int max_distance_i = *max_element(distance[i].begin(),distance[i].end());
        span[i] = 1.0 * max_distance_i;
        if(max_distance_i > max_distance)
            max_distance = max_distance_i;
    }
    for(vector<double>::size_type i=0; i<span.size(); ++i)
        span[i] /= max_distance;
    return span;    
}

// to construct molecular matrix constituted by A rows and 3 columns
// centring the molecule
static double **obmol2mm(OBMol &mol)
{
    unsigned A = mol.NumAtoms();
    double **m = (double**)malloc(sizeof(double*)*A);
    double center[3] = {0.,0.,0.};
    for(unsigned i=1; i<=A; ++i) {
        m[i-1] = (double*)malloc(sizeof(double)*3);
        double *coord = mol.GetAtom(i)->GetCoordinate();
        m[i-1][0] = coord[0];
        m[i-1][1] = coord[1];
        m[i-1][2] = coord[2];
        center[0] += coord[0];
        center[1] += coord[1];
        center[2] += coord[2];
    }
    center[0] /= A; center[1] /= A; center[2] /= A;
    for(unsigned i=1; i<=A; ++i) {
        m[i-1][0] -= center[0];
        m[i-1][1] -= center[1];
        m[i-1][2] -= center[2];
    }
    return m;
}
// value[]: H_GM, I_TH, I_SH, HIC
void getaway(OBMol &mol, double value[])
{
    double **M = obmol2mm(mol);
    double **H = leverage(M, mol.NumAtoms(), 3);
    int n;
    double D = 0.;  // matrix rank
    for(unsigned i=0; i<mol.NumAtoms(); ++i)
        D += H[i][i];
    // 1. H_GM
    value[0] = 1.;
    n = 0;
    for(unsigned i=0; i<mol.NumAtoms(); ++i) {
        if(H[i][i] <= 1e-6)
            continue;
        value[0] *= H[i][i];
        ++n;
    }
    value[0] = 100. * value[0] * pow(value[0], 1./n);
    // 2. I_TH, I_SH - total and standardized information content on the leverage equality
    map<double, int> Ng;
    int A0 = 0;
    for(unsigned i=0; i<mol.NumAtoms(); ++i) {
        if(mol.GetAtom(i+1)->IsHydrogen())
            continue;
        ++A0;
        pair<map<double,int>::iterator, bool> iter = Ng.insert(make_pair(H[i][i],1));
        if(!iter.second)
            ++iter.first->second;
    }
    double temp2 = 0.;
    for(map<double,int>::iterator iter=Ng.begin(); iter!=Ng.end(); ++iter)
        temp2 += (iter->first * log2(iter->first));
    double temp1 = A0 * log2(A0);
    value[1] = temp1 - temp2;      // I_TH
    value[2] = 1. - temp2 / temp1; // I_SH
    // 3. HIC - mean information content on the leverage magnitude
    value[3] = 0.;
    for(unsigned i=0; i<mol.NumAtoms(); ++i)
        value[3] += (-1. * H[i][i] / D * log2(H[i][i] / D));
    
}



