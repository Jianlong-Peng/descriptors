/*=============================================================================
#     FileName: constitutional.cpp
#         Desc: 
#       Author: jlpeng
#        Email: jlpeng1201@gmail.com
#     HomePage: 
#      Version: 0.0.1
#      Created: 2012-12-28 14:32:38
#   LastChange: 2014-01-08 13:36:25
#      History:
=============================================================================*/
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include "constitutional.h"
#include <openbabel/mol.h>
#include <openbabel/ring.h>

using std::map;
using std::string;
using std::vector;
using namespace OpenBabel;

/*
static string str_names[] = {"natoms","nCatoms","rnCatoms","nHatoms","rnHatoms","nOatoms","rnOatoms",
        "nNatoms","rnNatoms","nSatoms","rnSatoms","nFatoms","rnFatoms","nClatoms","rnClatoms",
        "nBratoms","rnBratoms","nIatoms","rnIatoms","nPatoms","rnPatoms", "MW", "rMW",
        "nbonds","nsbonds","rnsbond","ndbonds","rndbonds","ntbonds","rntbonds","narbonds","rnarbonds",
        "nrings","rnrings","nbrings","rnbrings"};
vector<string> constitute_names(str_names,str_names+36);
*/

vector<double> atoms(OBMol &mol)
{
    vector<double> features(23,0.);
    FOR_ATOMS_OF_MOL(atom,mol) {
        features[0] += 1.;
        int  i;
        if(atom->IsCarbon())
            i = 1;
        else if(atom->IsHydrogen())
            i = 3;
        else if(atom->IsOxygen())
            i = 5;
        else if(atom->IsNitrogen())
            i = 7;
        else if(atom->IsSulfur())
            i = 9;
        else if(atom->GetAtomicNum() == 9)
            i = 11;
        else if(atom->GetAtomicNum() == 17)
            i = 13;
        else if(atom->GetAtomicNum() == 35)
            i = 15;
        else if(atom->GetAtomicNum() == 53)
            i = 17;
        else if(atom->IsPhosphorus())
            i = 19;
        else
            i = -1;
        if(i != -1)
            features[i] += 1.;
    }
    features[21]  = mol.GetMolWt();
    for(int i=2; i<23; i+=2)
        features[i] = features[i-1] / features[0];
    return features;
}

vector<double> bonds(OBMol &mol)
{
    vector<double> features(9,0.);
    FOR_BONDS_OF_MOL(bond,mol) {
        features[0] += 1.;
        int i;
        if(bond->IsSingle())
            i = 1;
        else if(bond->IsDouble())
            i = 3;
        else if(bond->IsTriple())
            i = 5;
        else if(bond->IsAromatic())
            i = 7;
        else
            i = -1;
        if(i != -1)
            features[i] += 1.;
    }
    for(int i=2; i<9; i+=2)
        features[i] = features[i-1] / features[0];
    return features;
}

vector<double> rings(OBMol &mol)
{
    vector<double> features(4);
    vector<OBRing*> sssr = mol.GetSSSR();
    features[0] = static_cast<double>(sssr.size());
    features[1] = features[0] / mol.NumAtoms();
    features[2] = 0.;
    for(vector<OBRing*>::iterator iter=sssr.begin(); iter!=sssr.end(); ++iter) {
        if(strcmp((*iter)->GetType(),"benzene")==0)
            features[2] += 1.;
    }
    features[3] = features[2]/mol.NumAtoms();
    return features;
}

vector<double> hbonds(OBMol &mol)
{
    vector<double> features(2,0.);
    FOR_ATOMS_OF_MOL(atom,mol) {
        if(atom->MatchesSMARTS("[!$([#1,#6,F,Cl,Br,I,o,s,nx3,#7v5,#15v5,#16v4,#16v6])]"))
            features[0] += 1.;
        else if(atom->MatchesSMARTS("[#7,#8;!H0]"))
            features[1] += 1.;
        else
            continue;
    }
    return features;
}


