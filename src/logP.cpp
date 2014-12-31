/*=============================================================================
#     FileName: logP.cpp
#         Desc: to predict logP (Octanol-water partition coefficient)
#       Author: jlpeng
#        Email: jlpeng1201@gmail.com
#     HomePage: 
#      Version: 0.0.1
#      Created: 2013-01-24 21:21:48
#   LastChange: 2013-01-25 11:45:23
#      History:
=============================================================================*/
#include <iterator>
#include <cstdlib>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <numeric>
#include <map>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/graphsym.h>
#include "logP.h"

using std::string;
using std::vector;
using std::ifstream;
using std::cerr;
using std::endl;
using std::istringstream;
using std::accumulate;
using std::map;
using namespace OpenBabel;
/*
void AlogP::read_parameter(string para_file)
{
    frag_logP.clear();
    frag_logP.resize(6);

    ifstream inf(para_file.c_str());
    if(!inf) {
        cerr << "Error: can't open parameter file \"" << para_file << "\"" << endl;
        exit(EXIT_FAILURE);
    }

    int flag = -1;
    string line;
    while(getline(inf,line)) {
        if(line == "")
            continue;
        if(line[0] == '/') {
            if(line == "/hydrogen")
                flag = 0;
            else if(line == "/carbon")
                flag = 1;
            else if(line == "/nitrogen")
                flag = 2;
            else if(line == "/oxygen")
                flag = 3;
            else if(line=="/fluorine" || line=="/chlorine" 
                    || line=="/bromine" || line=="/iodine")
                flag = 4;
            else if(line=="/sulfur" || line=="/phosphorus")
                flag = 5;
            else
                continue;
            getline(inf,line);
        }
        while(line!="") {
            if(line[0] == '/') {
                getline(inf,line);
                continue;
            }
            istringstream is(line);
            string pattern;
            double value;
            is >> pattern >> value;
            frag_logP[flag].push_back(make_pair(pattern,value));
            getline(inf,line);
        }
    }
    inf.close();
    cerr << endl << "after reading parameter file" << endl
        << "number of atoms of each type:" << endl;
    for(int i=0; i<6; ++i) {
        cerr << " " << frag_logP[i].size() << endl;
    }
}
*/
AlogP::~AlogP()
{
    if(frag_logP != NULL) {
        for(int i=0; i<num; ++i) {
            if(frag_logP[i].last==NULL)
                continue;
            struct AlogPRecord *p = frag_logP[i].record;
            while(p != NULL) {
                struct AlogPRecord *q = p;
                p = p->next;
                delete q;
            }
            frag_logP[i].last = NULL;
        }
        delete[] frag_logP;
    }
}

void AlogP::read_parameter(string para_file)
{
    frag_logP = new AlogPHead[12];
    num = 12;
    for(int i=0; i<12; ++i) {
        frag_logP[i].last = NULL;
        frag_logP[i].record = NULL;
    }

    ifstream inf(para_file.c_str());
    if(!inf) {
        cerr << "Error: can't open parameter file \"" << para_file << "\"" << endl;
        exit(EXIT_FAILURE);
    }

    int flag = -1;
    string line;
    while(getline(inf,line)) {
        if(line == "")
            continue;
        if(line[0]=='/') {
            if (line.size()>6 && line.substr(0,5)=="/type") {
                line = line.substr(6);
                if(line == "hydrogen")
                    flag = 0;
                else if(line == "SP3 carbon")
                    flag = 1;
                else if(line == "SP2 carbon")
                    flag = 2;
                else if(line == "SP carbon")
                    flag = 3;
                else if(line == "aromatic carbon")
                    flag = 4;
                else if(line == "nitrogen")
                    flag = 5;
                else if(line == "oxygen")
                    flag = 6;
                else if(line == "fluorine")
                    flag = 7;
                else if(line == "chlorine")
                    flag = 8;
                else if(line == "bromine")
                    flag = 9;
                else if(line == "iodine")
                    flag = 10;
                else
                    flag = 11;
                getline(inf,line);
            } else
                continue;
        }
        while(line!="") {
            if(line[0] == '/') {
                getline(inf,line);
                continue;
            }
            istringstream is(line);
            struct AlogPRecord *new_record = new AlogPRecord;
            is >> new_record->pattern >> new_record->value;
            new_record->next = NULL;
            if(frag_logP[flag].last == NULL) {
                frag_logP[flag].record = new_record;
                frag_logP[flag].last = new_record;
            } else {
                (frag_logP[flag].last)->next = new_record;
                frag_logP[flag].last = new_record;
            }
            getline(inf,line);
        }
    }
    inf.close();
}

double AlogP::sub_predict(OBAtom &atom, int i)
{
    struct AlogPRecord *p = frag_logP[i].record;
    while(p != NULL) {
        if(atom.MatchesSMARTS((p->pattern).c_str()))
            return p->value;
        p = p->next;
    }
    return 0.;
}

double AlogP::predict(OBMol &mol)
{
    double *values = new double[mol.NumAtoms()];
    memset(values,0,sizeof(double)*(mol.NumAtoms()));
    FOR_ATOMS_OF_MOL(atom,mol) {
        int type = -1;
        unsigned atomic_number = atom->GetAtomicNum();
        if(atomic_number == 1)
            type = 0;
        else if(atomic_number == 6) {
            if(atom->MatchesSMARTS("[C^3]"))
                type = 1;
            else if(atom->MatchesSMARTS("[C^2]"))
                type = 2;
            else if(atom->MatchesSMARTS("[C^1]"))
                type = 3;
            else if(atom->IsAromatic())
                type = 4;
            else
                type = -1;
        }
        else if(atomic_number == 7)
            type = 5;
        else if(atomic_number == 8)
            type = 6;
        else if(atomic_number==9)
            type = 7;
        else if(atomic_number==17)
            type = 8;
        else if(atomic_number==35)
            type = 9;
        else if(atomic_number==53)
            type = 10;
        else
            type = 11;

        if(type != -1)
            values[atom->GetIdx()-1] = sub_predict(*atom,type);
    }

    /*
    cerr << endl << "predict ";
    copy(values,values+mol.NumAtoms(),std::ostream_iterator<double>(cerr," "));
    cerr << endl;
    */
    double value = accumulate(values,values+mol.NumAtoms(),0.);
    delete[] values;
    return value;
}

double AlogP::predict_sym(OBMol &mol)
{
    double *values = new double[mol.NumAtoms()];
    memset(values,0,sizeof(double)*(mol.NumAtoms()));

    // get topological symmetry
    vector<unsigned> topoSymClass;
    OBGraphSym sym(&mol);
    sym.GetSymmetry(topoSymClass);
    /*
    cerr << endl << "sym ";
    copy(topoSymClass.begin(),topoSymClass.end(),std::ostream_iterator<unsigned>(cerr," "));
    cerr << endl << "atom ";
    for(unsigned i=1; i<=mol.NumAtoms(); ++i)
        cerr << " " << (mol.GetAtom(i))->GetIdx() << "(" << (mol.GetAtom(i))->GetType() << ")";
    cerr << endl;
    */
    map<unsigned,vector<unsigned> > unique_topoSymClass;
    for(vector<unsigned>::size_type i=0; i<topoSymClass.size(); ++i) {
        if(unique_topoSymClass.find(topoSymClass[i]) == unique_topoSymClass.end())
            unique_topoSymClass.insert(make_pair(topoSymClass[i],vector<unsigned>()));
        unique_topoSymClass[topoSymClass[i]].push_back(i+1);
    }

    // to calculate atomic contribution to logP
    for(map<unsigned,vector<unsigned> >::size_type i=0; i<unique_topoSymClass.size(); ++i) {
        unsigned first_id = 0;
        for(vector<unsigned>::size_type j=0; j<unique_topoSymClass[i].size(); ++j) {
            unsigned atom_id = unique_topoSymClass[i][j];
            if(j == 0) {
                first_id = atom_id;
                OBAtom *atom = mol.GetAtom(atom_id);
                // find the atomic type -> 0~11
                int type = -1;
                unsigned atomic_number = atom->GetAtomicNum();
                if(atomic_number == 1)
                    type = 0;
                else if(atomic_number == 6) {
                    if(atom->MatchesSMARTS("[C^3]"))
                        type = 1;
                    else if(atom->MatchesSMARTS("[C^2]"))
                        type = 2;
                    else if(atom->MatchesSMARTS("[C^1]"))
                        type = 3;
                    else if(atom->IsAromatic())
                        type = 4;
                    else
                        type = -1;
                }
                else if(atomic_number == 7)
                    type = 5;
                else if(atomic_number == 8)
                    type = 6;
                else if(atomic_number==9)
                    type = 7;
                else if(atomic_number==17)
                    type = 8;
                else if(atomic_number==35)
                    type = 9;
                else if(atomic_number==53)
                    type = 10;
                else
                    type = 11;

                if(type != -1)
                    values[atom_id-1] = sub_predict(*atom,type);
            }
            else
                values[atom_id-1] = values[first_id-1];
        }
    }

    /*
    cerr << endl << "predict_sym ";
    copy(values,values+mol.NumAtoms(),std::ostream_iterator<double>(cerr," "));
    cerr << endl;
    */
    double value = accumulate(values,values+mol.NumAtoms(),0.);
    delete[] values;
    return value;
}

