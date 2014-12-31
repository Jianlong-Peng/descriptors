/*=============================================================================
#     FileName: qc.cpp
#         Desc: 
#       Author: jlpeng
#        Email: jlpeng1201@gmail.com
#     HomePage: 
#      Version: 0.0.1
#      Created: 2012-12-30 18:33:56
#   LastChange: 2014-02-13 20:26:09
#      History:
=============================================================================*/
#include <iterator>
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <cstdlib>
#include "qc.h"
#include "tools.h"
#include "geometrical.h"
using std::cout;
using std::endl;
using std::ostream;
using std::string;
using std::ifstream;
using std::vector;
using std::accumulate;
using std::find;

/* 
static string str_names[] = {"heat_formation","heat_formation/atoms","occ_electronicLevels",
    "occ_electronicLevels/atoms","HOMO-1","HOMO","LUMO","LUMO+1",
    "naCmin","naCmax","naCaver","naNmin","naNmax","naNaver","naOmin","naOmax","naOaver",
    "eaCmin","eaCmax","eaCaver","eaNmin","eaNmax","eaNaver","eaOmin","eaOmax","eaOaver",
    "raCmin","raCmax","raCaver","raNmin","raNmax","raNaver","raOmin","raOmax","raOaver",
    "chargeCmax","chargeCmin","chargeNmax","chargeNmin","chargeOmax","chargeOmin",
    "chargeHmax","chargeHmin","chargeMax","chargeMin",
    "dipole(point-chg.)","dipole(hybrid)","dipole(tot.)","minPopulation","maxPopulation",
    "maxSIGMA-SIGMA","maxSIGMA-PI","maxPI-PI",
    "totEN1","averEN1","totEE1","averEE1","totResonance2","averResonance2",
    "totExchange2","averExchange2","totEE2","averEE2","totEN2","averEN2","totNN2","averNN2",
    "chargeCmaxESP","chargeCminESP","chargeNmaxESP","chargeNminESP","chargeOmaxESP","chargeOminESP",
    "chargeHmaxESP","chargeHminESP","chargeMaxESP","chargeMinESP"};
vector<string> qc_names(str_names,str_names+77);
*/

// to calculate nucleophilic or electrophilic reactivity index
// maybe another function is better such that coefficients of all molecular orbitals were extracted!!!!!!
bool read_fmo_coeff(ifstream &inf, int homo_level, double fmo_energy[4], 
        vector<vector<double> > &homo_coeff, vector<vector<double> > &lumo_coeff)
{
    string line;
    while(getline(inf,line) && line.find("EIGENVECTORS")==string::npos)
        ;
    if(inf.eof()) {
        cout << endl << "<qc.cpp::read_fmo_coeff> Error: can't find \"EIGENVECTORS\"" << endl;
        return false;
    }
    getline(inf,line);  // empty line
    getline(inf,line);  // empty line
    getline(inf,line);
    vector<string> energies;
    vector<string> temp; // used to store returned value by 'split'
    bool get_lumo_next=false; // true only if when LUMO was the last column
    double done = false;
    while(!done && line.find("Root No.")!=string::npos) {
        temp = split(line);
        if(get_lumo_next) {
            getline(inf,line); // empty line
            getline(inf,line); // 1 A   2 A...
            getline(inf,line); // empty line
            getline(inf,line); // energies
            temp = split(line);
            for(vector<string>::iterator iter=temp.begin(); iter!=temp.end(); ++iter)
                energies.push_back(*iter);
            done = true;
        } else {
            vector<int> tmp_levels;
            for(vector<string>::size_type i=2; i<temp.size(); ++i)
                tmp_levels.push_back(atoi(temp[i].c_str()));
            vector<int>::iterator i_homo, i_lumo;
            i_homo = find(tmp_levels.begin(),tmp_levels.end(),homo_level);
            i_lumo = find(tmp_levels.begin(),tmp_levels.end(),homo_level+1);
            int homo_index=-1, lumo_index=-1;
            if(i_homo != tmp_levels.end())
                homo_index = i_homo-tmp_levels.begin();
            if(i_lumo != tmp_levels.end()) {
                lumo_index = i_lumo-tmp_levels.begin();
                // LUMO is the last column of the current block
                if(lumo_index == (int)(tmp_levels.size())-1)
                    get_lumo_next = true;
                else if(lumo_index < (int)(tmp_levels.size())-1)
                    done = true;
            }
            // to get energies
            for(int i=0; i<4; ++i)
                getline(inf,line);
            temp = split(line); // energies
            for(vector<string>::iterator iter=temp.begin(); iter!=temp.end(); ++iter)
                energies.push_back(*iter);
            // neither HOMO nor LUMO was in the current block
            if(homo_index==-1 && lumo_index==-1) {
                getline(inf,line);  // empty line
                while(getline(inf,line) && line!="")
                    ;
            }
            // to get HOMO or/and LUMO coefficients
            else {
                getline(inf,line);  // empty line
                getline(inf,line);
                int prev_id = 0;
                while(line != "") {
                    int current_id = atoi(line.substr(7,4).c_str());
                    if(current_id != prev_id) {
                        if(homo_index != -1)
                            homo_coeff.push_back(vector<double>());
                        if(lumo_index != -1)
                            lumo_coeff.push_back(vector<double>());
                        prev_id = current_id;
                    }
                    vector<string> coeff = split(line.substr(11));
                    if(homo_index != -1)
                        homo_coeff.back().push_back(atof(coeff[homo_index].c_str()));
                    if(lumo_index != -1)
                        lumo_coeff.back().push_back(atof(coeff[lumo_index].c_str()));
                    getline(inf,line);
                }
            }
            // go to the next block
            getline(inf,line); // empty line
            getline(inf,line); // supposet to be "Root No."
        }
    } // END of while-loop
    
    if(!done) {
        cout << endl << "<qc.cpp::read_fmo_coeff> Error: haven't collect all coefficients of HOMO and LUMO" << endl;
        return false;
    }
    // to fill fmo_energy
    for(int i=0; i<4; ++i)
        fmo_energy[i] = atof(energies[homo_level-2+i].c_str());

    // if has already read the line containing "NET ATOMIC CHARGE", then set the position to the begining!
    // otherwise, "netAtomicChargeMulliken" will fail.
    if(line.find("NET ATOMIC CHARGES") != string::npos)
        inf.seekg(0);
    return true;
}

// static method called by 'calcNA', 'calcEA'
// sum_\{i\}(Xi^2)
static double myfunc_sum(double x, double y)
{
    return x+pow(y,2);
}

vector<double> calcNA(vector<vector<double> > &homo_coeff, double energy)
{
    vector<double> result;
    for(vector<vector<double> >::iterator i=homo_coeff.begin(); i!=homo_coeff.end(); ++i) {
        double value = accumulate(i->begin(),i->end(),0.,myfunc_sum);
        value /= (1.-energy);
        result.push_back(value);
    }
    return result;
}

vector<double> calcEA(vector<vector<double> > &lumo_coeff, double energy)
{
    vector<double> result;
    for(vector<vector<double> >::iterator i=lumo_coeff.begin(); i!=lumo_coeff.end(); ++i) {
        double value = accumulate(i->begin(),i->end(),0.,myfunc_sum);
        value /= (10.+energy);
        result.push_back(value);
    }
    return result;
}

vector<double> calcRA(vector<vector<double> > &homo_coeff, vector<vector<double> > &lumo_coeff,
        double homo_energy, double lumo_energy)
{
    vector<double> result;
    for(vector<vector<double> >::size_type i=0; i<homo_coeff.size(); ++i) {
        double value = 0.;
        for(vector<double>::size_type j=0; j<homo_coeff[i].size(); ++j)
            value += (homo_coeff[i][j]*lumo_coeff[i][j]);
        value /= (lumo_energy - homo_energy);
        result.push_back(value);
    }
    return result;
}

bool netAtomicChargeMulliken(ifstream &inf, int natoms, vector<double> &result)
{
    result.clear();
    string line;
    while(getline(inf,line) && line.find("NET ATOMIC CHARGES")==string::npos)
        ;
    if(inf.eof()) {
        cout << endl << "<qc.cpp::netAtomicChargeMulliken>Error: can't find \"NET ATOMIC CHARGES\""
            << endl;
        return false;
    }
    vector<string> temp;
    result.clear();
    getline(inf,line); // empty line
    getline(inf,line); // ATOM NO. TYPE CHARGE ...
    for(int i=0; i<natoms; ++i) {
        getline(inf,line);
        temp = split(line);
        result.push_back(atof(temp[2].c_str()));
    }
    // next line: dipole !!!!!!!!!
    return true;
}

bool getElectronPopulation(ifstream &inf, vector<double> &result)
{
    result.clear();
    string line;
    vector<string> temp;
    while(getline(inf,line) && line.find("POPULATIONS")==string::npos)
        ;
    if(inf.eof()) {
        cout << endl << "<qc.cpp::calcPopulation>Error: can't find \"ATOMIC ORBITAL ELECTRON POPULATIONS\"" << endl;
        return false;
    }
    getline(inf,line); // empty line
    while(getline(inf,line) && line!="") {
        temp = split(line);
        for(vector<string>::iterator iter=temp.begin(); iter!=temp.end(); ++iter)
            result.push_back(atof(iter->c_str()));
    }

    return true;
}

static int orbitalSymbol2Index(const string &s1, const string &s2)
{
    int i,j;
    if(s1 == "S-SIGMA")
        i = 0;
    else if(s1 == "P-SIGMA")
        i = 1;
    else if(s1 == "P-PI")
        i = 2;
    else
        ; // never happen?!
    if(s2 == "S-SIGMA")
        j = 0;
    else if(s2 == "P-SIGMA")
        j = 1;
    else if(s2 == "P-PI")
        j = 2;
    else
        ; // never happen ?!
    return i*3+j;
}
bool sigmaPiBond(ifstream &inf, vector<string> &atom_symbol, double result[3][2])
{
    int natoms = static_cast<int>(atom_symbol.size());
    // initialize SIGMA-PI bond-order matrix
    // lower trianglar matrix
    vector<vector<vector<double> > > b_matrix;
    for(int i=0; i<natoms; ++i)
        b_matrix.push_back(vector<vector<double> >(i+1, vector<double>(9,0.)));
    string line;
    while(getline(inf,line) && line.find("SIGMA-PI BOND-ORDER MATRIX")==string::npos)
        ;
    if(inf.eof()) {
        cout << endl << "<qc.cpp::sigmaPiBond>Error: can't find \"SIGMA-PI BOND-ORDER MATRIX\"" << endl;
        return false;
    }
    getline(inf,line);  // empty line
    getline(inf,line);
    bool incomplete = false;
    int last_atom = natoms;
    int atoms_read = 0; // atoms appear in the column
    vector<int> y_atoms;
    vector<string> temp;
    while(line.find("SIGMA")!=string::npos || line.find("PI")!=string::npos) {
        // S-SIGMA    S-SIGMA    S-SIGMA    P-SIGMA      P-PI     S-SIGMA
        vector<string> y_types = split(line);
        // get atom ids in the column
        getline(inf,line);
        temp = split(line);  // ["C","1","C","1",...]
        if(!incomplete && last_atom==natoms) {
            y_atoms.clear();
            for(unsigned mm=0; mm<y_types.size(); ++mm) {
                if(y_types[mm] == "S-SIGMA")
                    ++atoms_read;
                y_atoms.push_back(atoms_read);
            }
        }
        /*
        cout << "y_atoms: ";
        copy(y_atoms.begin(), y_atoms.end(), std::ostream_iterator<int>(cout, " "));
        cout << endl;*/
        getline(inf,line); // ---------
        getline(inf,line);
        while(line != "") {
            if(line.size() == 1) {
                getline(inf,line);
                continue;
            }
            int current_atom = atoi(line.substr(11,4).c_str());
            last_atom = current_atom;
            vector<string>::size_type begin_i = (current_atom<=99)?3:2;
            string x_type = strip(line.substr(0,8));
            temp = split(line);
            for(vector<string>::size_type i=begin_i; i<temp.size(); ++i) {
                int index = orbitalSymbol2Index(x_type, y_types[i-begin_i]);
                b_matrix[current_atom-1][y_atoms[i-begin_i]-1][index] = atof(temp[i].c_str());
            }
            //appears to be end of the current block, but it's incomplete!!!!!
            if(current_atom==natoms && atom_symbol[current_atom-1]!="H" && x_type!="P-PI")
                incomplete = true;
            else
                incomplete = false;
            getline(inf,line);
        }
        getline(inf,line); // next block
    }
    // calculate SIGMA-SIGMA, SIGMA-PI, PI-PI
    result[0][0] = -100.; result[0][1] = 100.;
    result[1][0] = -100.; result[1][1] = 100.;
    result[2][0] = -100.; result[2][1] = 100.;
    double temp_value;
    for(int i=0; i<natoms; ++i) {
        for(int j=0; j<=i; ++j) {
            // SIGMA-SIGMA
            temp_value = b_matrix[i][j][0] + b_matrix[i][j][1] + b_matrix[i][j][3] + b_matrix[i][j][4];
            if(temp_value > result[0][0])
                result[0][0] = temp_value;
            if(temp_value < result[0][1])
                result[0][1] = temp_value;
            // SIGMA-PI
            temp_value = b_matrix[i][j][2] + b_matrix[i][j][5] + b_matrix[i][j][6] + b_matrix[i][j][7];
            if(temp_value > result[1][0])
                result[1][0] = temp_value;
            if(temp_value < result[1][1])
                result[1][1] = temp_value;
            // PI-PI
            temp_value = b_matrix[i][j][8];
            if(temp_value > result[2][0])
                result[2][0] = temp_value;
            if(temp_value < result[2][1])
                result[2][1] = temp_value;
        }
    }

    return true;
}

// 2014.02.13
// There is a bug when extracting SIGMA-PI BOND-ORDER MATRIX from file like
//     S-SIGMA    S-SIGMA    S-SIGMA    P-SIGMA      P-PI     S-SIGMA
//      H  1       H  2       N  3       N  3       N  3       C  4
//     ....
//     P-SIGMA      P-PI     S-SIGMA    P-SIGMA      P-PI 
//      C  4       C  4       N  5       N  5       N  5
//      ...
bool sigmaPiBond_back(ifstream &inf, vector<string> &atom_symbol, double result[3])
{
    int natom = static_cast<int>(atom_symbol.size());
    string line;
    while(getline(inf,line) && line.find("SIGMA-PI BOND-ORDER MATRIX")==string::npos)
        ;
    if(inf.eof()) {
        cout << endl << "<qc.cpp::sigmaPiBond>Error: can't find \"SIGMA-PI BOND-ORDER MATRIX\"" << endl;
        return false;
    }
    getline(inf,line);  // empty line
    getline(inf,line);
    vector<double> sigsig, sigpi, pipi;
    int last_atom = natom;
    vector<string> temp; // used to store returned value by 'split'
    int atoms_read = 0;
    vector<int> uniqueIDcol;  // atom IDs in the column (unique IDs)
    while(line.find("SIGMA")!=string::npos || line.find("PI")!=string::npos) {
        // 1. to get the two lines & atom ids in the column
        //     S-SGIMA P-SIGMA P-PI S-SIGMA P-SIGMA P-PI
        //     C  1     C  1   C  1 C   2   C   2   C  2
        vector<string> y_types = split(line);
        getline(inf,line);   // atoms  "C  1  C  1..."
        temp = split(line);  // ["C","1","C","1",...]
        //std::cout << "last_atom=" << last_atom << endl;
        // row is different from the row of last block
        if(last_atom == natom) {
            vector<string> y_atoms;
            unsigned nn=0;
            for(unsigned mm=0; mm<y_types.size(); ++mm) {
                if(y_types[mm] == "S-SIGMA")
                    ++atoms_read;
                y_atoms.push_back(atom_symbol[atoms_read-1]);
                if(atoms_read >= 100)
                    y_atoms.push_back(temp[nn].substr(atom_symbol[atoms_read-1].size()));
                else {
                    ++nn;
                    y_atoms.push_back(temp[nn]);
                }
                ++nn;
            }
            /*std::cout << endl << "y_atoms: ";
            copy(y_atoms.begin(),y_atoms.end(),std::ostream_iterator<string>(std::cout," "));
            std::cout << endl;*/

            uniqueIDcol.clear();  // atom IDs in the column (unique IDs)
            int prev_atom = 0;
            for(unsigned j=1; j<y_atoms.size(); j+=2) {
                int tmp_id = atoi(y_atoms[j].c_str());
                if(tmp_id != prev_atom) {
                    prev_atom = tmp_id;
                    uniqueIDcol.push_back(tmp_id);
                }
            }
        }
        /*
        std::cout << "uniqueIDcol: ";
        copy(uniqueIDcol.begin(), uniqueIDcol.end(), std::ostream_iterator<int>(std::cout," "));
        std::cout << endl;*/
        // begin to read orders...
        getline(inf,line);  // ------
        getline(inf,line);
        int current_atom = natom;
        while(line != "") {
            if(line.size() == 1) {
                getline(inf,line);
                continue;
            }
            current_atom = atoi(line.substr(11,4).c_str());
            vector<string>::size_type begin_i = (current_atom<=99)?3:2;
            // here we suppose that all heavy atoms have and only have "S-SIGMA, P-SIGMA, P-PI",
            // while, H has only "S-SGIMA"
            if(atom_symbol[current_atom-1] == "H") {
                temp = split(line);
                vector<double> tmp_order;
                for(vector<string>::size_type i=begin_i; i<temp.size(); ++i)
                    tmp_order.push_back(atof(temp[i].c_str()));
                int i=0;
                for(vector<int>::iterator iter=uniqueIDcol.begin(); 
                        iter!=uniqueIDcol.end() && *iter<=current_atom; ++iter) {
                    if(atom_symbol[*iter-1] == "H") { // sigsig
                        sigsig.push_back(tmp_order[i++]);
                        //std::cout << current_atom << " " << *iter << ", SIGMA-SIGMA=" << sigsig.back() << endl;
                    }
                    else {
                        sigsig.push_back(tmp_order[i]+tmp_order[i+1]); // sigsig
                        sigpi.push_back(tmp_order[i+2]); // sigpi
                        //std::cout << current_atom << " " << *iter << ", SIGMA-SIGMA=" << sigsig.back()
                        //    << ", SIGMA-PI=" << sigpi.back() << endl;
                        i += 3;
                    }
                }
                getline(inf,line);
            } else {
                vector<vector<double> > b_matrix;
                for(int ii=0; ii<3; ++ii) {
                    temp = split(line);
                    b_matrix.push_back(vector<double>());
                    for(vector<string>::size_type i=begin_i; i<temp.size(); ++i)
                        b_matrix.back().push_back(atof(temp[i].c_str()));
                    if(ii == 2)
                        break;
                    getline(inf,line);
                    temp = split(line);
                    if(temp.size() == 1)
                        getline(inf,line);
                    if(line == "") {
                        for(int j=0; j<4; ++j)
                            getline(inf,line);
                    }
                }
                /*
                cout << "bond_matrix:" << endl;
                for(vector<vector<double> >::iterator iter1=b_matrix.begin(); iter1!=b_matrix.end(); ++iter1) {
                    for(vector<double>::iterator iter2=iter1->begin(); iter2!=iter1->end(); ++iter2)
                        cout << " " << *iter2;
                    cout << endl;
                }
                */
                unsigned i=0;
                for(vector<int>::iterator iter=uniqueIDcol.begin();
                        iter!=uniqueIDcol.end() && *iter<=current_atom; ++iter) {
                    if(atom_symbol[*iter-1] == "H") { // sigsig
                        sigsig.push_back(b_matrix[0][i]+b_matrix[1][i]);
                        sigpi.push_back(b_matrix[2][i]);
                        //std::cout << current_atom << " " << *iter << ", SIGMA-SIGMA=" << sigsig.back() << ", SIGMPA-PI=" << sigpi.back() << endl;
                        i += 1;
                    }
                    else {
                        if(current_atom == *iter) {
                            sigsig.push_back(b_matrix[0][i]+b_matrix[1][i+1]);
                            sigpi.push_back(0.);
                            pipi.push_back(b_matrix[2][i+2]);
                            //std::cout << current_atom << " " << *iter << ", SIGMA-SIGMA=" << sigsig.back()
                            //    << ", SIGMA-PI=" << sigpi.back() << ", PI-PI=" << pipi.back() << endl;
                        } else {
                            sigsig.push_back(b_matrix[0][i]+b_matrix[0][i+1]+
                                    b_matrix[1][i]+b_matrix[1][i+1]);
                            sigpi.push_back(b_matrix[0][i+2]+b_matrix[1][i+2]+
                                    b_matrix[2][i]+b_matrix[2][i+1]);
                            pipi.push_back(b_matrix[2][i+2]);
                            //std::cout << current_atom << " " << *iter << ", SIGMA-SIGMA=" << sigsig.back()
                            //    << ", SIGMA-PI=" << sigpi.back() << ", PI-PI=" << pipi.back() << endl;
                        }
                        i += 3;
                    }
                }
                // next block of bond order matrix
                getline(inf,line);
            } // END of if-else
        } // END of while-loop
        last_atom = current_atom;
        // go to next block
        getline(inf,line);
    }

    result[0] = *max_element(sigsig.begin(),sigsig.end());
    result[1] = *max_element(sigpi.begin(),sigpi.end());
    result[2] = *max_element(pipi.begin(),pipi.end());
    return true;
}

// called by 'energyPartitionTerms'
static double getEnergy(string &line)
{
    int last2, last;
    last2 = 0;
    last = 0;
    for(string::size_type i=1; i<line.size(); ++i) {
        if(line[i] == ' ') {
            last2 = last;
            last = i;
        }
    }
    return atof(line.substr(last2,last-last2-1).c_str());
}
bool energyPartitionTerms(ifstream &inf, int natoms, double result[7][2])
{
    int i;
    string line;
    while(getline(inf,line) && line.find("SUMMARY OF ENERGY PARTITION")==string::npos)
        ;
    if(inf.eof()) {
        cout << "Error<qc.cpp::energyPartitionTerms>: can't find \"SUMMARY OF ENERGY PARTITION\"" << endl;
        return false;
    }
    for(i=0; i<3; ++i)
        getline(inf,line);
    getline(inf,line);  // ELECTRON-NUCLEAR
    result[0][0] = getEnergy(line);
    getline(inf,line);  // ELECTRON-ELECTRON
    result[1][0] = getEnergy(line);
    for(i=0; i<5; ++i)
        getline(inf,line);
    getline(inf,line);  // RESONANCE
    result[2][0] = getEnergy(line);
    getline(inf,line);  // EXCHANGE
    result[3][0] = getEnergy(line);
    for(i=0; i<3; ++i)
        getline(inf,line);
    for(i=0; i<3; ++i) {  // E-E, E_N, N-N (two-center)
        getline(inf,line);
        result[4+i][0] = getEnergy(line);
    }
    // result[][1]
    for(i=0; i<7; ++i)
        result[i][1] = result[i][0] / natoms;
    return true;
}

bool netAtomicChargeESP(ifstream &inf, int natoms, vector<double> &result)
{
    result.clear();
    string line;
    while(getline(inf,line) && line.find("ELECTROSTATIC POTENTIAL CHARGES")==string::npos)
        ;
    if(inf.eof()) {
        cout << endl << "Error<qc.cpp::netAtomicChargeESP>: can't find \"ELECTROSTATIC POTENTIAL CHARGES\"" << endl;
        return false;
    }
    getline(inf,line);  //empty line
    getline(inf,line);
    for(int i=0; i<natoms; ++i) {
        getline(inf,line);
        vector<string> temp = split(line);
        result.push_back(atof(temp[2].c_str()));
    }
    // now it's empty line
    return true;
}

// static method called by 'qcfeatures'
// result:
//   [i][]: each for C, N, and O
//   [][j]: each for min, max, and average
static void transform_reactivity(vector<double> &orig, vector<string> &atom_symbol, double result[3][3])
{
    vector<vector<double> > temp(3);  // each for C, N, and O
    for(vector<double>::size_type i=0; i<orig.size(); ++i) {
        if(atom_symbol[i] == "C")
            temp[0].push_back(orig[i]);
        else if(atom_symbol[i] == "N")
            temp[1].push_back(orig[i]);
        else if(atom_symbol[i] == "O")
            temp[2].push_back(orig[i]);
        else
            ;
    }
    for(int i=0; i<3; ++i) {
        if(temp[i].size() == 0) {
            result[i][0] = 1E+38;
            result[i][1] = 1E+38;
            result[i][2] = 1E+38;
        } else {
            result[i][0] = *min_element(temp[i].begin(),temp[i].end());
            result[i][1] = *max_element(temp[i].begin(),temp[i].end());
            result[i][2] = accumulate(temp[i].begin(),temp[i].end(),0.) / temp[i].size();
        }
    }
}
// static method called by 'qcfeatures'
// result[i][]: each for C, N, O, H, and all
// result[][j]: each for max and min
static void transform_charge(vector<double> &charges, vector<string> &atom_symbol, double result[5][2])
{
    for(int i=0; i<5; ++i) {
        result[i][0] = -1E+38;
        result[i][1] = 1E+38;
    }
    for(vector<string>::size_type i=0; i<atom_symbol.size(); ++i) {
        int j;
        if(atom_symbol[i] == "C")
            j = 0;
        else if(atom_symbol[i] == "N")
            j = 1;
        else if(atom_symbol[i] == "O")
            j = 2;
        else if(atom_symbol[i] == "H")
            j = 3;
        else
            j = -1;
        // max. and min. net atomic charge for C, N, O, or H atom
        if(j != -1) {
            if(charges[i] > result[j][0])
                result[j][0] = charges[i];
            if(charges[i] < result[j][1])
                result[j][1] = charges[i];
        }
        // max. and min. net atomic charge
        if(charges[i] > result[4][0])
            result[4][0] = charges[i];
        if(charges[i] < result[4][1])
            result[4][1] = charges[i];
    }
}

// to calculate QC features
bool qcfeatures(string mopout_file, QCValues &values,
        vector<double> &mulliken, vector<double> &esp)
{
    mulliken.clear();
    esp.clear();
    ifstream inf(mopout_file.c_str());
    if(!inf) {
        cout << endl << "<qc.cpp::qcfeatures>Error: can't open file " << mopout_file << endl;
        return false;
    }
    vector<string> temp;  // used to store return values of 'split'
    string line;
    // 1. get atom symbols
    while(getline(inf,line) && line.find("CARTESIAN COORDINATES")==string::npos)
        ;
    if(inf.eof()) {
        cout << endl << "<qc.cpp::qcfeatures>Error: can't find \"CARTESIAN COORDINATES\" in file "
            << mopout_file << endl;
        inf.close();
        return false;
    }
    getline(inf,line);  // empty line
    getline(inf,line);  // "  NO.  ATOM  X  Y  Z"
    getline(inf,line);  // empty line
    vector<string> atom_symbol;
    while(getline(inf,line) && line.size() != 0) {
        temp = split(line);
        atom_symbol.push_back(temp[1]);
    }
    // 2. get final heat of formation
    while(getline(inf,line) && line.find("FINAL HEAT OF FORMATION")==string::npos)
        ;
    if(inf.eof()) {
        cout << endl << "<qc.cpp::qcfeatures>Error: can't find \"FINAL HEAT OF FORMATION\" in file "
            << mopout_file << endl;
        inf.close();
        return false;
    }
    temp = split(line);
    values.heat_formation[0] = atof(temp[temp.size()-2].c_str());
    values.heat_formation[1] = values.heat_formation[0] / atom_symbol.size();
    // 3. get No. of filled levels
    while(getline(inf,line) && line.find("NO. OF FILLED LEVELS")==string::npos)
        ;
    if(inf.eof()) {
        cout << endl << "<qc.cpp::qcfeatures>Error: can't find \"NO. OF FILLED LEVELS\" in file "
            << mopout_file << endl;
        inf.close();
        return false;
    }
    temp = split(line);
    int homo_level = atoi(temp.back().c_str());
    getline(inf,line);
    if(line.find("OPEN LEVELS") != string::npos)
        homo_level += atoi(line.substr(line.find_last_of(" ")+1).c_str());
    values.occ_levels[0] = (double)homo_level;
    values.occ_levels[1] = values.occ_levels[0] / atom_symbol.size();
    // 4. to get front orbitals
    vector<vector<double> > homo_coeff, lumo_coeff;
    if(!read_fmo_coeff(inf,homo_level,values.fmo_energy,homo_coeff,lumo_coeff)) {
        inf.close();
        return false;
    }
    // 5. calculate nucleophilic reactivity index
    vector<double> na_result = calcNA(homo_coeff, values.fmo_energy[1]);
    transform_reactivity(na_result,atom_symbol,values.na);
    //calcReactivity(homo_coeff, atom_symbol, values.fmo_energy[1], values.na, funcNA);
    // 6. calculate electrophilic reactivity index
    vector<double> ea_result = calcEA(lumo_coeff, values.fmo_energy[2]);
    transform_reactivity(ea_result,atom_symbol,values.ea);
    //calcReactivity(lumo_coeff, atom_symbol, values.fmo_energy[2], values.ea, funcEA);
    // 7. calculate Fukui 1-electron reactivity index
    vector<double> ra_result = calcRA(homo_coeff, lumo_coeff, values.fmo_energy[1], values.fmo_energy[2]);
    transform_reactivity(ra_result,atom_symbol,values.ra);
    //calcRA(homo_coeff, lumo_coeff, atom_symbol, values.fmo_energy[1], values.fmo_energy[2], values.ra);
    // 8. net atomic charges
    //  - Mulliken
    if(!netAtomicChargeMulliken(inf, atom_symbol.size(),mulliken)) {
        inf.close();
        return false;
    }
    transform_charge(mulliken,atom_symbol,values.chargeMulliken);

    // 9. dipole moment
    getline(inf,line);
    for(int i=0; i<3; ++i) {
        getline(inf,line);
        temp = split(line);
        values.dipole[i] = atof(temp[4].c_str());
    }
    
    // 10. min/max atomic orbital electronic population
    vector<double> population;
    if(!getElectronPopulation(inf, population)) {
        inf.close();
        return false;
    }
    values.populationAO[0] = *min_element(population.begin(), population.end());
    values.populationAO[1] = *max_element(population.begin(), population.end());

    // 11. SIGMA/PI bond order
    //cout << "to sigmaPiBond" << endl;
    if(!sigmaPiBond(inf,atom_symbol,values.sigpiBond)) {
        inf.close();
        return false;
    }
    //cout << "after sigmaPiBond" << endl;

    // 12. energy partition terms in AM1  (one-center & two-center)
    if(!energyPartitionTerms(inf, atom_symbol.size(), values.energyTerms)) {
        inf.close();
        return false;
    }
    //cout << "after energyPartitionTerms" << endl;

    // 13. electrostatic potential charge
    if(!netAtomicChargeESP(inf, atom_symbol.size(), esp)) {
        inf.close();
        return false;
    }
    //cout << "to transform_charge" << endl;
    transform_charge(esp,atom_symbol,values.chargeESP);
    //cout << "after netAtomicChargeESP" << endl;

    // more features...
    // ...
    inf.close();

    return true;
}

static inline void printReactivity(ostream &os, const double value[3][3], string &sep)
{
    for(int i=0; i<3; ++i)
        for(int j=0; j<3; ++j)
            os << sep << value[i][j];
}
void print_qc_result(ostream &os, const QCValues &result, string sep)
{
    int i;
    os << result.heat_formation[0] << sep << result.heat_formation[1]
        << sep << result.occ_levels[0] << sep << result.occ_levels[1];
    for(i=0; i<4; ++i)
        os << sep << result.fmo_energy[i];
    printReactivity(os,result.na,sep);
    printReactivity(os,result.ea,sep);
    printReactivity(os,result.ra,sep);
    for(i=0; i<5; ++i)
        os << sep << result.chargeMulliken[i][0] << sep << result.chargeMulliken[i][1];
    os << sep << result.dipole[0] << sep << result.dipole[1] << sep << result.dipole[2]
        << sep << result.populationAO[0] << sep << result.populationAO[1]
        << sep << result.sigpiBond[0][0] << sep << result.sigpiBond[0][1]
        << sep << result.sigpiBond[1][0] << sep << result.sigpiBond[1][1]
        << sep << result.sigpiBond[2][0] << sep << result.sigpiBond[2][1];
    for(i=0; i<7; ++i)
        os << sep << result.energyTerms[i][0] << sep << result.energyTerms[i][1];
    for(i=0; i<5; ++i)
        os << sep << result.chargeESP[i][0] << sep << result.chargeESP[i][1];
}


