/*=============================================================================
#     FileName: main.cpp
#         Desc: 
#       Author: jlpeng
#        Email: jlpeng1201@gmail.com
#     HomePage: 
#      Version: 0.0.1
#      Created: 2012-12-29 09:18:38
#   LastChange: 2014-03-20 17:05:14
#      History:
=============================================================================*/

#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <fstream>
#include <utility>
#include <iterator>
#include <numeric>
#include <cstdlib>
#include <ctime>
#include "config.h"
#include "constitutional.h"
#include "geometrical.h"
#include "topological.h"
#include "qc.h"
#include "electrostatic.h"
#include "tools.h"
#include "physicochemical.h"
#include "names.h"
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>

using std::cout;
using std::cerr;
using std::endl;
using std::map;
using std::string;
using std::vector;
using std::ostream;
using std::ifstream;
using std::ofstream;
using std::pair;
using std::ostream_iterator;
using std::accumulate;
using OpenBabel::OBMol;
using OpenBabel::OBConversion;

bool calcFeatures(ofstream &outf, string &infile, int type, int density, bool codessa);
inline void print_names(ostream &os, vector<string> &names, string sep)
{
    for(vector<string>::const_iterator iter=names.begin(); iter!=names.end(); ++iter)
        os << sep << *iter;
}
inline void print_names_cpsa(ostream &os, string suffix, string sep)
{
    for(vector<string>::iterator iter=cpsa_names.begin(); iter!=cpsa_names.end(); ++iter)
        os << sep << *iter << suffix;
}

void exit_with_help(const char *name)
{
    //time_t rawtime;
    //time(&rawtime);
    cerr << endl
        << "To calculate descriptors for given molecule(s)" << endl
        << "                                                   www.dddc.ac.cn/admetus" << endl
        << "=========================================================================" << endl
        << "Usage: calcDescriptors.exe [options]" << endl << endl
        << "[options]" << endl
        << "Input:" << endl
        << "  -m molfile" << endl
        << "     support format: mol2, mol, sdf, mopout" << endl
        << "  -l mollist" << endl
        << "     each line should be a molecule, e.g. F:\\mols\\1.mol" << endl
        << endl << "Output:" << endl
        << "  -o output" << endl
        << "     specify the file to save the calculated descriptors" << endl
        << endl << "Descriptors:" << endl
        << "  -a: calculate all available descriptors" << endl
        << "  -c: constitutional descriptors" << endl
        << "  -e: electrostatic descriptors" << endl
        << "  -g: geometrical descriptors" << endl
        << "  -p: physico-chemical descriptors" << endl
        << "  -q: quantum-chemical descriptors" << endl
        << "      'mopout file' should be passed via '-m' or '-l'!" << endl
        << "  -t: topological descriptors" << endl
        << endl << "Additional:" << endl
        << "  -d density: <default: 10>" << endl
        << "     density used for some geometrical descriptors, like shadow and volume" << endl
        //<< "  --codessa" << endl
        //<< "    Only valid when '-g' is passed" << endl
        //<< "    shadow and volume will be rescaled only when density equals to 5" << endl
        //<< "    (relationship was learned using 903 compounds)" << endl
        << endl << "Information: (descriptors will not be calculated)" << endl
        << "  --names" << endl
        << "    to check available descriptors;" << endl
        << "    it can be used in combination with '-c,-g,-t,-q,-e' to show names of" << endl
        << "    specific kind of descriptors, or used alone to show all names." << endl
        << "  --ref" << endl
        << "    display reference" << endl
        << "  --version" << endl
        << "    show version" << endl
        << "  --help" << endl
        << "    display this information" << endl << endl
        << "Attention:" << endl
        << "  1. in case that value of a descriptor can't be calculated," << endl
        << "     it'll be given a value of +/-1E38" << endl
        << "=========================================================================" << endl
        << "                                                     jlpeng1201@gmail.com" << endl << endl;
        //<< "                                                 " << ctime(&rawtime) << endl;
    exit(EXIT_FAILURE);
}

void exit_with_help_simple(const char *name)
{
    cerr << endl
        << "    OBJ: To calculate descriptors for given molecule(s)" << endl << endl
        << "  Usage: calcDescriptors.exe [-ml input] [-o output] [-acegpqt]" << endl << endl
        << "  Try --help for more detailed information" << endl << endl;
    exit(EXIT_FAILURE);
}

int main(int argc, char *argv[])
{
    if(argc < 2)
        exit_with_help_simple(argv[0]);

    // get root dir
    char *home_path = getenv("DESCRIPTOR_HOME");
    string root_dir = "";
    if(home_path) {
        root_dir = home_path;
        if(root_dir.rfind(SEP) != root_dir.size()-1)
            root_dir += SEP;
    }
    else {
        string prog_name(argv[0]);
        string::size_type _i = prog_name.rfind(SEP);
        if(_i != string::npos)
            root_dir = prog_name.substr(0,_i+1);
    }

    // parse option
    string molfile="", mollist="", output="";
    int density = 10;
    bool codessa = false, names=false, help=false, all=false, ref=false, version=false;
    int type = 0;  // constitute, geometrical, topological, quantum-chemical, electrostatic
    for(int i=1; i<argc; ++i) {
        if(strlen(argv[i]) == 2) {
            if(argv[i][0] != '-') {
                cerr << endl << "Error: illegal option: " << argv[i] << endl;
                exit(EXIT_FAILURE);
            }
            if(argv[i][1] == 'm')
                molfile = argv[++i];
            else if(argv[i][1] == 'l')
                mollist = argv[++i];
            else if(argv[i][1] == 'o')
                output = argv[++i];
            else if(argv[i][1] == 'a')
                all = true;
            else if(argv[i][1] == 'c')
                type |= 1;
            else if(argv[i][1] == 'e')
                type |= (1<<1);
            else if(argv[i][1] == 'g')
                type |= (1<<2);
            else if(argv[i][1] == 'd')
                density = atoi(argv[++i]);
            else if(argv[i][1] == 'p')
                type |= (1<<3);
            else if(argv[i][1] == 'q')
                type |= (1<<4);
            else if(argv[i][1] == 't')
                type |= (1<<5);
            else {
                cerr << endl << "Error: illegal option: " << argv[i] << endl;
                exit(EXIT_FAILURE);
            }
        } else {
            if(strcmp(argv[i],"--codessa")==0)
                codessa = true;
            else if(strcmp(argv[i],"--names")==0)
                names = true;
            else if(strcmp(argv[i],"--ref")==0)
                ref = true;
            else if(strcmp(argv[i],"--version")==0)
                version = true;
            else if(strcmp(argv[i],"--help")==0)
                help = true;
            else {
                cerr << endl << "Error: illegal option: " << argv[i] << endl;
                exit(EXIT_FAILURE);
            }
        }
    }

    // to display help information
    if(help)
        exit_with_help(argv[0]);

    // to display reference
    if(ref) {
        cout << reference << endl;
        return 0;
    }

    // to display version
    if(version) {
        cout << endl << VERSION << endl;
        return 0;
    }

    // for calculating all available descriptors
    if(all) {
        type = 0;
        type |= (1<<6);
        type -= 1;
    }

    // print descriptor names
    if(names) {
        unsigned cpsa_num = cpsa_names.size();
        unsigned count = 0;
        cout << "###names of available descriptors###";
        if(type==0 || (type&1)) {
            cout << endl << ">>> constitutional (" << constitute_names.size() << ")";
            print_names(cout, constitute_names, "\n");
            count += constitute_names.size();
        }
        if(type==0 || (type&(1<<1))) {
            cout << endl << ">>> electrostatic (" << electrostatic_names.size()+cpsa_num << ")";
            print_names(cout, electrostatic_names, "\n");
            print_names_cpsa(cout, " [Zefirov]","\n");
            count += (electrostatic_names.size()+cpsa_num);
        }
        if(type==0 || (type&(1<<2))) {
            cout << endl << ">>> geometrical (" << geometrical_names.size() << ")";
            print_names(cout, geometrical_names, "\n");
            count += geometrical_names.size();
        }
        if(type==0 || (type&(1<<3))) {
            cout << endl << ">>> physico-chemical (" << physicochemical_names.size() << ")";
            print_names(cout, physicochemical_names, "\n");
            count += physicochemical_names.size();
        }
        if(type==0 || (type&(1<<4))) {
            cout << endl << ">>> quantum-chemical (" << qc_names.size()+2*cpsa_num << ")";
            print_names(cout, qc_names, "\n");
            print_names_cpsa(cout, " [Mulliken]","\n");
            print_names_cpsa(cout, " [ESP]","\n");
            count += (qc_names.size()+2*cpsa_num);
        }
        if(type==0 || (type&(1<<5))) {
            cout << endl << ">>> topological (" << topological_names.size() << ")";
            print_names(cout, topological_names,"\n");
            count += topological_names.size();
        }
        // more features here...
        cout << endl 
            << "==============================================" << endl
            << "totally " << count << " descriptors" << endl << endl;
        return 0;
    }

    // check options
    if(codessa && density!=5)
        cerr << endl << "Warning: `--codessa` will be ignored when `-d density` is not 5!" << endl;
    if(molfile=="" && mollist=="") {
        cerr << endl << "Error: at least one of '-m -l' should be given!" << endl;
        exit(EXIT_FAILURE);
    }
    if(output == "") {
        cerr << endl << "Error: miss option '-o'" << endl;
        exit(EXIT_FAILURE);
    }
    if(type == 0) {
        cerr << endl << "Error: at least one of '-a -c -g -t -q' should be given!" << endl;
        exit(EXIT_FAILURE);
    }
    if(!(type&(1<<1)) && codessa)
        cerr << endl << "option '--codessa' not used" << endl;

    // open the output file
    ofstream outf(output.c_str());
    if(!outf) {
        cerr << endl << "Error: failed to open file " << output << endl;
        exit(EXIT_FAILURE);
    }
    // print simplified names of descriptors to be calculated
    outf << "mol";
    if(type&1) // constitute
        print_names(outf, constitute_names, "\t");
    if(type&(1<<1)) {
        print_names(outf, electrostatic_names, "\t");
        print_names_cpsa(outf, " [Zefirov]","\t");
    }
    if(type&(1<<2))
        print_names(outf, geometrical_names, "\t");
    if(type&(1<<3)) {
        print_names(outf, physicochemical_names, "\t");
    }
    if(type&(1<<4)) {
        print_names(outf, qc_names, "\t");
        print_names_cpsa(outf, " [Mulliken]","\t");
        print_names_cpsa(outf, " [ESP]", "\t");
    }
    if(type&(1<<5))
        print_names(outf, topological_names, "\t");
    // add more feature names here...
    outf << endl;

    // to calculate descriptors
    if(molfile != "")
        calcFeatures(outf, molfile,type,density,codessa);
    if(mollist != "") {
        ifstream inf(mollist.c_str());
        if(!inf) {
            cerr << endl << "Error: can't open file " << mollist << endl;
            exit(EXIT_FAILURE);
        }
        int num_fail = 0;
        string line;
        while(getline(inf,line)) {
            if(line.size() == 0)
                continue;
            if(line[0] == '#')
                continue;
            if(!calcFeatures(outf, line,type,density,codessa))
                ++num_fail;
        }
        inf.close();
        if(num_fail != 0)
            cout << "Failed to calculate descriptors for " << num_fail << " molecule(s)" << endl;
    }

    outf.close();

    return 0;
}

// function to calculate descriptors for each molecule
// infile: in case to calculate QC descriptors
static void predict_each(ofstream &outf, string &infile, OBMol &mol, int type, int density, bool codessa)
{
    if(type&1) {  // constitute
        //cout << "constitute" << endl;
        vector<double> result = atoms(mol);
        for(vector<double>::iterator iter=result.begin(); iter!=result.end(); ++iter)
            outf << "\t" << *iter;
        result.clear();
        result = bonds(mol);
        for(vector<double>::iterator iter=result.begin(); iter!=result.end(); ++iter)
            outf << "\t" << *iter;
        result.clear();
        result = rings(mol);
        for(vector<double>::iterator iter=result.begin(); iter!=result.end(); ++iter)
            outf << "\t" << *iter;
        result.clear();
        result = hbonds(mol);
        for(vector<double>::iterator iter=result.begin(); iter!=result.end(); ++iter)
            outf << "\t" << *iter;
    }
    if(type&(1<<1)) { // electrostatic
        //cout << "electrostatic" << endl;
        vector<double> zefirov_pc = calcZefirovPC(mol);
        double result[12][2];
        calcMaxMinZefirovPC(mol,zefirov_pc,result);
        for(int i=0; i<12; ++i)
            outf << "\t" << result[i][0] << "\t" << result[i][1];
        double result1[2];
        calcTopolElectronicIndex(mol,zefirov_pc,result1);
        outf << "\t" << result1[0] << "\t" << result1[1];
        double result2[38];
        calcCPSA(mol,zefirov_pc,result2);
        for(int i=0; i<38; ++i)
            outf << "\t" << result2[i];
    }
    if(type&(1<<2)) { // geometry
        //cout << "geometry" << endl;
        // - principle moment of inertia
        //   if codessa was set to be true, then transformation will be
        //   carried out in order to be consistent with those of CODESSA
        vector<pair<int,double> > moi = calcMOI(mol);
        int axes[3];
        for(vector<pair<int,double> >::size_type i=0; i<moi.size(); ++i) {
            axes[i] = moi[i].first;
            if(codessa) {
                if(moi[i].second == 0.)
                    outf << "\t1E+38";
                else
                    outf << "\t" << 16.86/moi[i].second;
            } else
                outf << "\t" << moi[i].second;
        }
        // - shadow area
        //   if codessa was set to be true and density equals to 5, then
        //   transformation will be carried out in order to be consistent
        //   with those of CODESSA.
        double shadow_result[6];
        shadow(mol,axes,density,shadow_result);
        if(codessa && density==5) {
            outf << "\t" << shadow_result[0]/1.0335 << "\t" << shadow_result[1]/1.0425 
                << "\t" << shadow_result[2]/1.0434 << "\t" << shadow_result[3]/1.0262
                << "\t" << shadow_result[4]/1.0351 << "\t" << shadow_result[5]/1.0358;
        } else {
            for(int i=0; i<6; ++i)
                outf << "\t" << shadow_result[i];
        }
        // - molecular volume
        //   if codessa was set to be true and density equals to 5, then
        //   transformation will be carried out in order to be consistent
        //   with those of CODESSA.
        double volume_result[2];
        volume(mol,axes,density,volume_result);
        if(codessa && density==5)
            outf << "\t" << volume_result[0]/1.1196 << "\t" << volume_result[1]/1.1019;
        else
            outf << "\t" << volume_result[0] << "\t" << volume_result[1];
        // - gravitation index
        double grav_result[2];
        gravitation(mol,grav_result);
        outf << "\t" << grav_result[0] << "\t" << grav_result[1];
        // - TPSA topological polar surface area
        double tpsa = calcTPSA(mol);
        outf << "\t" << tpsa;
    }
    if(type&(1<<3)) { // physicochemical
        //cout << "physicochemical" << endl;
        double logp = calcLogP(mol);
        outf << "\t" << logp;
    }
    if(type&(1<<4)) { // quantum-chemical
        //cout << "QC" << endl;
        QCValues qc_result;
        vector<double> mulliken, esp;
        bool success = qcfeatures(infile, qc_result, mulliken, esp);
        if(!success) {
            for(unsigned k=0; k<qc_names.size(); ++k)
                outf << "\t" << 1E+38;
            for(unsigned k=0; k<2*cpsa_names.size(); ++k)
                outf << "\t" << 1E+038;
        } else {
            outf << "\t";
            print_qc_result(outf, qc_result, "\t");
            double result[38];
            // CPSA[quantum-chemical charge]
            calcCPSA(mol,mulliken,result);
            for(int i=0; i<38; ++i)
                outf << "\t" << result[i];
            // CPSA[ESP]
            calcCPSA(mol,esp,result);
            for(int i=0; i<38; ++i)
                outf << "\t" << result[i];
        }
    }
    if(type&(1<<5)) { // topological
        //cout << "topological" << endl;
        // - Wiener index
        outf << "\t" << wiener(mol);
        double result[4];
        // - Randic index
        randic(mol,result);
        for(int i=0; i<4; ++i)
            outf << "\t" << result[i];
        // - Kier&Hall index
        kier_hall(mol,result);
        for(int i=0; i<4; ++i)
            outf << "\t" << result[i];
        // - Kier shape index
        kier_shape(mol,result);
        for(int i=0; i<4; ++i)
            outf << "\t" << result[i];
        // - Balaban index
        outf << "\t" << balaban(mol);
        // - infoContent
        double ic_result[3][8];
        infoContent(mol,ic_result);
        for(int i=0; i<3; ++i)
            for(int j=0; j<8; ++j)
                outf << "\t" << ic_result[i][j];
    }
    // calculate more descriptors here...
    outf << endl;
}

bool calcFeatures(ofstream &outf, string &infile, int type, int density, bool codessa)
{
    cout << ">>> " << infile << endl;
    string::size_type i = infile.rfind(".");
    if(i == string::npos) {
        cout << "Error: no file format given: " << infile << endl;
        return false;
    }
    string _format = infile.substr(i+1);
    OBConversion conv;
    OBMol mol;
    if(!conv.SetInFormat(_format.c_str())) {
        cout << "Error: illegal format " << _format << endl;
        return false;
    }
    bool success = conv.ReadFile(&mol,infile);
    int count = 0;
    while(success) {
        ++count;
        const char *title = mol.GetTitle();
        if(title[0] == '\0') {
            cout << "    title: " << infile;
            outf << infile;
            if(_format=="mol2" || _format=="sdf") {
                outf << "_" << count;
                cout << "_" << count;
            }
        }
        else {
            outf << title;
            cout << "    title: " << title;
        }
        cout << endl;
        predict_each(outf, infile, mol, type, density, codessa);
        success = conv.Read(&mol);
    }
    
    if(count == 0) {
        cout << "Error: failed to read molecule from " << infile << " OR there is no molecule at all!" << endl;
        return false;
    } else
        return true;
}


