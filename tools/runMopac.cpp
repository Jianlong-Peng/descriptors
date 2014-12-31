/*=============================================================================
#     FileName: runMopac.cpp
#         Desc: 
#       Author: jlpeng
#        Email: jlpeng1201@gmail.com
#     HomePage: 
#      Created: 2014-02-19 18:50:29
#   LastChange: 2014-03-20 17:10:19
#      History:
=============================================================================*/
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <cstring>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/builder.h>
#include <openbabel/forcefield.h>


using namespace std;
using namespace OpenBabel;

vector<string> runMopac(string &infile, string &mopac7exe, string &keywords);

void exit_with_help(char *name)
{
    cerr << endl << "Usage"<< endl
        << "  " << name << " [options]" << endl
        << endl << "[options]" << endl
        << "  -m inmol" << endl
        << "     all formats OpenBabel supported" << endl
        << "  -l inlist" << endl
        << "     file, each line should be a molecule, e.g. F:\\mols\\1.mol" << endl
        << "  --mopac dir" << endl
        << "     directory where mopac7.exe could be found" << endl
        << "     e.g. F:\\software" << endl
        << "     default: mopac7.exe could be found in the current dir or in PATH" << endl
        << "  --key words" << endl
        << "     additional keywords passed to mopac7. separated by ';'" << endl
        << "     e.g. DOUBLET;BOND" << endl
        << endl << "Attention:" << endl
        << "  1. output files will be named as 'basename_i.mopout'" << endl
        << "  2. output files will be generated in the same directory as its input" << endl
        << "  3. keywords being used by default:" << endl
        << "     ESP MMOK VECTORS BONDS PI PRECISE ENPART AM1 GNORM=0.01" << endl
        << "     EF XYZ MULLIK Charge" << endl
        << endl;
    exit(EXIT_FAILURE);
}

int main(int argc, char *argv[])
{
    if(argc<3 || argc>9)
        exit_with_help(argv[0]);

    string molfile("");
    string mollist("");
    string mopacdir("");
    string keywords("");
    string tmp_keywords("");
    for(int i=1; i<argc; ++i) {
        if(strcmp(argv[i], "-m") == 0)
            molfile = argv[++i];
        else if(strcmp(argv[i],"-l") == 0)
            mollist = argv[++i];
        else if(strcmp(argv[i],"--mopac") == 0)
            mopacdir = argv[++i];
        else if(strcmp(argv[i],"--key") == 0)
            tmp_keywords = argv[++i];
        else {
            cerr << "Error: invalid option " << argv[i] << endl;
            exit(EXIT_FAILURE);
        }
    }

    if(molfile=="" && mollist=="") {
        cerr << "Error: at least one of '-m -l' is needed!" << endl;
        exit(EXIT_FAILURE);
    }
    
    string mopac7exe;
    if(mopacdir == "")
        mopac7exe = "mopac7.exe";
    else {
        if(mopacdir[mopacdir.size()-1] == '\\')
            mopac7exe = mopacdir + "mopac7.exe";
        else
            mopac7exe = mopacdir + "\\mopac7.exe";
    }

    string::size_type i=0;
    while(i < tmp_keywords.size()) {
        int count = 0;
        int start_pos = i;
        while(i<tmp_keywords.size() && tmp_keywords[i]!=';') {
            ++count;
            ++i;
        }
        keywords += (" " + tmp_keywords.substr(start_pos, count));
    }

    vector<string> fail_list;

    if(molfile != "") {
        vector<string> tmp = runMopac(molfile, mopac7exe, keywords);
        for(vector<string>::iterator i=tmp.begin(); i!=tmp.end(); ++i)
            fail_list.push_back(*i);
    }
    if(mollist != "") {
        ifstream inf(mollist.c_str());
        if(!inf) {
            cerr << "Error: failed to open file " << mollist << endl;
            exit(EXIT_FAILURE);
        }
        string line;
        while(getline(inf,line)) {
            if(line.size() == 0)
                continue;
            if(line[0] == '#')
                continue;
            vector<string> tmp = runMopac(line, mopac7exe, keywords);
            for(vector<string>::iterator i=tmp.begin(); i!=tmp.end(); ++i)
                fail_list.push_back(*i);
        }
        inf.close();
    }

    if(!fail_list.empty()) {
        cout << endl << "failed to run mopac7 for the following files:" << endl;
        for(vector<string>::iterator i=fail_list.begin(); i!=fail_list.end(); ++i)
            cerr << "  " << *i << endl;
    } else
        cout << endl << "All finished!" << endl;

    return 0;
}

static void make3D(OBMol &mol)
{
    // reference: pybel.py::Molecule::make3D
    OBBuilder builder;
    builder.Build(mol);
#if defined(_MSC_VER)
    OBForceField* ff = (OBForceField *)OBPlugin::GetPlugin("forcefields", "MMFF94");
#else
    OBForceField *ff = OBForceField::FindType("mmff94");
#endif
    if(ff) {
        ff->SetLogLevel(OBFF_LOGLVL_LOW);
        bool success = ff->Setup(mol);
        if(!success) {
            cout << "   can't set the forcefield for molecule " << mol.GetTitle() << endl;
            exit(EXIT_FAILURE);
        }
        ff->SteepestDescent(50);
        ff->GetCoordinates(mol);
    } else {
        cout << "   Can't load force field: MMFF94!!!" << endl;
        exit(EXIT_FAILURE);
    }
}


vector<string> runMopac(string &infile, string &mopac7exe, string &keywords)
{
    cout << ">>> to process file " << infile;
    vector<string> fail_list;
    // 1. get format and base name
    string::size_type i = infile.rfind(".");
    if(i == string::npos) {
        cout << endl << "Error: can't recognize the format of " << infile << endl;
        fail_list.push_back(infile);
        return fail_list;
    }
    string format = infile.substr(i+1);
    string basic_name = infile.substr(0,i);
    // 2. read and process molecles
    OBMol mol;
    OBConversion conv;
    conv.SetInAndOutFormats(format.c_str(), "mop");
    bool next = conv.ReadFile(&mol, infile.c_str());
    int count = 0;
    while(next) {
        ++count;
        cout << endl << "    " << count << "th molecule... ";
        ostringstream os;
        os << basic_name << "_" << count << ".mopout";
        string mopout_name = os.str();
        // 2.1 add hydrogens or generate 3D coordinates
        if(!mol.HasHydrogensAdded())
            mol.AddHydrogens();
        if(!mol.Has3D())
            make3D(mol);
        // 2.2 generate FOR005
        char keys[256];
        sprintf(keys,
                "ESP MMOK VECTORS BONDS PI PRECISE ENPART AM1 GNORM=0.01 EF XYZ MULLIK +\nCharge=%d%s",
                mol.GetTotalCharge(), keywords.c_str());
        conv.AddOption("k",OBConversion::OUTOPTIONS, keys);
        ofstream outf("FOR005");
        outf << conv.WriteString(&mol);
        outf.close();
        // 2.3 run Mopac
        string cmd;
        cmd = mopac7exe + "> " + mopout_name;
        int status = system(cmd.c_str());
        if(status) {
            fail_list.push_back(mopout_name);
            cout << "failure";
        }
        else {
            ifstream inf(mopout_name.c_str());
            string line, lastline;
            while(getline(inf,line)) lastline=line;
            if(lastline.find("DONE") == string::npos) {
                fail_list.push_back(mopout_name);
                cout << "failure";
            } else
                cout << "success";
            inf.close();
        }
        system("del FOR005 FOR009 FOR010 FOR011 FOR012 SHUTDOWN");
        // next molecule
        next = conv.Read(&mol);
    }
    if(count == 0) {
        cout << "NO molecule being found!";
        fail_list.push_back(infile);
    }
    cout << endl;

    return fail_list;
}

