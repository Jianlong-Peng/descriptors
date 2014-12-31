/*=============================================================================
#     FileName: logP.h
#         Desc: to predict logP (Octanol-water partition coefficient)
#       Author: jlpeng
#        Email: jlpeng1201@gmail.com
#     HomePage: 
#      Version: 0.0.1
#      Created: 2013-01-24 21:20:15
#   LastChange: 2013-02-20 19:29:11
#      History:
=============================================================================*/


#ifndef  LOGP_H
#define  LOGP_H
//#include <vector>
//#include <utility>
#include <string>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include "config.h"

struct DESCRIPTOR_API AlogPRecord {
    std::string pattern;
    double value;
    struct AlogPRecord *next;
};
struct DESCRIPTOR_API AlogPHead {
    struct AlogPRecord *last; // point to the currently last record
    struct AlogPRecord *record;
};

class DESCRIPTOR_API AlogP
{
public:
    AlogP(): frag_logP(NULL),num(0) {}
    ~AlogP();
    explicit AlogP(std::string para_file) {read_parameter(para_file);}
    double predict(OpenBabel::OBMol &mol);
    double predict_sym(OpenBabel::OBMol &mol);
    void read_parameter(std::string para_file);
private:
    // atomic contribution for a given type.
    double sub_predict(OpenBabel::OBAtom &atom, int i);
private:
    // Default: H, SP3 C, SP2 C, SP C, ar C, nitrogen, oxygen, F, Cl, Br, I, others
    AlogPHead *frag_logP;
    int num;
    // hydrogen,carbon,nitrogen,oxygen,halogen,others
    //std::vector<std::vector<std::pair<std::string,double> > > frag_logP;
};

#endif   /* ----- #ifndef LOGP_H  ----- */
