/*=============================================================================
#     FileName: physicochemical.cpp
#         Desc: 
#       Author: jlpeng
#        Email: jlpeng1201@gmail.com
#     HomePage: 
#      Created: 2013-11-30 13:59:45
#   LastChange: 2013-12-27 16:52:43
#      History:
=============================================================================*/
#include <iostream>
#include <vector>
#include <string>
#include <openbabel/mol.h>
#include <openbabel/groupcontrib.h>
#include "config.h"
#include "physicochemical.h"

using std::vector;
using std::string;
using namespace OpenBabel;


double calcLogP(OBMol &mol)
{
    double logp;
#if defined(_MSC_VER)
    OBDescriptor* obLogP = (OBDescriptor *)OBPlugin::GetPlugin("descriptors", "logP");
#else
    OBDescriptor *obLogP = OBDescriptor::FindType("logP");
#endif
    if(obLogP)
        logp = obLogP->Predict(&mol);
    else {
        std::cerr << "Error: can't calculate logP for molecule " << mol.GetTitle() << std::endl;
        logp = 1E+38;
    }
    return logp;
}

