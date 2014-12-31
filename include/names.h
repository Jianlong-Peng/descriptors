/*=============================================================================
#     FileName: names.h
#         Desc: 
#       Author: jlpeng
#        Email: jlpeng1201@gmail.com
#     HomePage: 
#      Created: 2013-12-27 12:34:05
#   LastChange: 2013-12-27 16:41:18
#      History:
=============================================================================*/

#ifndef  DES_NAMES_H
#define  DES_NAMES_H
#include <vector>
#include <string>
#include "config.h"

extern std::vector<std::string> DESCRIPTOR_API constitute_names;
extern std::vector<std::string> DESCRIPTOR_API electrostatic_names;
extern std::vector<std::string> DESCRIPTOR_API geometrical_names;
extern std::vector<std::string> DESCRIPTOR_API topological_names;
extern std::vector<std::string> DESCRIPTOR_API qc_names;
extern std::vector<std::string> DESCRIPTOR_API physicochemical_names;
// names for charged partial surface area (CPSA)
// different kind of atomic charge could yield various such descriptors
extern std::vector<std::string> DESCRIPTOR_API cpsa_names;

extern std::string DESCRIPTOR_API reference;

#endif   /* ----- #ifndef DES_NAMES_H  ----- */

