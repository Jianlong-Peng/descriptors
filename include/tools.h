/*=============================================================================
#     FileName: tools.h
#         Desc: some functions used by shadow and volume
#       Author: jlpeng
#        Email: jlpeng1201@gmail.com
#     HomePage: 
#      Version: 0.0.1
#      Created: 2012-12-06 21:06:49
#   LastChange: 2014-01-08 12:33:02
#      History:
=============================================================================*/


#ifndef  TOOLS_H
#define  TOOLS_H

#include <openbabel/babelconfig.h>
#include <openbabel/atom.h>
#include <openbabel/mol.h>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include "config.h"

// used to join words in [first,last) by 'sep'
template <typename InputIterator>
std::string joinBy(InputIterator first, InputIterator last, std::string sep)
{
    if(first == last) return std::string();
    std::ostringstream os;
    os << *first;
    while(++first != last)
        os << sep << *first;
    return os.str();
}
// used to split a given string by space
std::vector<std::string> DESCRIPTOR_API split(const std::string &s);
// overloaded method to split a given string by a specified separator
// continous separators will be removed!
std::vector<std::string> DESCRIPTOR_API split(const std::string &s, std::string sep);
std::string DESCRIPTOR_API strip(const std::string &s);
// get Pauling electronegativity
// if invalid atom type, return 0.
double DESCRIPTOR_API paulingEN(OpenBabel::OBAtom &atom);
// get Sanderson electronegativity
// if invalid atom type, return 0.
double DESCRIPTOR_API sandersonEN(OpenBabel::OBAtom &atom);
// to get Van der Waals radius
double DESCRIPTOR_API getVWRadius(OpenBabel::OBAtom &atom);
// to get atomic radius
// if invalid atom type, return 0.
double DESCRIPTOR_API getCovalentRadius(OpenBabel::OBAtom &atom);
// to get valence electrons of a given atom
// if invalid atom type, return -1
int DESCRIPTOR_API getValElectrons(OpenBabel::OBAtom &atom);
// to get the principal quantum number of a given atom;
// if invalid atom type, return -1
int DESCRIPTOR_API getPrincipalQuantumNumber(OpenBabel::OBAtom &atom);
// to get number of hydrogen neighbors of a given atom;
inline int DESCRIPTOR_API getNumOfHNbors(OpenBabel::OBAtom &atom) {
    return static_cast<int>(atom.GetImplicitValence() - atom.GetHvyValence());
}
// calculate Euclidean distance
inline double DESCRIPTOR_API pos_distance_sq(const double *p1, const double *p2) {
    return (pow(p1[0]-p2[0],2)+pow(p1[1]-p2[1],2)+pow(p1[2]-p2[2],2));
}
inline double DESCRIPTOR_API pos_distance(const double *p1, const double *p2) {
    return sqrt(pos_distance_sq(p1,p2));
}
// Returns list of 3d coordinates of points on a sphere using the
// Golden Section Spiral algorithm.
std::vector<std::vector<double> > DESCRIPTOR_API generate_sphere_points(int n);
// OBJ: to see if part of or all the rectangle is within the circle
//      used by 'shadow'
// Attention: the size of rectangle is smaller than that of circle!!!
bool DESCRIPTOR_API overlap(double rectangle[][2], double circle[3]);
// OBJ: to see if part of or all the cube is within the sphere
//      used by 'volume'
// Attention: the size of cube is smaller than that of sphere!!
//            cube's surfaces are parallel to axes!!
bool DESCRIPTOR_API overlap3D(double cube[4], double sphere[4]);

#endif   /* ----- #ifndef TOOLS_H  ----- */

