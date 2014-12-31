/*=============================================================================
#     FileName: tools.cpp
#         Desc: 
#       Author: jlpeng
#        Email: jlpeng1201@gmail.com
#     HomePage: 
#      Version: 0.0.1
#      Created: 2012-12-06 21:08:11
#   LastChange: 2014-03-13 20:55:24
#      History:
=============================================================================*/
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <numeric>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <openbabel/atom.h>
#include <openbabel/mol.h>
#include "tools.h"
using std::cout;
using std::cerr;
using std::endl;
using std::ifstream;
using std::vector;
using std::string;
using std::accumulate;
using namespace OpenBabel;

// used to split a given string by space
vector<string> split(const string &s)
{
    vector<string> ss;
    string::size_type i=0;
    while(i < s.size()) {
        while(i<s.size() && s[i]==' ') ++i;
        if(i == s.size()) break;
        string::size_type j = s.find(" ",i);
        if(j != string::npos) {
            ss.push_back(s.substr(i,j-i));
            i = j+1;
        } else {
            ss.push_back(s.substr(i));
            i = s.size();
        }
    }
    return ss;
}

vector<string> split(const string &s, string sep)
{
    vector<string> ss;
    string::size_type i = 0;
    while(i < s.size()) {
        string::size_type j = s.find(sep,i);
        if(j != string::npos) {
            if(j != i)
                ss.push_back(s.substr(i,j-i));
            i = j + sep.size();
        } else {
            ss.push_back(s.substr(i));
            i = s.size();
        }
    }
    return ss;
}

// remove spaces both in the begining and end of the string
string strip(const string &s)
{
    string::size_type begin_i, end_i;
    string::size_type i = 0;
    while(s[i]==' ' || s[i]=='\t' || s[i]=='\n')
        ++i;
    begin_i = i;
    i = s.size()-1;
    while(s[i]==' ' || s[i]=='\t' || s[i]=='\n')
        --i;
    end_i = i;
    return s.substr(begin_i,end_i-begin_i+1);
}

// ref: http://en.wikipedia.org/wiki/Electronegativity
double paulingEN(OBAtom &atom)
{
    double en[55] = {0.,
        2.20, 0.,
        0.98, 1.57, 2.04, 2.55, 3.04, 3.44, 3.98, 0.,
        0.93, 1.31, 1.61, 1.90, 2.19, 2.58, 3.16, 0.,
        0.82, 1.00
    };
    en[26] = 1.83;  // Fe
    en[30] = 1.65;  // Zn
    en[35] = 2.96;  // Br
    en[53] = 2.66;  // I
    int atomic_number = atom.GetAtomicNum();
    if(atomic_number>55 || en[atomic_number]==0.) {
        cout << "<tools.cpp::paulingEN>Error: unsupport atom type " << atom.GetType() << endl;
        return 0.;
    } else
        return en[atomic_number];
}

// ref: R. T. Sanderson. J. Am. Chem. Soc. 1983, 105, 2259-2261
double sandersonEN(OBAtom &atom)
{
    double en[55] = {0.,
        2.592, 0.,
        0.670, 1.810, 2.275, 2.746, 3.194, 3.654, 4.000, 0.,
        0.560, 1.318, 1.714, 2.138, 2.515, 2.957, 3.475, 0.,
        0.445, 0.946
    };
    en[30] = 2.223;  // Zn
    en[35] = 3.219;  // Br
    en[53] = 2.778;  // I
    int atomic_number = atom.GetAtomicNum();
    if(atomic_number>55 || en[atomic_number]==0.) {
        cout << "<tools.cpp::sandersonEN>Error: unsupport atom type " << atom.GetType() << endl;
        return 0.;
    } else
        return en[atomic_number];
}
// alpha value used for calculating Kier shape index (Kappa shape index)
// source: http://www.edusoft-lc.com/molconn/manuals/400/appII.html
// NOT USED ANY MORE
double kierAlpha(OBAtom &atom)
{
    if(atom.IsHydrogen())
        return -0.52;
    if(atom.MatchesSMARTS("C"))
        return 0.00;
    if(atom.MatchesSMARTS("[C^2,c]"))
        return -0.13;
    if(atom.MatchesSMARTS("[$(C#*),$(C(=*)=*)]"))
        return -0.22;
    if(atom.MatchesSMARTS("N"))
        return -0.04;
    if(atom.MatchesSMARTS("[n,$(N=*),$(N[#6,#7,#8]=,:,#*)]"))
        return -0.20;
    if(atom.MatchesSMARTS("[$(N#*),$([ND2](=*)=*)]"))
        return -0.29;
    if(atom.MatchesSMARTS("O"))
        return -0.04;
    if(atom.MatchesSMARTS("[o,$(O=*),$(O[#6,#7,#8]=,:*)]"))
        return -0.20;
    if(atom.GetAtomicNum() == 9) // F
        return -0.07;
    if(atom.MatchesSMARTS("P"))  // sp3 P
        return 0.43;
    if(atom.MatchesSMARTS("[#15;$([PD1]=*)]")) // sp2 P
        return 0.30;
    // sp P
    if(atom.MatchesSMARTS("S")) // sp3 S
        return 0.35;
    if(atom.MatchesSMARTS("[#16;s,$([SD1]=*)]")) //sp2 S
        return 0.22;
    // sp S
    if(atom.GetAtomicNum() == 17) // Cl
        return 0.29;
    if(atom.GetAtomicNum() == 35) // Br
        return 0.48;
    if(atom.GetAtomicNum() == 53) // I
        return 0.73;

    cerr << "tools.cpp::KierAlpha: unsupport atom type: " << atom.GetType() << endl;
    return 0;
    
}

double getVWRadius(OBAtom &atom)
{
    // reference:
    // - descriptors of molecular shape applied in studies of
    //   structure/activity and structure/property relationships, 1986
    int atomic_number = atom.GetAtomicNum();
    
    if(atomic_number == 1)
        return 1.1;
    if(atomic_number == 6) {
        if(atom.MatchesSMARTS("C"))    // SP3 C
            return 1.70;
        if(atom.MatchesSMARTS("[C^2]"))  // SP2 C
            return 1.70;
        if(atom.MatchesSMARTS("[$(C#*)]")) // SP C
            return 1.77;
        if(atom.MatchesSMARTS("[$(C(=*)=*)]")) // allenyl
            return 1.70;
        if(atom.MatchesSMARTS("c"))       // aromatic
            return 1.77;
    }
    if(atomic_number == 7) {
        if(atom.MatchesSMARTS("N"))       // SP3 nitrogen
            return 1.55;
        if(atom.MatchesSMARTS("[$(N=*),$(N[#6,#7,#8]=,:,#*)]"))  // SP2 nitrogen
            return 1.55;
        if(atom.MatchesSMARTS("[$(N#*),$([ND2](=*)=*)]"))   // SP nitrogen
            return 1.60;
        if(atom.MatchesSMARTS("[n]"))       // aromatic nitrogen
            return 1.60;
    }
    if(atomic_number == 8) {
        if(atom.MatchesSMARTS("[$([#8]*)]")) //singly bonded oxygen
            return 1.52;
        if(atom.MatchesSMARTS("[$(O=*)]"))   // doubly donded oxygen
            return 1.50;
    }
    
    if(atomic_number == 16)      //sulfur
        return 1.80;
    if(atomic_number == 9)       // F
        return 1.50;
    if(atomic_number == 17)      // Cl
        return 1.75;
    if(atomic_number == 35)      // Br
        return 1.85;
    if(atomic_number == 53)      // I
        return 1.97;

    // - if the atom doesn't match any of above SMARTS strings, then the following values will be used.
    // - the following values(angstrom, A) are from OBDATA_DIR/element.txt, which is from
    //   "Consistent van der Waals Radii for the Whole Main Group, 2009"
    double value[55] = {0.,
        1.1, 1.4,
        1.81, 1.53, 1.92, 1.70, 1.55, 1.52, 1.47, 1.54,
        2.27, 1.73, 1.84, 2.10, 1.80, 1.80, 1.75, 1.88,
        2.75, 2.31
    };
    value[26] = 2.05;
    value[30] = 2.10;
    value[35] = 1.83;
    value[53] = 1.98;
    if(atomic_number>55 || value[atomic_number]==0.) {
        cout << "<tools.cpp::getVWRadius>Error: unsupport atom type " << atom.GetType() << endl;
        return 2.0;
    }
    else
        return value[atomic_number];
}

/*
// ref: Cordero B. et al. Dalton. Trans., 2008: 2832-2838.
double getCovalentRadius(OBAtom &atom)
{
    int hyb = static_cast<int>(atom.GetHyb());
    int atomic_number = static_cast<int>(atom.GetAtomicNum());
    if(atomic_number == 1)
        return 0.31;
    if(atomic_number == 6) {
        if(hyb == 3)
            return 0.76;
        if(hyb == 2)
            return 0.73;
        if(hyb == 1)
            return 0.69;
    }
    if(atomic_number == 7)
        return 0.71;
    if(atomic_number == 8)
        return 0.66;
    if(atomic_number == 9)
        return 0.57;
    if(atomic_number == 14)
        return 1.11;
    if(atomic_number == 15)
        return 1.07;
    if(atomic_number == 16)
        return 1.05;
    if(atomic_number == 17)
        return 1.02;
    if(atomic_number == 35)
        return 1.20;
    if(atomic_number == 53)
        return 1.39;
    cerr << "<tools.cpp::getCovalentRadius>Error: unsupport atomic type: " << atom.GetType() << endl;
    return 0.;
}
*/
// parameters are from page 430 of "molecular descriptors for chemoinformatics. 2 edn."
double getCovalentRadius(OBAtom &atom)
{
    int hyb = static_cast<int>(atom.GetHyb());
    int atomic_number = static_cast<int>(atom.GetAtomicNum());
    if(atomic_number == 6) {
        if(hyb == 3)
            return 0.77;
        if(hyb == 2)
            return 0.67;
        if(hyb == 1)
            return 0.60;
    }
    if(atomic_number == 7) {
        if(hyb == 3)
            return 0.74; // 0.725
        if(hyb == 2)
            return 0.62; // 0.625
        if(hyb == 1)
            return 0.55;
    }
    if(atomic_number == 8) {
        if(hyb == 3)
            return 0.74;
        if(hyb == 2)
            return 0.62; // 0.61
    }
    if(atomic_number == 9)
        return 0.72; // 0.71
    if(atomic_number == 15) {
        // return 1.06;
        if(hyb == 3)
            return 1.10;
        else if(hyb == 2)
            return 1.00;
        else
            return 1.06;
    }
    if(atomic_number == 16) {
        if(hyb == 3)
            return 1.04; // 1.035
        if(hyb == 2)
            return 0.94; // 0.945
    }
    if(atomic_number == 17)
        return 0.99;
    if(atomic_number == 35)
        return 1.14;
    if(atomic_number == 53)
        return 1.33;
    // the following are taken from "Cordero B. et al. Dalton. Trans., 2008: 2832-2838"
    double value[55] = {0.,
        0.31, 0.28,
        1.28, 0.96, 0.84, 0.00, 0.71, 0.66, 0.57, 0.58,
        1.66, 1.41, 1.21, 1.11, 1.07, 1.05, 1.02, 1.06,
        2.03, 1.76
    };
    value[30] = 1.22; // Zn
    value[35] = 1.20; // Br
    value[53] = 1.39; // I
    if(atomic_number>55 || value[atomic_number]==0.) {
        cout << "<tools.cpp::getCovalentRadius>Error: unsupport atom type: " << atom.GetType() << endl;
        return 0.;
    } else
        return value[atomic_number];
}

int getValElectrons(OBAtom &atom)
{
    int value[55] = {0,
        1, 0,
        1, 2, 3, 4, 5, 6, 7, 0,
        1, 2, 3, 4, 5, 6, 7, 0,
        1, 2
    };
    value[35] = 7;
    value[53] = 7;
    int atomic_number = atom.GetAtomicNum();
    if(atomic_number>55 || value[atomic_number]==0) {
        cout << "<tools.cpp::getValElectrons>Error: unsupport atom type " << atom.GetType() << endl;
        return -1;
    } else
        return value[atomic_number];
}

vector<vector<double> > generate_sphere_points(int n)
{
    vector<vector<double> > points;
    double inc = PI * (3 - sqrt(5.));
    double offset = 2.0 / n;
    for(int k=0; k<n; ++k) {
        double y = k * offset - 1 + (offset / 2);
        double r = sqrt(1 - y*y);
        double phi = k * inc;
        vector<double> tmp;
        tmp.push_back(cos(phi)*r);
        tmp.push_back(y);
        tmp.push_back(sin(phi)*r);
        points.push_back(tmp);
    }
    return points;
}

// static method called by 'overlap'
static bool intersect(double point1[2], double point2[2], double circle[3])
{
    //1. if p is the endpoint of the line
    if((fabs(circle[0]-point1[0])<=1e-4 && fabs(circle[1]-point1[1])<=1e-4) ||
            (fabs(circle[0]-point2[0])<=1e-4 && fabs(circle[1]-point2[1])<=1e-4))
        return true;
    //calculate the distance from the center of circle to the line segment
    double pa = sqrt(pow(circle[0]-point1[0],2)+pow(circle[1]-point1[1],2));
    double pb = sqrt(pow(circle[0]-point2[0],2)+pow(circle[1]-point2[1],2));
    double ab = sqrt(pow(point1[0]-point2[0],2)+pow(point1[1]-point2[1],2));
    double cos_pab = (pa*pa+ab*ab-pb*pb)/(2*pa*ab);
    double cos_pba = (ab*ab+pb*pb-pa*pa)/(2*ab*pb);
    double d;
    if(cos_pab <= 0)
        d = pa;
    else if(cos_pba <= 0)
        d = pb;
    else {
        double a = point2[1] - point1[1];
        double b = point1[0] - point2[0];
        double c = point2[0]*point1[1] - point1[0]*point2[1];
        d = fabs(a*circle[0]+b*circle[1]+c)/sqrt(a*a+b*b);
    }

    if(circle[2] - d > 1e-4)
        return true;
    else
        return false;
}
// the four points of rectangle should be sequential!!
// a--b
// |  |
// d--c
bool overlap(double rectangle[][2], double circle[3])
{
    // 1. get the center of the rectangle
    double center_x = (rectangle[0][0]+rectangle[2][0])/2.0;
    double center_y = (rectangle[0][1]+rectangle[2][1])/2.0;
    double edge01 = sqrt(pow(rectangle[0][0]-rectangle[1][0],2)+
            pow(rectangle[0][1]-rectangle[1][1],2));
    double edge03 = sqrt(pow(rectangle[0][0]-rectangle[3][0],2)+
            pow(rectangle[0][1]-rectangle[3][1],2));
    double edge_short = edge01<edge03 ? edge01 : edge03;
    double diagonal = sqrt(pow(rectangle[0][0]-rectangle[2][0],2)+
            pow(rectangle[0][1]-rectangle[2][1],2));
    // 2. calculate the distance between the center of the circle and the center of the rectangle
    double d = sqrt(pow(center_x-circle[0],2)+pow(center_y-circle[1],2));
    // if center of the rectangle is within the range of circle (x0,y0,r+edge_short/2.0) ==> True
    if((circle[2]+edge_short/2.0)-d > 1e-4)
        return true;
    // if center of the rectangle is out of the range of circle (x0,y0,r+diagonal/2.0)   ==> False
    if((circle[2]+diagonal/2.0)-d < 1e-4)
        return false;
    // otherwise, for each eadge, to see if it intersect the circle or in the circle
    // 0---1, 1---2, 2---3
    for(int i=0; i<3; i++) {
        if(intersect(rectangle[i],rectangle[i+1],circle))
            return true;
    }
    // 0---3
    if(intersect(rectangle[0],rectangle[3],circle))
        return true;

    return false;
}

bool overlap3D(double cube[4], double sphere[4])
{
    /*
     * cube: (x,y,z,d)
     *     cube whose surfaces are parallel to axes
     *     a_____b
     *   d/____c/|  a(x,y,z), b(x+d,y,z), c(x+d,y+d,z), d(x,y+d,z)
     *    |____|/
     * sphere: (x,y,z,r)
     */
    double points[4][2];
    double circle[3];
    double radius[2];
    // center of the cube
    double x0 = (2*cube[0]+cube[3])/2.0;
    double y0 = (2*cube[1]+cube[3])/2.0;
    double z0 = cube[2] - cube[3]/2.0;
    // distance from center of the cube to that of sphere
    double d = sqrt(pow(x0-sphere[0],2)+pow(y0-sphere[1],2)+pow(z0-sphere[2],2));
    // if the center of the cube is within the sphere ==> true
    if ((cube[3]/2.0+sphere[3]) - d > 1e-4)
        return true;
    // half of the diagonal length
    double diagonal_2 = sqrt(pow(x0-cube[0],2)+pow(y0-cube[1],2)+pow(z0-cube[2],2));
    if((diagonal_2+sphere[3]) - d < 1e-4)
        return false;

    // projected face <<if overlap??>> projected cross section of sphere by the face
    /*
     * Type 1 -- upper/bottom face
     * y  d----c  known a(x,y)
     * ^  |    |  ==> b(x+d,y), c(x+d,y+d), z(x,d+y)
     * |  a----b
     * +---->x
     */
    radius[0]=0.; radius[1]=0.;
    // radius of cross section of sphere by upper face  (z_sphere < z < z_sphere+r)
    if (((sphere[2]+sphere[3]) - cube[2] > 1e-4) && (cube[2] - (sphere[2]-sphere[3]) > 1e-4))
        radius[0] = sqrt(pow(sphere[3],2)-pow(sphere[2]-cube[2],2));
    // radius of cross section of sphere by bottom face
    if (((sphere[2]+sphere[3])-(cube[2]-cube[3]) > 1e-4) &&
           ((cube[2]-cube[3])-(sphere[2]-sphere[3]) > 1e-4))
        radius[1] = sqrt(pow(sphere[3],2)-pow(sphere[2]-(cube[2]-cube[3]),2));
    for(int i=0; i<2; ++i) {
        if(radius[i] == 0.)
            continue;
        points[0][0] = cube[0]; points[0][1] = cube[1];
        points[1][0] = cube[0] + cube[3]; points[1][1] = cube[1];
        points[2][0] = cube[0] + cube[3]; points[2][1] = cube[1] + cube[3];
        points[3][0] = cube[0]; points[3][1] = cube[1] + cube[3];
        circle[0] = sphere[0]; circle[1] = sphere[1]; circle[2] = radius[i];
        if(overlap(points,circle))
            return true;
    }

    /*
     * Type 2 -- left/right face
     * z a----d  known a(y,z)
     * ^ |    |  ==> b(y,z-d), c(y+d,z-d), d(y+d,z)
     * | b----c
     * +--->y
     */
    radius[0]=0.; radius[1]=0.;
    // radius of cross section of sphere by left face
    if (((sphere[0]+sphere[3]) - cube[0] > 1e-4) && (cube[0] - (sphere[0]-sphere[3]) > 1e-4))
        radius[0] = sqrt(pow(sphere[3],2)-pow(sphere[0]-cube[0],2));
    // radius of cross section of sphere by right face
    if (((sphere[0]+sphere[3]) - (cube[0]+cube[3]) > 1e-4) &&
           ((cube[0]+cube[3]) - (sphere[0]-sphere[3]) > 1e-4))
        radius[1] = sqrt(pow(sphere[3],2)-pow(sphere[0]-(cube[0]+cube[3]),2));
    for(int i=0; i<2; ++i) {
        if(radius[i] == 0.)
            continue;
        points[0][0] = cube[1]; points[0][1] = cube[2];
        points[1][0] = cube[1]; points[1][1] = cube[2] - cube[3];
        points[2][0] = cube[1] + cube[3]; points[2][1] = cube[2] - cube[3];
        points[3][0] = cube[1] + cube[3]; points[3][1] = cube[2];
        circle[0] = sphere[1]; circle[1] = sphere[2]; circle[2] = radius[i];
        if(overlap(points,circle))
            return true;
    }

    /*
     * Type 3 -- back/front face
     * z  a----b   known that a (x,z)
     * ^  |    |   => b(x+d,z), c(x+d,z-d), d(x,z-d)
     * |  d----c
     * +------>x
     */
    radius[0]=0.; radius[1]=0.;
    // radius of cross section of sphere by back face  (y_sphere-r < y < y_sphere+r)
    if (((sphere[1]+sphere[3]) - cube[1] > 1e-4) && (cube[1] - (sphere[1]-sphere[3]) > 1e-4))
        radius[0] = sqrt(pow(sphere[3],2)-pow(sphere[1]-cube[1],2));
    // radius of cross section of sphere by front face
    if (((sphere[1]+sphere[3]) - (cube[1]+cube[3]) > 1e-4) &&
            ((cube[1]+cube[3]) - (sphere[1]-sphere[3]) > 1e-4))
        radius[1] = sqrt(pow(sphere[3],2)-pow(sphere[1]-(cube[1]+cube[3]),2));
    for(int i=0; i<2; ++i) {
        if(radius[i] == 0.)
            continue;
        points[0][0] = cube[0]; points[0][1] = cube[2];
        points[1][0] = cube[0] + cube[3]; points[1][1] = cube[2];
        points[2][0] = cube[0] + cube[3]; points[2][1] = cube[2] - cube[3];
        points[3][0] = cube[0]; points[3][1] = cube[2] - cube[3];
        circle[0] = sphere[0]; circle[1] = sphere[2]; circle[2] = radius[i];
        if(overlap(points,circle))
            return true;
    }

    // all test failed
    return false;
}


int getPrincipalQuantumNumber(OBAtom &atom)
{
    int value[55] = {0,
        1, 1,
        2, 2, 2, 2, 2, 2, 2, 2,
        3, 3, 3, 3, 3, 3, 3, 3,
        4, 4
    };
    value[35] = 4;
    value[53] = 5;
    int atomic_number = atom.GetAtomicNum();
    if(atomic_number>55 || value[atomic_number]==0) {
        cout << "<tools.cpp::getPrincipalQuantumNumber>Error: unsupport atom type " << atom.GetType() << endl;
        return -1;
    } else
        return value[atomic_number];
}


