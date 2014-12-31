/*=============================================================================
#     FileName: geometrical.cpp
#         Desc: 
#       Author: jlpeng
#        Email: jlpeng1201@gmail.com
#     HomePage: 
#      Version: 0.0.1
#      Created: 2012-12-28 20:38:17
#   LastChange: 2013-12-27 14:29:30
#      History:
=============================================================================*/
#include <iostream>
#include <vector>
#include <utility>
#include <algorithm>
#include <numeric>
#include <string>
#include <cassert>
#include <openbabel/mol.h>
#include <openbabel/obiter.h>
#include "geometrical.h"
#include "tools.h"
using std::cerr;
using std::endl;
using std::vector;
using std::pair;
using std::make_pair;
using std::sort;
using std::string;
using std::accumulate;
using namespace OpenBabel;

/*
static string str_names[] = {"MOIA","MOIB","MOIC",
    "XYShadow","ZXShadow","YZShadow","XYShadow/XYR","ZXShadow/ZXR","YZShadow/YZR",
    "MV","MV/XYZBox",
    "grav-pairs","grav-bonds",
    "SASA_sum", "TPSA"};
vector<string> DESCRIPTOR_API geometry_names(str_names,str_names+15);
*/

// used by 'calcMOI'
bool comp(const pair<int,double> &v1, const pair<int, double> &v2)
{
    return v1.second < v2.second;
}
vector<pair<int,double> > calcMOI(OBMol &mol)
{
    //std::cout << "ToInertialFrame()...";
    mol.ToInertialFrame();
    double MOIA=0., MOIB=0., MOIC=0.;
    double mi,x,y,z;
    //std::cout << std::endl << "to calculate MOI..." << std::endl;
    FOR_ATOMS_OF_MOL(atom, mol) {
        mi = atom->GetAtomicMass();
        x = atom->GetX();
        y = atom->GetY();
        z = atom->GetZ();
        //std::cout << x << " " << y << " " << z << std::endl;
        MOIA += mi*(y*y+z*z);
        MOIB += mi*(x*x+z*z);
        MOIC += mi*(x*x+y*y);
    }
    vector<pair<int, double> > result;
    result.push_back(make_pair(0,MOIA));
    result.push_back(make_pair(1,MOIB));
    result.push_back(make_pair(2,MOIC));
    sort(result.begin(),result.end(),comp);
    return result;
}

// static methods called by 'shadow' and volume'
static double getCoord(OBAtom &atom, int index)
{
    if(index==0)
        return atom.GetX();
    else if(index==1)
        return atom.GetY();
    else if(index==2)
        return atom.GetZ();
    else {
	cerr << "shadow.cpp::get_coord: index should be one of (0, 1, 2), but " 
            << index << " is passed!" << endl;
	exit(EXIT_FAILURE);
    }
}
// static methods called by 'shadow' and volume'
static int toInt(double value)
{
    if(value - int(value) > 0.5)
        return int(value)+1;
    else
        return int(value);
}
// static methods called by 'shadow' and volume'
static void refineBorder(double border[], int n, int density)
{
    for(int m=0; m<n; ++m) {
        if(border[m] < 0)
            border[m] = 1.0*floor(border[m]*density)/density;
        else
            border[m] = 1.0*ceil(border[m]*density)/density;
    }
}
// static mehtod called by 'shadow'
static void getBorderShadow(OBMol &mol, int i, int j, int density, double border[4])
{
    border[0] = 1e8; border[1] = -1e-8;
    border[2] = 1e8; border[3] = -1e-8;
    FOR_ATOMS_OF_MOL(atom,mol) {
        double x = getCoord(*atom,i);
        double y = getCoord(*atom,j);
        double r = getVWRadius(*atom);
        double Xmin = x-r;
        double Xmax = x+r;
        double Ymin = y-r;
        double Ymax = y+r;
        if(border[0] > Xmin)
            border[0] = Xmin;
        if(border[1] < Xmax)
            border[1] = Xmax;
        if(border[2] > Ymin)
            border[2] = Ymin;
        if(border[3] < Ymax)
            border[3] = Ymax;
    }
    
    refineBorder(border,4,density);
}
// int axes[3]: principle inertia axes (MOI: axes[0] < axes[1] < axes[2])
//   0: X-axis; 1: Y-axis; 2: Z-axis
// if axes[0]==-1, then mol will be aligned to principle inertia
// otherwise, using XYZ-axes indicated by axes
// return: [XYshadow,ZXShadow,YZShadow,XYShadow/XYR,ZXShadow/ZXR,YZShadow/YZR]
void shadow(OBMol &mol, int axes[3], int density, double result[6])
{
    //initialize_radius();
    int index[3];  // principle axes
    //align the molecule
    if(axes[0] == -1) {
        vector<pair<int,double> > moi = calcMOI(mol);
        index[0] = moi[0].first;
        index[1] = moi[1].first;
        index[2] = moi[2].first;
    } else {
        index[0] = axes[0];
        index[1] = axes[1];
        index[2] = axes[2];
    }
    //XYShadow, ZXShadow, YZShadow
    bool *cover = new bool[200*200];
    int total = 200*200;
    double area_grid = pow(1.0/density,2);
    double border[4], border_atom[4], circle[3], coords[4][2];
    for(int i=0; i<2; i++) {
        for(int j=i+1; j<3; j++) {
            //1st: get the border (xmin,xmax,ymin,ymax) of the molecule's shadow
            getBorderShadow(mol,index[i],index[j],density,border);
            double rectangle_area = (border[1]-border[0])*(border[3]-border[2]);
            //2nd: for each atom, to find possible grids within the atom's external rectangle
            int num_grids_x = toInt((border[1]-border[0])*density);
            int num_grids_y = toInt((border[3]-border[2])*density);
            if(num_grids_x*num_grids_y > total) {
                total = num_grids_x*num_grids_y;
                delete[] cover;
                cover = new bool[num_grids_x*num_grids_y];
                if(cover == NULL) {
                    cerr << "shadow.cpp::calcShadow: out of memory!" << endl;
                    exit(1);
                }
            }
            memset(cover,0,sizeof(bool)*total);
            FOR_ATOMS_OF_MOL(atom,mol) {
                //2.1. to get the external rectangle of the atom
                circle[0] = getCoord(*atom, index[i]);
                circle[1] = getCoord(*atom, index[j]);
                circle[2] = getVWRadius(*atom);
                border_atom[0] = circle[0] - circle[2]; border_atom[1] = circle[0] + circle[2];
                border_atom[2] = circle[1] - circle[2]; border_atom[3] = circle[1] + circle[2];
                refineBorder(border_atom,4,density);
                //2.2 to check all grids within the external rectangle
                int num_grids_atom_x = toInt((border_atom[1]-border_atom[0])*density);
                int num_grids_atom_y = toInt((border_atom[3]-border_atom[2])*density);
                for(int m=0; m<num_grids_atom_y; ++m) {
                    double y = border_atom[3] - 1.0*m/density;
                    int cover_y = toInt((border[3]-y)*density);
                    for(int n=0; n<num_grids_atom_x; ++n) {
                        double x = border_atom[0] + 1.0*n/density;
                        int cover_x = toInt((x-border[0])*density);
                        assert(num_grids_y*cover_x+cover_y < total);
                        if(!cover[num_grids_y*cover_x+cover_y]) {
                            coords[0][0] = x; coords[0][1] = y;
                            coords[1][0] = x+1.0/density; coords[1][1] = y;
                            coords[2][0] = x+1.0/density; coords[2][1] = y-1.0/density;
                            coords[3][0] = x; coords[3][1] = y-1.0/density;
                            if(overlap(coords,circle))
                                cover[num_grids_y*cover_x+cover_y] = true;
                        }
                    }
                } // END of checking all grids within the external rectangle of the atom
            } // checking all atoms
            int total_covered = 0;
            for(int m=0; m<total; ++m)
                if(cover[m]) ++total_covered;
            result[i+j-1] = total_covered*area_grid;
            result[i+j-1+3] = result[i+j-1]/rectangle_area;
        }
    }
    delete[] cover;
}

// static method called by 'volume'
// to get the external parallelipide box for a given molecule
static void getBorderVolume(OBMol &mol, int index[3], int density, double border[6])
{
    double x,y,z,r,Xmin,Xmax,Ymin,Ymax,Zmin,Zmax;
    border[0] = 1e8; border[1] = -1e8;
    border[2] = 1e8; border[3] = -1e8;
    border[4] = 1e8; border[5] = -1e8;
    FOR_ATOMS_OF_MOL(atom,mol) {
        x = getCoord(*atom,index[0]);
        y = getCoord(*atom,index[1]);
        z = getCoord(*atom,index[2]);
        r = getVWRadius(*atom);
        Xmin = x-r; Xmax = x+r;
        Ymin = y-r; Ymax = y+r;
        Zmin = z-r; Zmax = z+r;
        if(Xmin < border[0])
            border[0] = Xmin;
        if(Xmax > border[1])
            border[1] = Xmax;
        if(Ymin < border[2])
            border[2] = Ymin;
        if(Ymax > border[3])
            border[3] = Ymax;
        if(Zmin < border[4])
            border[4] = Zmin;
        if(Zmax > border[5])
            border[5] = Zmax;
    }
    refineBorder(border,6,density);
}
void volume(OBMol &mol, int axes[3], int density, double result[2])
{
    double cube[4];
    double sphere[4];
    double grid_volume = pow(1.0/density,3);
    int index[3];
    if(axes[0] == -1) {
        vector<pair<int,double> > moi = calcMOI(mol);
        index[0] = moi[0].first;
        index[1] = moi[1].first;
        index[2] = moi[2].first;
    } else {
        index[0] = axes[0];
        index[1] = axes[1];
        index[2] = axes[2];
    }

    // step 1. to get the external parallelipide box
    double border[6];
    getBorderVolume(mol,index,density,border);
    double box_volume = (border[1]-border[0])*(border[3]-border[2])*(border[5]-border[4]);

    // step 2. for each atom, to check cubes within the atom's external box
    int num_grids_x = toInt((border[1]-border[0])*density);
    int num_grids_y = toInt((border[3]-border[2])*density);
    int num_grids_z = toInt((border[5]-border[4])*density);
    bool *cover = new bool[num_grids_x*num_grids_y*num_grids_z];
    if(cover == NULL) {
        cerr << "volume.cpp:calcVolume: out of memory!" << endl;
        exit(EXIT_FAILURE);
    }
    memset(cover,0,sizeof(bool)*num_grids_x*num_grids_y*num_grids_z);
    FOR_ATOMS_OF_MOL(atom,mol) {
        double border_atom[6];
        sphere[0] = getCoord(*atom,index[0]);
        sphere[1] = getCoord(*atom,index[1]);
        sphere[2] = getCoord(*atom,index[2]);
        sphere[3] = getVWRadius(*atom);
        border_atom[0] = sphere[0]-sphere[3]; border_atom[1] = sphere[0]+sphere[3];
        border_atom[2] = sphere[1]-sphere[3]; border_atom[3] = sphere[1]+sphere[3];
        border_atom[4] = sphere[2]-sphere[3]; border_atom[5] = sphere[2]+sphere[3];
        refineBorder(border_atom,6,density);
        int num_grids_atom_x = toInt((border_atom[1]-border_atom[0])*density);
        int num_grids_atom_y = toInt((border_atom[3]-border_atom[2])*density);
        int num_grids_atom_z = toInt((border_atom[5]-border_atom[4])*density);
        for(int k=0; k<num_grids_atom_z; ++k) {
            double z = border_atom[5] - 1.0*k/density;
            int cover_z = toInt((border[5]-z)*density);
            for(int j=0; j<num_grids_atom_y; ++j) {
                double y = border_atom[2] + 1.0*j/density;
                int cover_y = toInt((y-border[2])*density);
                for(int i=0; i<num_grids_atom_x; ++i) {
                    double x = border_atom[0] + 1.0*i/density;
                    int cover_x = toInt((x-border[0])*density);
                    int _index = cover_z*(num_grids_x*num_grids_y)+cover_x*num_grids_y+cover_y;
                    assert(_index < num_grids_x*num_grids_y*num_grids_z);
                    if(!cover[_index]) {
                        cube[0] = x; cube[1] = y; cube[2] = z; cube[3] = 1.0/density;
                        if(overlap3D(cube,sphere))
                            cover[_index] = true;
                    }
                }
            }
        }
    } // END of visting all atoms
    int total_covered=0;
    for(int i=0; i<num_grids_x*num_grids_y*num_grids_z; ++i)
        if(cover[i]) ++total_covered;
    result[0] = total_covered * grid_volume;
    result[1] = result[0] / box_volume;
    delete[] cover;
}

// Molecular surface area
double surfaceArea(OBMol &mol)
{
    cerr << "geometrical.cpp::surfaceArea: haven't been implemented!!" << endl;
    return 0;
}


// - static method called by 'calcSASA'
// - from radii.txt of asa.zip (http://boscoh.com/protein/asapy)
static double getRadius(OBAtom *atom)
{
    int atomic_number = atom->GetAtomicNum();
    if(atomic_number == 1)
        return 1.20;
    if(atomic_number == 6)
        return 1.70;
    if(atomic_number == 7)
        return 1.55;
    if(atomic_number == 8)
        return 1.52;
    if(atomic_number == 9)
        return 1.47;
    if(atomic_number == 15)
        return 1.80;
    if(atomic_number == 16)
        return 1.80;
    if(atomic_number == 17)
        return 1.75;
    if(atomic_number == 35)
        return 1.85;
    if(atomic_number == 53)
        return 1.98;
    return 1.8;
}
// static method called by 'calcSASA'
static vector<unsigned> find_neighbor_indices(OBMol &mol, double probe, unsigned k)
{
    vector<unsigned> neighbor_indices;
    OBAtom *atom_k = mol.GetAtom(k);
    double radius = getRadius(atom_k) + probe + probe;
    for(unsigned i=1; i<=mol.NumAtoms(); ++i) {
        if(i == k)
            continue;
        OBAtom *atom_i = mol.GetAtom(i);
        double dist = pos_distance(atom_k->GetCoordinate(), atom_i->GetCoordinate());
        if(dist < radius+getRadius(atom_i))
            neighbor_indices.push_back(i);
    }
    return neighbor_indices;
}
// - solvent accessible surface area
// - modified based on 'asa.py' which can be downloaded from 'http://boscoh.com/protein/asapy'
// - ref
//   Shrake, A., and J. A. Rupley. "Environment and Exposure to Solvent
//   of Protein Atoms. Lysozyme and Insulin." JMB (1973) 79:351-371.
// 
vector<double> calcSASA(OBMol &mol)
{
    double probe = 1.4;
    int n_sphere_point = 960;
    vector<vector<double> > sphere_points = generate_sphere_points(n_sphere_point);
    double _const = 4.0 * PI / sphere_points.size();
    double *test_point = new double[3];
    vector<double> areas;
    for(unsigned i=1; i<=mol.NumAtoms(); ++i) {
        vector<unsigned> neighbor_indices = find_neighbor_indices(mol, probe, i);
        unsigned n_neighbor = neighbor_indices.size();
        unsigned j_closest_neighbor = 0;
        OBAtom *atom_i = mol.GetAtom(i);
        double radius = probe + getRadius(atom_i);
        unsigned *cycled_indices = new unsigned[n_neighbor];

        int n_accessible_point = 0;
        for(vector<vector<double> >::iterator point=sphere_points.begin();
                point!=sphere_points.end(); ++point) {
            bool is_accessible = true;
            test_point[0] = (*point)[0]*radius + atom_i->x();
            test_point[1] = (*point)[1]*radius + atom_i->y();
            test_point[2] = (*point)[2]*radius + atom_i->z();

            unsigned k=0;
            for(unsigned j=j_closest_neighbor; j<n_neighbor; ++j)
                cycled_indices[k++] = j;
            for(unsigned j=0; j<j_closest_neighbor; ++j)
                cycled_indices[k++] = j;

            for(k=0; k<n_neighbor; ++k) {
                unsigned j = cycled_indices[k];
                OBAtom *atom_j = mol.GetAtom(neighbor_indices[j]);
                double r = getRadius(atom_j) + probe;
                double diff_sq = pos_distance_sq(atom_j->GetCoordinate(), test_point);
                if(diff_sq < r*r) {
                    j_closest_neighbor = j;
                    is_accessible = false;
                    break;
                }
            }
            if(is_accessible)
                ++n_accessible_point;
        }
        delete[] cycled_indices;
        //std::cout << "_const=" << _const << ", n_accessible_point=" << n_accessible_point
        //    << ", radius=" << radius << endl;
        double area = _const * n_accessible_point * radius * radius;
        areas.push_back(area);
    }
    delete[] test_point;
    return areas;
}


/*
static string tmp_names[] = {"TMSA",
    "PPSA-1","PNSA-1","DPSA-1","FPSA-1","FNSA-1","WPSA-1","WNSA-1",
    "PPSA-2","PNSA-2","DPSA-2","FPSA-2","FNSA-2","WPSA-2","WNSA-2",
    "PPSA-3","PNSA-3","DPSA-3","FPSA-3","FNSA-3","WPSA-3","WNSA-3",
    "RPCG","RPCS","RNCG","RNCS",
    "HDSA","FHDSA","HASA","FHASA","HBSA","FHBSA",
    "HDCA","FHDCA","HACA","FHACA","HBCA","FHBCA"};
vector<string> cpsa_names(tmp_names,tmp_names+38);
*/

// result: TMSA, (PPSA, PNSA, DPSA, FPSA, FNSA, WPSA, WNSA) (1~3), RPCG, RPCS, RNCG, RNCS
// reference:
//   Stanton, D.T., Jurs, P.C. Analytical Chemistry. 1990, 62(21): 2323-2329.
//   DOI: 10.1021/ac00220a013
bool calcCPSA(OpenBabel::OBMol &mol, vector<double> &pc, double result[38])
{
    if(pc.empty())
        return false;
    
    memset(result,0,sizeof(double)*38);
    vector<double> area = calcSASA(mol);
    result[0] = accumulate(area.begin(),area.end(),0.);  // TMSA
    double pos_charge=0., neg_charge=0.;
    double max_pos_charge=-1., min_neg_charge=1.;
    double tot_pos_charge=0., tot_neg_charge=0.;
    int max_pos_i=-1, min_neg_i=-1;
    for(vector<double>::size_type i=0; i<pc.size(); ++i) {
        OBAtom *atom = mol.GetAtom(i+1);
        if(atom->MatchesSMARTS("[!$([#1,#6,F,Cl,Br,I,o,s,nx3,#7v5,#15v5,#16v4,#16v6])]")) {
            result[28] += area[i];            // HASA
            result[34] += (area[i] * fabs(pc[i]));  // HACA
        }
        if(atom->MatchesSMARTS("[#7,#8;!H0]")) {
            //result[26] += area[i];            // HDSA
            //result[32] += (area[i] * fabs(pc[i]));  // HDCA
            // including all attached hydrogens
            FOR_NBORS_OF_ATOM(nbor,atom) {
                if(!nbor->IsHydrogen())
                    continue;
                int indx = nbor->GetIdx();
                result[26] += area[indx-1];
                result[32] += (area[indx-1] * fabs(pc[indx-1]));
            }
        }

        if(pc[i] > 0.) {
            tot_pos_charge += pc[i];
            if(pc[i] > max_pos_charge) {
                max_pos_charge = pc[i];
                max_pos_i = i;
            }
            pos_charge += pc[i];
            result[1] += area[i];          // PPSA-1
            result[15] += (area[i]*pc[i]);  // PPSA-3
        }
        else if(pc[i] < 0.) {
            tot_neg_charge += pc[i];
            if(pc[i] < min_neg_charge) {
                min_neg_charge = pc[i];
                min_neg_i = i;
            }
            neg_charge += pc[i];
            result[2] += area[i];            // PNSA-1
            result[16] += (area[i]*pc[i]);   // PNSA-3
        }
    }
    result[8] = result[1] * pos_charge; // PPSA-2
    result[9] = result[2] * neg_charge; // PNSA-2
    result[22] = max_pos_charge / tot_pos_charge; // relative positive charge (RPCG)
    result[24] = min_neg_charge / tot_neg_charge; // relative negative charge (RNCG)
    result[23] = result[22] * area[max_pos_i]; // relative positive charged surface area (RPCS)
    result[25] = result[24] * area[min_neg_i]; // relative negative charged surface area (RNCS)
    for(int j=0; j<3; ++j) {
        result[7*j+3] = result[7*j+1] - result[7*j+2];    // DPSA-i
        result[7*j+4] = result[7*j+1] / result[0];        // FPSA-i
        result[7*j+5] = result[7*j+2] / result[0];        // FNSA-i
        result[7*j+6] = result[7*j+1] * result[0] / 1000; // WPSA-i
        result[7*j+7] = result[7*j+2] * result[0] / 1000; // WNSA-i
    }

    result[27] = result[26] / result[0];  // FHDSA
    result[29] = result[28] / result[0];  // FHASA
    result[30] = result[26] + result[28]; // HBSA
    result[31] = result[30] / result[0];  // FHBSA
    result[33] = result[32] / result[0];  // FHBCA
    result[35] = result[34] / result[0];  // FHACA
    result[36] = result[32] + result[34]; // HBCA
    result[37] = result[36] / result[0];  // FHACA
    
    return true;
}


// totally 26 SMARTSs for N
static string tpsa_n_smarts[] = {"[NH0X3](-*)(-*)-*","[NH0X2](-*)=*","[NH0X1]#*","[NH0X3](-*)(=*)=*",
    "[NH0X2](=*)#*","[NH0X3]1(-*)-*-*1","[NH1X3](-*)-*","[NH1X3]1-*-*1","[NH1X2]=*","[NH2X3]-*",
    "[NH0X4+](-*)(-*)(-*)-*","[NH0X3+](-*)(-*)=*","[NH0X2+](-*)#*","[NH1X4+](-*)(-*)-*",
    "[NH1X3+](-*)=*","[NH2X4+](-*)-*","[NH2X3+]=*","[NH3X4+]-*","[nH0X2](:*):*","[nH0X3](:*)(:*):*",
    "[nH0X3](-*)(:*):*","[nH0X3](=*)(:*):*","[nH1X3](:*):*","[nH0X3+](:*)(:*):*","[nH0X3+](-*)(:*):*",
    "[nH1X2+](:*):*"};
static double tpsa_n_value[] = {3.24,12.36,23.79,11.68,13.60,3.01,12.03,21.94,23.85,26.02,0.00,
    3.01,4.36,4.44,14.97,16.61,25.59,27.64,12.89,4.41,4.93,8.39,15.79,4.10,3.88,14.14};
// totally 6 SMARTSs for O
static string tpsa_o_smarts[] = {"[OH0X2](-*)-*","[OH0X2]1-*-*1","[OH0]=*","[OH]-*","[OH0-]-*",
    "[oH0X2](:*):*"};
static double tpsa_o_values[] = {9.23,12.53,17.07,20.23,23.06,13.14};
// totally 7 SMARTSs for S
static string tpsa_s_smarts[] = {"[SH0X2](-*)-*","[SH0X1]=*","[SH0X3](-*)(-*)=*",
    "[SH0X4](-*)(-*)(=*)=*","[SH1X2]-*","[sH0X2](:*):*","[sH0X3](=*)(:*):*"};
static double tpsa_s_values[] = {25.30,32.09,19.21,8.38,38.80,28.24,21.70};
// totally 4 SMARTSs for P
static string tpsa_p_smarts[] = {"[PH0X3](-*)(-*)-*","[PH0X2](-*)=*","[PH0X4](-*)(-*)(-*)=*",
    "[PH1X3](-*)(-*)=*"};
static double tpsa_p_values[] = {13.59,34.14,9.81,23.47};
// topologiial polar surface area
// ref: Ertl P. etc., J. Med. Chem. 2000, 43(20): 3714-3717
double calcTPSA(OBMol &mol)
{
    double value = 0.;
    FOR_ATOMS_OF_MOL(atom,mol) {
        unsigned atomic_number = atom->GetAtomicNum();
        if(atomic_number == 7) {
            for(int i=0; i<26; ++i) {
                if(atom->MatchesSMARTS(tpsa_n_smarts[i].c_str()))
                    value += tpsa_n_value[i];
            }
        }
        else if(atomic_number == 8) {
            for(int i=0; i<6; ++i) {
                if(atom->MatchesSMARTS(tpsa_o_smarts[i].c_str()))
                    value += tpsa_o_values[i];
            }
        }
        else if(atomic_number == 16) {
            for(int i=0; i<7; ++i) {
                if(atom->MatchesSMARTS(tpsa_s_smarts[i].c_str()))
                    value += tpsa_s_values[i];
            }
        }
        else if(atomic_number == 15) {
            for(int i=0; i<4; ++i) {
                if(atom->MatchesSMARTS(tpsa_p_smarts[i].c_str()))
                    value += tpsa_p_values[i];
            }
        }
        else
            ;
    }
    return value;
}


// Gravitation index
// value[0]: summation over all pairs of atoms
// value[1]: summation over all bonded pairs of atoms
void gravitation(OBMol &mol, double value[])
{
    double r2;
    // 1. summation over all pairs of atoms
    value[0] = 0.;
    for(unsigned i=1; i<mol.NumAtoms(); ++i) {
        for(unsigned j=i+1; j<=mol.NumAtoms(); ++j) {
            OBAtom *atom1 = mol.GetAtom(i);
            OBAtom *atom2 = mol.GetAtom(j);
            r2 = pow(atom1->GetX()-atom2->GetX(),2)+pow(atom1->GetY()-atom2->GetY(),2)+
                pow(atom1->GetZ()-atom2->GetZ(),2);
            value[0] += (atom1->GetAtomicMass() * atom2->GetAtomicMass() / r2);
        }
    }
    // 2. summation over all pairs of bonds
    value[1] = 0;
    FOR_BONDS_OF_MOL(bond,mol) {
        OBAtom *atom1 = bond->GetBeginAtom();
        OBAtom *atom2 = bond->GetEndAtom();
        r2 = pow(atom1->GetX()-atom2->GetX(),2)+pow(atom1->GetY()-atom2->GetY(),2)+
            pow(atom1->GetZ()-atom2->GetZ(),2);
        value[1] += (atom1->GetAtomicMass() * atom2->GetAtomicMass() / r2);
    }
}

