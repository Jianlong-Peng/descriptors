/*=============================================================================
#     FileName: graph.cpp
#         Desc: 
#       Author: jlpeng
#        Email: jlpeng1201@gmail.com
#     HomePage: 
#      Version: 0.0.1
#      Created: 2012-12-28 15:53:06
#   LastChange: 2014-02-14 11:11:32
#      History:
=============================================================================*/
//#include <iterator>
#include <iostream>
#include <vector>
#include <bitset>
#include <queue>
#include <sstream>
#include <string>
#include <cassert>
#include <climits>
#include "graph.h"
#include "tools.h"
#include <openbabel/atom.h>
#include <openbabel/mol.h>
using std::queue;
using std::bitset;
using std::cerr;
using std::endl;
using std::vector;
using std::ostringstream;
using std::string;
using namespace OpenBabel;

vector<double> shortestDist(OBMol &mol, unsigned startId)
{
    if(mol.NumAtoms() >300) {
        cerr << "(graph.cpp::shortestDist)Warning: there are more than 300 atoms!";
        return vector<double>();
    }

    //pathLength[i]: shortest path length from `startId` to `i+1`
    vector<double> pathLength(mol.NumAtoms(),0.);
    bitset<300> visited;
    //bool *visited = (bool*)malloc(sizeof(bool)*(mol.NumAtoms()));
    queue<unsigned int> toBeVisit;

    //push atom `startId` to the queue
    toBeVisit.push(startId);
    visited.set(startId-1);

    //BFS
    while(!toBeVisit.empty())
    {
        unsigned int currentId = toBeVisit.front();

        //to visit all neighbors that are not hydrogen and haven't been visited
        FOR_NBORS_OF_ATOM(neighbor, mol.GetAtom(currentId))
        {
            unsigned neighId = neighbor->GetIdx();
            //if not visited, then set pathLength[neighId-1]
            if(!visited.test(neighId-1) && !neighbor->IsHydrogen()) {
                pathLength[neighId-1] = pathLength[currentId-1]+1.;
                visited.set(neighId-1);
                toBeVisit.push(neighId);
            }
        }
        //delete the currently visited atom from `toBeVisit`
        toBeVisit.pop();
    }

    return pathLength;
}


// called by 'pathOfLenN'
// no duplicated vertices ?????????????
static void pathDFS(OBAtom &atom, vector<vector<unsigned> > &path, vector<unsigned> tmp_path,
        int maxDepth, int currDepth)
{
    // atom already visited in the current path
    if(find(tmp_path.begin(),tmp_path.end(),atom.GetIdx()) != tmp_path.end())
        return;

    tmp_path.push_back(atom.GetIdx());
    if(currDepth == maxDepth) {
        path.push_back(tmp_path);
        return;
    }
    FOR_NBORS_OF_ATOM(nbor,atom) {
        if(nbor->IsHydrogen())
            continue;
        pathDFS(*nbor,path,tmp_path,maxDepth,currDepth+1);
    }
}

static void refinePath(vector<unsigned> &path)
{
    vector<unsigned>::size_type n = path.size();
    if(path[0] <= path[n-1])
        return;
    for(unsigned i=0; i<n/2; ++i) {
        unsigned temp = path[i];
        path[i] = path[n-1-i];
        path[n-1-i] = temp;
    }
}
// 
vector<vector<unsigned> > pathOfLenN(OBMol &mol, int maxDepth)
{
    vector<vector<unsigned> > path;
    vector<string> unique_path;
    
    FOR_ATOMS_OF_MOL(atom,mol) {
        if(atom->IsHydrogen())
            continue;
        //std::cout << "atom " << atom->GetIdx() << endl;
        vector<vector<unsigned> > spath;
        vector<unsigned> tmp_path;
        pathDFS(*atom,spath,tmp_path,maxDepth,0);
        for(vector<vector<unsigned> >::iterator iter=spath.begin(); iter!=spath.end(); ++iter) {
            if(static_cast<int>(iter->size()) != maxDepth+1)
                continue;
            //copy(iter->begin(),iter->end(),std::ostream_iterator<unsigned>(std::cout," "));
            //std::cout << endl;
            //sort(iter->begin(),iter->end());
            refinePath(*iter);
            string path_str = joinBy(iter->begin(),iter->end(),"_");
            if(find(unique_path.begin(),unique_path.end(),path_str) == unique_path.end()) {
                path.push_back(*iter);
                unique_path.push_back(path_str);
            }
        }
    }
    return path;
}

// reference:  http://en.wikipedia.org/wiki/Floyd_algorithm
vector<vector<int> > floyd(OBMol &mol)
{
    // initialization
    vector<vector<int> > dis;
    vector<bool> searched(mol.NumAtoms(), false);
    for(unsigned i=0; i<mol.NumAtoms(); ++i)
        dis.push_back(vector<int>(mol.NumAtoms(), INT_MAX));
    for(vector<bool>::size_type i=0; i<searched.size(); ++i) {
        dis[i][i] = 0;
        searched[i] = true;
        FOR_NBORS_OF_ATOM(nbor, mol.GetAtom(i+1)) {
            unsigned j = nbor->GetIdx()-1;
            if(searched[j])
                continue;
            dis[i][j] = 1;
            dis[i][j] = 1;
        }
    }
    //Floyd-Warshall algorithm
    unsigned natoms = mol.NumAtoms();
    for(unsigned k=0; k<natoms; ++k) {
        for(unsigned i=0; i<natoms; ++i) {
            for(unsigned j=0; j<natoms; ++j) {
                if(dis[i][k]==INT_MAX || dis[k][j]==INT_MAX)
                    continue;
                if(dis[i][k]+dis[k][j] < dis[i][j])
                    dis[i][j] = dis[i][k] + dis[k][j];
            }
        }
    }
    return dis;
}



vector<double> weightedShortestDist(OBMol &mol, unsigned startId,
        double (*getWeight)(const OBBond *bond))
{
    if(mol.NumAtoms() > 200) {
        std::cerr << "(graph.cpp::weightedShortestDist)Warning: there are more than 200 atoms!";
        return vector<double>();
    }

    //pathLength[i]: shortest path length from `startId` to `i+1`
    vector<double> pathLength(mol.NumAtoms(),1e8);
    pathLength[startId-1] = 0.;
    bitset<200> visited;
    //bool *visited = (bool*)malloc(sizeof(bool)*(mol.NumAtoms()));
    queue<unsigned int> toBeVisit;

    //push atom `startId` to the queue
    toBeVisit.push(startId);
    visited.set(startId-1);

    //BFS
    while(!toBeVisit.empty())
    {
        unsigned int currentId = toBeVisit.front();
        //unsigned current_atom_num = mol.GetAtom(currentId)->GetAtomicNum();

        //to visit all neighbors that are not hydrogen and haven't been visited
        FOR_NBORS_OF_ATOM(neighbor, mol.GetAtom(currentId))
        {
            unsigned neighId = neighbor->GetIdx();
            //unsigned neigh_atom_num = mol.GetAtom(neighId)->GetAtomicNum();
            if(neighId!=startId && !neighbor->IsHydrogen()) {
                OBBond *bond = mol.GetBond(mol.GetAtom(currentId),mol.GetAtom(neighId));
                double weight;
                if(getWeight == NULL) {
                    if(bond->IsAromatic())
                        weight = 2./3;
                    else if(bond->IsTriple())
                        weight = 1./3;
                    else if(bond->IsDouble())
                        weight = 1./2;
                    else
                        weight = 1.;
                }
                else
                    weight = getWeight(bond);
                if(pathLength[currentId-1]+weight < pathLength[neighId-1]) {
                    pathLength[neighId-1] = pathLength[currentId-1]+weight;
                    visited.set(neighId-1);
                    toBeVisit.push(neighId);
                }
            }
        }
        //delete the currently visited atom from `toBeVisit`
        toBeVisit.pop();
    }
    FOR_ATOMS_OF_MOL(atom,mol) {
        if(atom->IsHydrogen())
            pathLength[atom->GetIdx()-1] = 0.;
    }
    //pathLength[startId-1] = 1. - 6. / mol.GetAtom(startId)->GetAtomicNum();

    return pathLength;
}


