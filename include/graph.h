/*=============================================================================
#     FileName: graph.h
#         Desc: some generic graph methods
#       Author: jlpeng
#        Email: jlpeng1201@gmail.com
#     HomePage: 
#      Version: 0.0.1
#      Created: 2012-12-28 15:52:13
#   LastChange: 2014-01-02 21:01:53
#      History:
=============================================================================*/

#ifndef  GRAPH_H
#define  GRAPH_H
#include <vector>
#include <openbabel/mol.h>
#include <openbabel/bond.h>
#include "config.h"
// to get length of shortest path from 'startId' to other atoms of 'mol'
// "exclude hydrogen"
std::vector<double> DESCRIPTOR_API shortestDist(OpenBabel::OBMol &mol, unsigned startId);
// to get all possible and nonredundant paths of length 'maxLen'
// "exclude hydrogen"
std::vector<std::vector<unsigned> > DESCRIPTOR_API pathOfLenN(OpenBabel::OBMol &mol, int maxLen);
// to calculate shortest disctance from each atom to any other atoms
// Floyd-Warshall algorithm
// reference:  http://en.wikipedia.org/wiki/Floyd_algorithm
std::vector<std::vector<int> > DESCRIPTOR_API floyd(OpenBabel::OBMol &mol);


// to calculate weighted shortest distance between 'startId' to any other heavy atoms
// getWeight: user defined function to assign weight for specific bond type
//            if NULL, the following weight values will be applied:
//            [single: 1.; double: 1./2; triple: 1./3; aromatic: 2./3]
std::vector<double> DESCRIPTOR_API weightedShortestDist(OpenBabel::OBMol &mol, 
        unsigned startId, double (*getWeight)(const OpenBabel::OBBond *bond)=NULL);

#endif   /* ----- #ifndef GRAPH_H  ----- */
