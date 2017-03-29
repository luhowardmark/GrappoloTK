// **************************************************************************************************
// Grappolo: A C++ library for parallel graph community detection
// Hao Lu, Ananth Kalyanaraman (hao.lu@wsu.edu, ananth@eecs.wsu.edu) Washington State University
// Mahantesh Halappanavar (hala@pnnl.gov) Pacific Northwest National Laboratory
//
// For citation, please cite the following paper:
// Lu, Hao, Mahantesh Halappanavar, and Ananth Kalyanaraman. 
// "Parallel heuristics for scalable community detection." Parallel Computing 47 (2015): 19-37.
//
// **************************************************************************************************
// Copyright (c) 2016. Washington State University ("WSU"). All Rights Reserved.
// Permission to use, copy, modify, and distribute this software and its documentation
// for educational, research, and not-for-profit purposes, without fee, is hereby
// granted, provided that the above copyright notice, this paragraph and the following
// two paragraphs appear in all copies, modifications, and distributions. For
// commercial licensing opportunities, please contact The Office of Commercialization,
// WSU, 280/286 Lighty, PB Box 641060, Pullman, WA 99164, (509) 335-5526,
// commercialization@wsu.edu<mailto:commercialization@wsu.edu>, https://commercialization.wsu.edu/

// IN NO EVENT SHALL WSU BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL,
// OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF
// THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF WSU HAS BEEN ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

// WSU SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE AND
// ACCOMPANYING DOCUMENTATION, IF ANY, PROVIDED HEREUNDER IS PROVIDED "AS IS". WSU HAS NO
// OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
// **************************************************************************************************


#ifndef __CLUSTERING__FUNCTIONS__
#define __CLUSTERING__FUNCTIONS__

#include "defs.h"

using namespace std;

void sumVertexDegree(edge* vtxInd, long* vtxPtr, double* vDegree, long NV, Comm* cInfo);

double calConstantForSecondTerm(double* vDegree, long NV);

void initCommAss(long* pastCommAss, long* currCommAss, long NV);

void initCommAssOpt(long* pastCommAss, long* currCommAss, long NV, 
		    mapElement* clusterLocalMap, long* vtxPtr, edge* vtxInd,
		    Comm* cInfo, double constant, double* vDegree );

double buildLocalMapCounter(long adj1, long adj2, map<long, long> &clusterLocalMap, 
						  vector<double> &Counter, edge* vtxInd, long* currCommAss, long me);

double buildLocalMapCounterNoMap(long v, mapElement* clusterLocalMap, long* vtxPtr, edge* vtxInd,
                               long* currCommAss, long &numUniqueClusters);

long max(map<long, long> &clusterLocalMap, vector<double> &Counter, 
		 double selfLoop, Comm* cInfo, double degree, long sc, double constant ) ;

long maxNoMap(long v, mapElement* clusterLocalMap, long* vtxPtr, double selfLoop, Comm* cInfo, double degree,
              long sc, double constant, long numUniqueClusters );

double computeMerkinMetric(long* C1, long N1, long* C2, long N2);
double computeVanDongenMetric(long* C1, long N1, long* C2, long N2);

#endif
