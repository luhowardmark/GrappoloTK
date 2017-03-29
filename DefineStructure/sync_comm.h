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

#ifndef __sync_comm__
#define __sync_comm__

#include "basic_util.h"
#include "utilityClusteringFunctions.h"

void runMultiPhaseSyncType(graph *G, long *C_orig, int syncType, long minGraphSize,
			double threshold, double C_threshold, int numThreads, int threadsOpt);

double parallelLouvainMethodFullSyncEarly(graph *G, long *C, int nThreads, double Lower,
				double thresh, double *totTime, int *numItr,int ytype, int freedom);
				
double parallelLouvainMethodFullSync(graph *G, long *C, int nThreads, double Lower,
				double thresh, double *totTime, int *numItr,int ytype, int freedom);
				
double parallelLouvianMethodEarlyTerminate(graph *G, long *C, int nThreads, double Lower,
				double thresh, double *totTime, int *numItr);
				
// Define in fullSyncUtility.cpp
double buildAndLockLocalMapCounter(long v, mapElement* clusterLocalMap, long* vtxPtr, edge* vtxInd,
                               long* currCommAss, long &numUniqueClusters, omp_lock_t* vlocks, omp_lock_t* clocks, int ytype, double& eix, int freedom);

void maxAndFree(long v, mapElement* clusterLocalMap, long* vtxPtr, edge* vtxInd, double selfLoop, Comm* cInfo, long* CA, 
							double constant, long numUniqueClusters, omp_lock_t* vlocks, omp_lock_t* clocks, int ytype, double eix, double* vDegree);

#endif
