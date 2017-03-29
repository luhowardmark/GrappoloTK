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

#include "defs.h"
#include "sync_comm.h"
#include <algorithm>
using namespace std;
bool byVertexId(edge v1, edge v2){
	return (v1.tail<v2.tail);
}

bool byCommId(mapElement c1,mapElement c2){
	return (c1.cid<c2.cid);
}


//Build the local-map data structure using vectors
double buildAndLockLocalMapCounter(long v, mapElement* clusterLocalMap, long* vtxPtr, edge* vtxInd,
                               long* currCommAss, long &numUniqueClusters, omp_lock_t* vlocks, omp_lock_t* clocks, int ytype, double& eix, int freedom) {
  double selfLoop = 0;
	long adj1  = vtxPtr[v];
	long adj2  = vtxPtr[v+1];


	// Has to be one over to act like array.end()
	std::sort(&vtxInd[adj1],&vtxInd[adj2],byVertexId);

	/*********** Calculate eii ***************/
	// Lock all neighbors to make sure no move is performed // This Lock is to protect Data: Allow to have error
	for(long j=adj1; j<adj2; j++){
   // long nn = adj2-adj1;
   // long fracFree = (nn*freedom/10);
   // if( rand()%nn >= fracFree)
  		omp_set_lock(&vlocks[vtxInd[j].tail]);
	}
	

	// Aggregate the neighbors and lock their Community
	long sPosition = vtxPtr[v]+v; //Starting position of local map for v
	long storedAlready = 0;

	for(long j=adj1; j<adj2; j++) {
		if(vtxInd[j].tail == v) {	// SelfLoop need to be recorded
			selfLoop += vtxInd[j].weight;
		}
		
		bool storedAlready = false; //Initialize to zero
		for(long k=0; k<numUniqueClusters; k++) { //Check if it already exists
			if(currCommAss[vtxInd[j].tail] ==  clusterLocalMap[sPosition+k].cid) {
				storedAlready = true;
				clusterLocalMap[sPosition + k].Counter += vtxInd[j].weight; //Increment the counter with weight
				break;
			}
		}
		if( storedAlready == false ) {	//Does not exist, add to the map
			clusterLocalMap[sPosition + numUniqueClusters].cid     = currCommAss[vtxInd[j].tail];
			clusterLocalMap[sPosition + numUniqueClusters].Counter = vtxInd[j].weight; //Initialize the count
			numUniqueClusters++;
		}
	}//End of for(j)*/
	eix = clusterLocalMap[sPosition].Counter - selfLoop;

	if(ytype == 1){
		// Locking the community information for all neighbors // This lock is to protect Data: Allow to have error
		std::sort(&clusterLocalMap[sPosition],&clusterLocalMap[sPosition+numUniqueClusters],byCommId);  
		for(long j=sPosition; j<sPosition + numUniqueClusters; j++){
			omp_set_lock(&clocks[clusterLocalMap[j].cid]);
		}
	}

	return selfLoop;
}//End of buildLocalMapCounter()


void maxAndFree(long v, mapElement* clusterLocalMap, long* vtxPtr, edge* vtxInd, double selfLoop, Comm* cInfo, long* CA, 
							double constant, long numUniqueClusters, omp_lock_t* vlocks, omp_lock_t* clocks, int ytype, double eix, double* vDegree) {
                                                                                
	long maxIndex = CA[v];	//Assign the initial value as the current community
	long sc = CA[v];		
	double curGain = 0;
	double maxGain = 0;
	long sPosition = vtxPtr[v]+v; //Starting position of local map for v
	double degree = vDegree[v];
	double ax  = cInfo[sc].degree - degree;
	double eiy = 0;
	double ay  = 0;

		
	/*********** Calculate DeltaQ using aii ***************/    
	for(long k=0; k<numUniqueClusters; k++) {
		if(sc != clusterLocalMap[sPosition + k].cid) {
			ay = cInfo[clusterLocalMap[sPosition + k].cid].degree; // degree of cluster y
			eiy = clusterLocalMap[sPosition + k].Counter; 	//Total edges incident on cluster y
			curGain = 2*(eiy - eix) - 2*degree*(ay - ax)*constant;
			if( (curGain > maxGain) ||
					((curGain==maxGain) && (curGain != 0) && (clusterLocalMap[sPosition + k].cid < maxIndex)) ) {
				maxGain  = curGain;
				maxIndex = clusterLocalMap[sPosition + k].cid;
			}
		}
	}//End of for()

	if (sc != maxIndex){
		if(ytype == 1){
			CA[v] = maxIndex;
			cInfo[maxIndex].degree += vDegree[v];
			cInfo[maxIndex].size += 1;
			cInfo[sc].degree -= vDegree[v];
			cInfo[sc].size -=1;
		
		}else{
			CA[v] = maxIndex;
			#pragma omp atomic update
			cInfo[maxIndex].degree += vDegree[v];
			#pragma omp atomic update
			cInfo[maxIndex].size += 1;
			#pragma omp atomic update
			cInfo[sc].degree -= vDegree[v];
			#pragma omp atomic update
			cInfo[sc].size -=1;
		}
	}

	if(ytype == 1){
		// unLock all neighbors community 	
		for(long j=sPosition; j<sPosition + numUniqueClusters; j++){
				omp_unset_lock(&clocks[clusterLocalMap[j].cid]);
		}
	}

	// Free Neighbor
	long adj1  = vtxPtr[v];
	long adj2  = vtxPtr[v+1];
	/*********** Calculate eii ***************/
	// unLock all neighbors vertex
	for(long j=adj1; j<adj2; j++){
		omp_unset_lock(&vlocks[vtxInd[j].tail]);
	}


	
}//End maxNoMap()
