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
#include "utilityClusteringFunctions.h"
#include "sync_comm.h"

using namespace std;

double parallelLouvainMethodFullSyncEarly(graph *G, long *C, int nThreads, double Lower,
				double thresh, double *totTime, int *numItr,int ytype, int freedom) {
#ifdef PRINT_DETAILED_STATS_  
  printf("Within parallelLouvianMethodNoMap()\n");
#endif
  if (nThreads < 1)
    omp_set_num_threads(1);
  else
    omp_set_num_threads(nThreads);
  int nT;
#pragma omp parallel
  {
    nT = omp_get_num_threads();
  }
#ifdef PRINT_DETAILED_STATS_
  printf("Actual number of threads: %d (requested: %d)\n", nT, nThreads);
#endif
  double time1, time2, time3, time4; //For timing purposes  
  double total = 0, totItr = 0;
  
  long    NV        = G->numVertices;
  long    NS        = G->sVertices;      
  long    NE        = G->numEdges;
  long    *vtxPtr   = G->edgeListPtrs;
  edge    *vtxInd   = G->edgeList;
 
  /* Variables for computing modularity */
  long totalEdgeWeightTwice;
  double constantForSecondTerm;
  double prevMod=-1;
  double currMod=-1;
  double thresMod = thresh; //Input parameter
  int numItrs = 0;
  
  /********************** Initialization **************************/
  time1 = omp_get_wtime();
  //Store the degree of all vertices
  double* vDegree = (double *) malloc (NV * sizeof(double)); assert(vDegree != 0);
  //Community info. (ai and size)
  Comm *cInfo = (Comm *) malloc (NV * sizeof(Comm)); assert(cInfo != 0);
  omp_lock_t* vlocks = (omp_lock_t*) malloc (NV*sizeof(*vlocks));
  omp_lock_t* clocks = (omp_lock_t*) malloc (NV*sizeof(*clocks));

  //use for Modularity calculation (eii)
  double* clusterWeightInternal = (double*) malloc (NV*sizeof(double)); assert(clusterWeightInternal != 0);

  sumVertexDegree(vtxInd, vtxPtr, vDegree, NV , cInfo);	// Sum up the vertex degree
  
  /*** Compute the total edge weight (2m) and 1/2m ***/
  constantForSecondTerm = calConstantForSecondTerm(vDegree, NV); // 1 over sum of the degree
     
  //Vectors used in place of maps: Total size = |V|+2*|E| -- The |V| part takes care of self loop
  mapElement* clusterLocalMap = (mapElement *) malloc ((NV + 2*NE) * sizeof(mapElement)); assert(clusterLocalMap != 0);

  
  //Store previous iteration's community assignment
  long* pastCommAss = (long *) malloc (NV * sizeof(long)); assert(pastCommAss != 0);
  //Store current community assignment
  long* currCommAss = (long *) malloc (NV * sizeof(long)); assert(currCommAss != 0);  
  //Store the target of community assignment  
    
  //Initialize each vertex to its own cluster
  initCommAss(C, C, NV); 
  initCommAss(pastCommAss, currCommAss, NV); 
  
  
  // Store the termination node
  bool* verT = (bool *) malloc (NV * sizeof(bool)); assert(verT != 0);	
  #pragma omp parallel for
  for (long i=0; i<NV; i++) {
    verT[i] = false;
  }
  long termNodes = 0;
  
  time2 = omp_get_wtime();
  printf("Time to initialize: %3.3lf\n", time2-time1);
	

  // Set up locks for full sync
  #pragma omp parallel for
  for (long i=0; i<NV; i++) {
    omp_init_lock(&vlocks[i]);
    omp_init_lock(&clocks[i]);
  }


#ifdef PRINT_DETAILED_STATS_
  printf("========================================================================================================\n");
  printf("Itr      E_xx            A_x2           Curr-Mod         Time-1(s)       Time-2(s)        T/Itr(s)\n");
  printf("========================================================================================================\n");
#endif
#ifdef PRINT_TERSE_STATS_
  printf("=====================================================\n");
  printf("Itr      Curr-Mod         T/Itr(s)      T-Cumulative\n");
  printf("=====================================================\n");
#endif
  //Start maximizing modularity
  while(true) {
    numItrs++;    
    time1 = omp_get_wtime();
    /* Re-initialize datastructures */
    
	long totalEdgeTravel= 0;
	long totalUniqueComm = 0;
	
	#pragma omp parallel for reduction(+:totalEdgeTravel), reduction(+:totalUniqueComm)
    for (long i=0; i<NV; i++) {
		if(verT[i])
			continue;
      long adj1 = vtxPtr[i];
      long adj2 = vtxPtr[i+1];
      long selfLoop = 0;
	  totalEdgeTravel += (adj2-adj1);
      long numUniqueClusters = 0;
	    //Add v's current cluster:
	    if(adj1 != adj2){
			//Add the current cluster of i to the local map
			long sPosition = vtxPtr[i]+i; //Starting position of local map for i
			double eix;        
			clusterLocalMap[sPosition].Counter = 0;          //Initialize the counter to ZERO (no edges incident yet)
			clusterLocalMap[sPosition].cid = C[i]; //Initialize with current community
			numUniqueClusters++; //Added the first entry
          
			//Find unique cluster ids and #of edges incident (eicj) to them
			selfLoop = buildAndLockLocalMapCounter(i, clusterLocalMap, vtxPtr, vtxInd, C, numUniqueClusters, vlocks, clocks, ytype, eix, freedom);
			// Update delta Q calculation
			//Calculate the max
			maxAndFree(i, clusterLocalMap, vtxPtr, vtxInd, selfLoop, cInfo, C, constantForSecondTerm, numUniqueClusters, vlocks, clocks, ytype, eix, vDegree);
            //assert((targetCommAss[i] >= 0)&&(targetCommAss[i] < NV));
			  
			if(numItrs > 2 && C[i] == currCommAss[i] && pastCommAss[i]==currCommAss[i]){
				//Swaping!!!
				verT[i] = true;
				termNodes++;
			}
			else{
				pastCommAss[i] = currCommAss[i];
				currCommAss[i] = C[i];
			}
      } else {

      }
	  totalUniqueComm += numUniqueClusters;
    }//End of for(i)
    time2 = omp_get_wtime();
    
		time3 = omp_get_wtime();    
    double e_xx = 0;
    double a2_x = 0;	


		// Calculate Modularity
		#pragma omp parallel for  //Parallelize on each vertex
		for (long i =0; i<NV;i++){
			clusterWeightInternal[i] = 0;
		}
		#pragma omp parallel for  //Parallelize on each vertex
		for (long i=0; i<NV; i++) {
			long adj1 = vtxPtr[i];
			long adj2 = vtxPtr[i+1];
			for(long j=adj1; j<adj2; j++) {
				if(C[vtxInd[j].tail] == C[i]){
					clusterWeightInternal[i] += vtxInd[j].weight;
				}
			}
		}		
		#pragma omp parallel for reduction(+:e_xx) reduction(+:a2_x)
    for (long i=0; i<NV; i++) {
      e_xx += clusterWeightInternal[i];
      a2_x += (cInfo[i].degree)*(cInfo[i].degree);
    }
    time4 = omp_get_wtime();

    currMod = (e_xx*(double)constantForSecondTerm) - (a2_x*(double)constantForSecondTerm*(double)constantForSecondTerm);
    totItr = (time2-time1) + (time4-time3);
    total += totItr;

#ifdef PRINT_DETAILED_STATS_
    //printf("%d \t %g \t %g \t %lf \t %3.3lf \t %3.3lf  \t %3.3lf\n",numItrs, e_xx, a2_x, currMod, (time2-time1), (time4-time3), totItr );
	printf("%d %d %d %d %d %3.5lf\n",numItrs, NV, termNodes, totalEdgeTravel, totalUniqueComm, currMod);
#endif
#ifdef PRINT_TERSE_STATS_
   printf("%d \t %lf \t %3.3lf  \t %3.3lf\n",numItrs, currMod, totItr, total);
#endif

    //Break if modularity gain is not sufficient
    if((currMod - prevMod) < thresMod) {
      break;
    }
    prevMod = currMod;
  }//End of while(true)
  *totTime = total; //Return back the total time for clustering
  *numItr  = numItrs;

#ifdef PRINT_DETAILED_STATS_
  printf("========================================================================================================\n");
  printf("Total time for %d iterations is: %lf\n",numItrs, total);  
  printf("========================================================================================================\n");
#endif  
#ifdef PRINT_TERSE_STATS_
  printf("========================================================================================================\n");
  printf("Total time for %d iterations is: %lf\n",numItrs, total);  
  printf("========================================================================================================\n");
#endif

  //Cleanup
  free(vDegree);
  free(cInfo);
  free(clusterWeightInternal);
  free(clusterLocalMap);

  return currMod;
}