// **************************************************************************************************
// GrappoloTK: A C++ library for parallel graph coloring
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
#include "basic_comm.h"
#include "color_comm.h"
using namespace std;
//WARNING: This will overwrite the original graph data structure to 
//         minimize memory footprint
// Return: C_orig will hold the cluster ids for vertices in the original graph
//         Assume C_orig is initialized appropriately
//WARNING: Graph G will be destroyed at the end of this routine
void runMultiPhaseColoring(graph *G, long *C_orig, int coloring, long minGraphSize,
			double threshold, double C_threshold, int numThreads, int threadsOpt)
{
  double totTimeClustering=0, totTimeBuildingPhase=0, totTimeColoring=0, tmpTime;
  int tmpItr=0, totItr = 0;  
  long NV = G->numVertices;

  //long minGraphSize = 100000; //Need at least 100,000 vertices to turn coloring on

  int *colors;
  int numColors = 0;
  
	// Coloring Steps
	if(coloring >= 1) {
	  colors = (int *) malloc (G->numVertices * sizeof(int)); assert (colors != 0);
    #pragma omp parallel for
	  for (long i=0; i<G->numVertices; i++) {
		  colors[i] = -1;
	  }
	  numColors = algoDistanceOneVertexColoringOpt(G, colors, numThreads, &tmpTime)+1;
	  totTimeColoring += tmpTime;
	  //printf("Number of colors used: %d\n", numColors);
		
  }
	if(coloring == 2){
		vBaseRedistribution(G, colors, numColors, 0);
	}
	
  /* Step 3: Find communities */
  double prevMod = -1;
  double currMod = -1;
  long phase = 1;

  graph *Gnew; //To build new hierarchical graphs
  long numClusters;
  long *C = (long *) malloc (NV * sizeof(long));
  assert(C != 0);
  #pragma omp parallel for
  for (long i=0; i<NV; i++) {
  	C[i] = -1;
  }	
  
	bool nonColor = false; //Make sure that at least one phase with lower threshold runs
  while(1){
    printf("===============================\n");
	  printf("Phase %ld\n", phase);
    printf("===============================\n");
   	prevMod = currMod;
	  //Compute clusters
	  if((G->numVertices > minGraphSize)&&(nonColor == false)) {
		  // No Map is not constructed yet
			// currMod = algoLouvainWithDistOneColoringNoMap(G, C, numThreads, colors, numColors, currMod, C_threshold, &tmpTime, &tmpItr);
      currMod = algoLouvainWithDistOneColoring(G, C, numThreads, colors, numColors, currMod, C_threshold, &tmpTime, &tmpItr);
		  totTimeClustering += tmpTime;
      totItr += tmpItr;
			if(phase == 1)
				nonColor = true;
	  }else {
      parallelLouvianMethodNoMap(G, C, numThreads, currMod, threshold, &tmpTime, &tmpItr);
      totTimeClustering += tmpTime;
      totItr += tmpItr;
	  } 
  
    //Renumber the clusters contiguiously
  	numClusters = renumberClustersContiguously(C, G->numVertices);
  	printf("Number of unique clusters: %ld\n", numClusters);
  
    //printf("About to update C_orig\n");
	  //Keep track of clusters in C_orig
	  if(phase == 1) {
      #pragma omp parallel for
	    for (long i=0; i<NV; i++) {
	  	  C_orig[i] = C[i]; //After the first phase
	    } 	
	  } else {
      #pragma omp parallel for
	    for (long i=0; i<NV; i++) {
        assert(C_orig[i] < G->numVertices);
        if (C_orig[i] >=0)
	  		  C_orig[i] = C[C_orig[i]]; //Each cluster in a previous phase becomes a vertex
	    }	
	  }
    printf("Done updating C_orig\n");
	  //Break if too many phases or iterations
	  if((phase > 200)||(totItr > 10000)) {
	  	break;
	  }
    
    //Check for modularity gain and build the graph for next phase
	  //In case coloring is used, make sure the non-coloring routine is run at least once

    if( (currMod - prevMod) > threshold ) {
		  Gnew = (graph *) malloc (sizeof(graph)); assert(Gnew != 0);
		  tmpTime =  buildNextLevelGraphOpt(G, Gnew, C, numClusters, numThreads);
		  totTimeBuildingPhase += tmpTime;
		  //Free up the previous graph		
		  free(G->edgeListPtrs);	
		  free(G->edgeList);
		  free(G);
		  G = Gnew; //Swap the pointers
      G->edgeListPtrs = Gnew->edgeListPtrs;
      G->edgeList = Gnew->edgeList;
		  
      //Free up the previous cluster & create new one of a different size
		  free(C);
		  C = (long *) malloc (numClusters * sizeof(long)); assert(C != 0);
      
      #pragma omp parallel for
		  for (long i=0; i<numClusters; i++) {
			  C[i] = -1;
		  }
		  phase++; //Increment phase number
		  //If coloring is enabled & graph is of minimum size, recolor the new graph
		  /*if((coloring == 1)&&(G->numVertices > minGraphSize)&&(nonColor = false)){
        #pragma omp parallel for
			  for (long i=0; i<G->numVertices; i++){
			  	colors[i] = -1;
			  }
			  numColors = algoDistanceOneVertexColoringOpt(G, colors, numThreads, &tmpTime)+1;
			  totTimeColoring += tmpTime;
		  }*/
	  } else { //To force another phase with coloring again
		  if ( nonColor ) {
		  	nonColor = false;
				if(coloring == 1){
					#pragma omp parallel for
					for (long i=0; i<G->numVertices; i++){
						colors[i] = -1;
					}
					numColors = algoDistanceOneVertexColoringOpt(G, colors, numThreads, &tmpTime)+1;
					totTimeColoring += tmpTime;
				}else if(coloring == 2){
					// Space for balacing
				}
		  }
		  else {
		    break; //Modularity gain is not enough. Exit.
		  }
	  } 	
  } //End of while(1)
 
  printf("********************************************\n"); 
  printf("*********    Compact Summary   *************\n");
  printf("********************************************\n");
  printf("Number of threads              : %ld\n", numThreads);
  printf("Total number of phases         : %ld\n", phase);
  printf("Total number of iterations     : %ld\n", totItr);
  printf("Final number of clusters       : %ld\n", numClusters);
  printf("Final modularity               : %lf\n", prevMod);
  printf("Total time for clustering      : %lf\n", totTimeClustering);
  printf("Total time for building phases : %lf\n", totTimeBuildingPhase);
  if(coloring == 1) {
     printf("Total time for coloring        : %lf\n", totTimeColoring);
  }
  printf("********************************************\n");
  printf("TOTAL TIME                     : %lf\n", (totTimeClustering+totTimeBuildingPhase+totTimeColoring) );
  printf("********************************************\n");

  //Clean up:
  free(C);
  if(G != 0) {
    free(G->edgeListPtrs);
    free(G->edgeList);
    free(G);
  }

  if(coloring > 0) {
    if(colors != 0) free(colors);
  }

}//End of runMultiPhaseLouvainAlgorithm()
