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
#include "basic_comm.h"
using namespace std;

long vertexFollowing(graph *G, long *C)
{
	long    NV        = G->numVertices;
	long    *vtxPtr   = G->edgeListPtrs;
	edge    *vtxInd   = G->edgeList;
	long numNode = 0;
	double time1 = omp_get_wtime();
// Initialize the Communities
#pragma omp parallel for  //Parallelize on the outer most loop
  	for (long i=0; i<NV; i++) {
		C[i] = i; //Initialize each vertex to its own cluster
	}

// Remove Isolated and degree-one vertices
#pragma omp parallel for
    	for (long i=0; i<NV; i++) {
		long adj1 = vtxPtr[i];
		long adj2 = vtxPtr[i+1];
		if(adj1 == adj2) {	// Isolated vertex
			__sync_fetch_and_add(&numNode, 1);
			C[i] = -1;
		} else {
			if( (adj2-adj1) == 1 ) { //Degree one
			    //Check if the tail has degree greater than one:
			    long tail = vtxInd[adj1].tail;
			    long adj11 = vtxPtr[tail];
			    long adj12 = vtxPtr[tail+1];
                            if( ((adj12-adj11) > 1)||(i > tail) ) { //Degree of tail greater than one
				__sync_fetch_and_add(&numNode, 1);
				C[i] = tail;
			    } //else don't do anything
			}//End of if(degree one)
		}//End of else
	}//End of for(i)

        time1 = omp_get_wtime() - time1;
#ifdef PRINT_DETAILED_STATS_
        printf("Time to determine number of vertices (numNode) to fix: %lf\n", time1);	
#endif
	return numNode; //These are nodes that need to be removed
}//End of vertexFollowing()

//WARNING: Will assume that the cluster id have been renumbered contiguously
//Return the total time for building the next level of graph
//This will not add any self-loops
double buildNewGraphVF(graph *Gin, graph *Gout, long *C, long numUniqueClusters) {
  int nT;
#pragma omp parallel
  {
    nT = omp_get_num_threads();
  }
#ifdef PRINT_DETAILED_STATS_
  printf("Within buildNewGraphVF(): # of unique clusters= %ld\n",numUniqueClusters);
  printf("Actual number of threads: %d \n", nT);
#endif

  double time1, time2, TotTime=0; //For timing purposes  
  double total = 0, totItr = 0;  
  //Pointers into the input graph structure:
  long    NV_in        = Gin->numVertices;
  long    NE_in        = Gin->numEdges;
  long    *vtxPtrIn    = Gin->edgeListPtrs;
  edge    *vtxIndIn    = Gin->edgeList;
  
  time1 = omp_get_wtime();
  // Pointers into the output graph structure
  long NV_out = numUniqueClusters;
  long NE_self = 0; //Not all vertices get self-loops
  long NE_out = 0;  //Cross edges
  long *vtxPtrOut = (long *) malloc ((NV_out+1)*sizeof(long));
  assert(vtxPtrOut != 0);
  vtxPtrOut[0] = 0; //First location is always a zero
  /* Step 1 : Regroup the node into cluster node */
  map<long,long>** cluPtrIn = (map<long,long>**) malloc(numUniqueClusters*sizeof(map<long,long>*));
  assert(cluPtrIn != 0);	

#pragma omp parallel for
  for (long i=0; i<numUniqueClusters; i++) {
	cluPtrIn[i] = new map<long,long>();
	//Do not add self-loops
        //(*(cluPtrIn[i]))[i] = 0; //Add for a self loop with zero weight
  }
#pragma omp parallel for
  for (long i=1; i<=NV_out; i++)
	vtxPtrOut[i] = 0; 

  //Create an array of locks for each cluster
  omp_lock_t *nlocks = (omp_lock_t *) malloc (numUniqueClusters * sizeof(omp_lock_t));
  assert(nlocks != 0);
#pragma omp parallel for
  for (long i=0; i<numUniqueClusters; i++) {
    omp_init_lock(&nlocks[i]); //Initialize locks
  }
  time2 = omp_get_wtime();
  TotTime += (time2-time1);
#ifdef PRINT_DETAILED_STATS_
  printf("Time to initialize: %3.3lf\n", time2-time1);
#endif
  time1 = omp_get_wtime();
#pragma omp parallel for
  for (long i=0; i<NV_in; i++) {
	if((C[i] < 0)||(C[i]>numUniqueClusters))
		continue; //Not a valid cluster id
	long adj1 = vtxPtrIn[i];
	long adj2 = vtxPtrIn[i+1];
	map<long, long>::iterator localIterator;
        assert(C[i] < numUniqueClusters);
  	//Now look for all the neighbors of this cluster
	for(long j=adj1; j<adj2; j++) {
		long tail = vtxIndIn[j].tail; 
		assert(C[tail] < numUniqueClusters);			
		//Add the edge from one endpoint	
		if(C[i] >= C[tail]) {
                        omp_set_lock(&nlocks[C[i]]);  // Locking the cluster
	
			localIterator = cluPtrIn[C[i]]->find(C[tail]); //Check if it exists			
			if( localIterator != cluPtrIn[C[i]]->end() ) {	//Already exists
                                localIterator->second += (long)vtxIndIn[j].weight;
			} else {
				(*(cluPtrIn[C[i]]))[C[tail]] = (long)vtxIndIn[j].weight; //Self-edge
				__sync_fetch_and_add(&vtxPtrOut[C[i]+1], 1);
				if(C[i] == C[tail])
                                	__sync_fetch_and_add(&NE_self, 1); //Keep track of self #edges 
				if(C[i] > C[tail]) {
					__sync_fetch_and_add(&NE_out, 1); //Keep track of non-self #edges
					__sync_fetch_and_add(&vtxPtrOut[C[tail]+1], 1); //Count edge j-->i
				}
			}
                        
                        omp_unset_lock(&nlocks[C[i]]); // Unlocking the cluster
		} //End of if
	}//End of for(j)
  }//End of for(i)  
  //Prefix sum:
  for(long i=0; i<NV_out; i++) {
	vtxPtrOut[i+1] += vtxPtrOut[i];
  }
  
  time2 = omp_get_wtime();
  TotTime += (time2-time1);
  printf("NE_out= %ld   NE_self= %ld\n", NE_out, NE_self);
  printf("These should match: %ld == %ld\n",(2*NE_out + NE_self), vtxPtrOut[NV_out]);
#ifdef PRINT_DETAILED_STATS_
  printf("Time to count edges: %3.3lf\n", time2-time1);
#endif
  assert(vtxPtrOut[NV_out] == (NE_out*2+NE_self)); //Sanity check
  
  time1 = omp_get_wtime();
  // Step 3 : build the edge list:
  long numEdges   = vtxPtrOut[NV_out];
  long realEdges  = NE_out + NE_self; //Self-loops appear once, others appear twice
  edge *vtxIndOut = (edge *) malloc (numEdges * sizeof(edge));
  assert (vtxIndOut != 0);
  long *Added = (long *) malloc (NV_out * sizeof(long)); //Keep track of what got added
  assert (Added != 0);

#pragma omp parallel for
  for (long i=0; i<NV_out; i++) {
	Added[i] = 0;
  }  
  //Now add the edges in no particular order
#pragma omp parallel for
  for (long i=0; i<NV_out; i++) {
	long Where;
	map<long, long>::iterator localIterator = cluPtrIn[i]->begin();
	//Now go through the other edges:
	while ( localIterator != cluPtrIn[i]->end()) {
		Where = vtxPtrOut[i] + __sync_fetch_and_add(&Added[i], 1);
		vtxIndOut[Where].head = i; //Head
		vtxIndOut[Where].tail = localIterator->first; //Tail
		vtxIndOut[Where].weight = localIterator->second; //Weight
		if(i != localIterator->first) {		
			Where = vtxPtrOut[localIterator->first] + __sync_fetch_and_add(&Added[localIterator->first], 1);
			vtxIndOut[Where].head = localIterator->first;
			vtxIndOut[Where].tail = i; //Tail
			vtxIndOut[Where].weight = localIterator->second; //Weight
			//printf("%d\n",localIterator->first);
		}
		localIterator++;
	}	
  }//End of for(i)  
  time2 = omp_get_wtime();
  TotTime += (time2-time1);
#ifdef PRINT_DETAILED_STATS_
  printf("Time to build the graph: %3.3lf\n", time2-time1);
  printf("Total time: %3.3lf\n", TotTime);
#endif
#ifdef PRINT_TERSE_STATS_
  printf("Total time to build next phase: %3.3lf\n", TotTime);
#endif
  // Set the pointers
  Gout->numVertices  = NV_out;
  Gout->sVertices    = NV_out;
  //Note: Self-loops are represented ONCE, but others appear TWICE
  Gout->numEdges     = realEdges; //Add self loops to the #edges
  Gout->edgeListPtrs = vtxPtrOut;
  Gout->edgeList     = vtxIndOut;
	
  //Clean up
  free(Added);
#pragma omp parallel for
  for (long i=0; i<numUniqueClusters; i++)
	delete cluPtrIn[i];	
  free(cluPtrIn);

#pragma omp parallel for
  for (long i=0; i<numUniqueClusters; i++) {
    omp_destroy_lock(&nlocks[i]); 
  }
  free(nlocks);
  
  return TotTime;
}//End of buildNextLevelGraph2()


