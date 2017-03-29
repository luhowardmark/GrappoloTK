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
#include "input_output.h"
#include "defs.h"
#include "sstream"
#include "utilityStringTokenizer.hpp"

void parse_DirectedEdgeList(graph * G, char *fileName) {
  printf("Parsing a DoulbedEdgeList formatted file as a general graph...\n");
  printf("WARNING: Assumes that the graph is undirected -- an edge is stored twince.\n");
  int nthreads = 0;

#pragma omp parallel
  {
    nthreads = omp_get_num_threads();
  }
  printf("parse_EdgeListUndirectional: Number of threads: %d\n ", nthreads);
  
  double time1, time2;
  FILE *file = fopen(fileName, "r");
  if (file == NULL) {
    printf("Cannot open the input file: %s\n",fileName);
    exit(1);
  }  
	
	long NV=0, NE=0;
	long nv1, nv2;
	while(!feof(file))
	{
		fscanf(file, "%ld %ld", &nv1, &nv2);
		if(nv1 > NV)
			NV = nv1;
		if(nv2 > NV)
			NV = nv2;
		NE++;
	}
	NV++; NE--;
	fclose(file);
  
	file = fopen(fileName, "r");


  printf("|V|= %ld, |E|= %ld \n", NV, NE);  
  printf("Weights will be converted to positive numbers.\n");
  /*---------------------------------------------------------------------*/
  /* Read edge list: U V W                                             */
  /*---------------------------------------------------------------------*/  
  edge *tmpEdgeList = (edge *) malloc( NE * sizeof(edge)); //Every edge stored ONCE
  assert( tmpEdgeList != NULL);
  long Si, Ti;
  double Twt;
  time1 = omp_get_wtime();
  for (long i = 0; i < NE; i++) {
    fscanf(file, "%ld %ld", &Si, &Ti);
    assert((Si >= 0)&&(Si < NV));
    assert((Ti >= 0)&&(Ti < NV));
    tmpEdgeList[i].head   = Si;       //The S index 
    tmpEdgeList[i].tail   = Ti;    //The T index: Zero-based indexing
    tmpEdgeList[i].weight = 1; //Make it positive and cast to Double
		//printf("%d %d\n",Si,Ti);
  }//End of outer for loop
  
  fclose(file); //Close the file
  time2 = omp_get_wtime(); 
  printf("Done reading from file: NE= %ld. Time= %lf\n", NE, time2-time1);
  
  ///////////
  time1 = omp_get_wtime();
  long *edgeListPtr = (long *)  malloc((NV+1) * sizeof(long));
  assert(edgeListPtr != NULL);
  edge *edgeList = (edge *) malloc( NE * sizeof(edge)); //Every edge stored twice
  assert( edgeList != NULL);
  time2 = omp_get_wtime();
  printf("Time for allocating memory for storing graph = %lf\n", time2 - time1);

#pragma omp parallel for
  for (long i=0; i <= NV; i++)
    edgeListPtr[i] = 0; //For first touch purposes
  
  //////Build the EdgeListPtr Array: Cumulative addition 
  time1 = omp_get_wtime();
//#pragma omp parallel for
  for(long i=0; i<NE; i++) {
    __sync_fetch_and_add(&edgeListPtr[tmpEdgeList[i].head+1], 1); //Leave 0th position intact
  }
  for (long i=0; i<NV; i++) {
    edgeListPtr[i+1] += edgeListPtr[i]; //Prefix Sum:
  }
  //The last element of Cumulative will hold the total number of characters
  time2 = omp_get_wtime();
  printf("Done cumulative addition for edgeListPtrs:  %9.6lf sec.\n", time2 - time1);
  printf("Sanity Check: |E| = %ld, edgeListPtr[NV]= %ld\n", NE, edgeListPtr[NV]);
  printf("*********** (%ld)\n", NV);
  
  //time1 = omp_get_wtime();
  //Keep track of how many edges have been added for a vertex:
  printf("About to allocate for added vector: %ld\n", NV);
  long  *added  = (long *)  malloc( NV  * sizeof(long));
  printf("Done allocating memory fors added vector\n");
  assert( added != NULL);
	
	/*for(long i= 0; i<NV;i++)
	{
		printf("%d\n",edgeListPtr[i+1]);
	}*/
	
#pragma omp parallel for
  for (long i = 0; i < NV; i++) 
    added[i] = 0;
  
  printf("About to build edgeList...\n");
  //Build the edgeList from edgeListTmp:
//#pragma omp parallel for
  for(long i=0; i<NE; i++) {
    long head      = tmpEdgeList[i].head;
    long tail      = tmpEdgeList[i].tail;
    double weight  = tmpEdgeList[i].weight;
    
    long Where = edgeListPtr[head] + __sync_fetch_and_add(&added[head], 1);   
    edgeList[Where].head = head;
    edgeList[Where].tail = tail;
    edgeList[Where].weight = weight;
  }
  //time2 = omp_get_wtime();
  printf("Time for building edgeList = %lf\n", time2 - time1);

  G->sVertices    = NV;
  G->numVertices  = NV;
  G->numEdges     = NE;
  G->edgeListPtrs = edgeListPtr;
  G->edgeList     = edgeList;
  
  //Clean up*/
  free(tmpEdgeList);
  free(added);

}//End of parse_DirectedEdgeList()




void parse_UndirectedEdgeList(graph * G, char *fileName) {
  printf("Parsing a SingledEdgeList formatted file as a general graph...\n");
  printf("WARNING: Assumes that the graph is undirected -- an edge is stored twince.\n");
  int nthreads = 0;

	#pragma omp parallel
  {
    nthreads = omp_get_num_threads();
  }
  
  double time1, time2;
  FILE *file = fopen(fileName, "r");
  if (file == NULL) {
    printf("Cannot open the input file: %s\n",fileName);
    exit(1);
  }  
	
	long NV=0, NE=0;
	long nv1, nv2;
	while(!feof(file))
	{
		fscanf(file, "%ld %ld", &nv1, &nv2);
		if(nv1 > NV)
			NV = nv1;
		if(nv2 > NV)
			NV = nv2;
		NE++;
	}
	NV++; NE--;
	NE*=2;
	fclose(file);
  
	file = fopen(fileName, "r");
  printf("|V|= %ld, |E|= %ld \n", NV, NE);  
  printf("Weights will be converted to positive numbers.\n");
  /*---------------------------------------------------------------------*/
  /* Read edge list: a U V W                                             */
  /*---------------------------------------------------------------------*/  
  edge *tmpEdgeList = (edge *) malloc( NE * sizeof(edge)); //Every edge stored ONCE
  assert( tmpEdgeList != NULL);
  long Si, Ti;
  double Twt;
  time1 = omp_get_wtime();
  for (long i = 0; i < NE; i++) {
    fscanf(file, "%ld %ld", &Si, &Ti);
    assert((Si >= 0)&&(Si < NV));
    assert((Ti >= 0)&&(Ti < NV));
    tmpEdgeList[i].head   = Si;       //The S index 
    tmpEdgeList[i].tail   = Ti;    //The T index: Zero-based indexing
    tmpEdgeList[i].weight = 1; //Make it positive and cast to Double
    i++;
    tmpEdgeList[i].head = Ti;
    tmpEdgeList[i].tail = Si;
    tmpEdgeList[i].weight = 1;


  }//End of outer for loop
  printf("%d %d\n",Si,Ti);
  fclose(file); //Close the file
  time2 = omp_get_wtime(); 
  printf("Done reading from file: NE= %ld. Time= %lf\n", NE, time2-time1);
  
  ///////////
  time1 = omp_get_wtime();
  long *edgeListPtr = (long *)  malloc((NV+1) * sizeof(long));
  assert(edgeListPtr != NULL);
  edge *edgeList = (edge *) malloc( NE * sizeof(edge)); //Every edge stored twice
  assert( edgeList != NULL);
  time2 = omp_get_wtime();
  printf("Time for allocating memory for storing graph = %lf\n", time2 - time1);

#pragma omp parallel for
  for (long i=0; i <= NV; i++)
    edgeListPtr[i] = 0; //For first touch purposes
  
  //////Build the EdgeListPtr Array: Cumulative addition 
  time1 = omp_get_wtime();
//#pragma omp parallel for
  for(long i=0; i<NE; i++) {
    __sync_fetch_and_add(&edgeListPtr[tmpEdgeList[i].head+1], 1); //Leave 0th position intact
  }
  for (long i=0; i<NV; i++) {
    edgeListPtr[i+1] += edgeListPtr[i]; //Prefix Sum:
		//printf("%d ",edgeListPtr[i]);
  }
  //The last element of Cumulative will hold the total number of characters
  time2 = omp_get_wtime();
  printf("Done cumulative addition for edgeListPtrs:  %9.6lf sec.\n", time2 - time1);
  printf("Sanity Check: |E| = %ld, edgeListPtr[NV]= %ld\n", NE, edgeListPtr[NV]);
  printf("*********** (%ld)\n", NV);
  
  //time1 = omp_get_wtime();
  //Keep track of how many edges have been added for a vertex:
  printf("About to allocate for added vector: %ld\n", NV);
  long  *added  = (long *)  malloc( NV  * sizeof(long));
  printf("Done allocating memory fors added vector\n");
  assert( added != NULL);
#pragma omp parallel for
  for (long i = 0; i < NV; i++) 
    added[i] = 0;
  
  printf("About to build edgeList...\n");
  //Build the edgeList from edgeListTmp:
//#pragma omp parallel for
  for(long i=0; i<NE; i++) {
    long head      = tmpEdgeList[i].head;
    long tail      = tmpEdgeList[i].tail;
    double weight  = tmpEdgeList[i].weight;
    
    long Where = edgeListPtr[head] + __sync_fetch_and_add(&added[head], 1);   
    edgeList[Where].head = head;
    edgeList[Where].tail = tail;
    edgeList[Where].weight = weight;

  }
  //time2 = omp_get_wtime();
  printf("Time for building edgeList = %lf\n", time2 - time1);

  G->sVertices    = NV;
  G->numVertices  = NV;
  G->numEdges     = NE/2;
  G->edgeListPtrs = edgeListPtr;
  G->edgeList     = edgeList;
  
  free(tmpEdgeList);
  free(added);

}//End of parse_UndirectedEdgeList()
