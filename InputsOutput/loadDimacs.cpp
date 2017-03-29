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

void parse_Dimacs9FormatDirectedNewD(graph * G, char *fileName) {
  printf("Parsing a DIMACS-9 formatted file as a general graph...\n");
  printf("WARNING: Assumes that the graph is directed -- an edge is stored only once.\n");
  printf("       : Graph will be stored as undirected, each edge appears twice.\n");
  int nthreads = 0;
#pragma omp parallel
  {
    nthreads = omp_get_num_threads();
  }
  printf("parse_Dimacs9FormatDirectedNewD: Number of threads: %d\n ", nthreads);
  
  double time1, time2;
  FILE *file = fopen(fileName, "r");
  if (file == NULL) {
    printf("Cannot open the input file: %s\n",fileName);
    exit(1);
  }  
  char line[1024], LS1[25], LS2[25];
  //Ignore the Comment lines starting with "c"
  do {
    fgets(line, 1024, file);
  } while ( line[0] == 'c');  
  //Expecting a problem line here:  p sp n m  
  long NV = 0, NE=0;  
  //Parse the problem line:
  if(line[0] == 'p') {
    if (sscanf(line, "%s %s %ld %ld", LS1, LS2, &NV, &NE) != 4) {
      printf("parse_Dimacs9(): bad file format - 01");
      exit(1);   
    }
  } else {
    printf("parse_Dimacs9(): bad file format, expecting a line starting with p\n");
    printf("%s\n", line);
    exit(1);
  }
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
    //Vertex lines:  degree vlabel xcoord ycoord
    fscanf(file, "%s %ld %ld %lf", LS1, &Si, &Ti, &Twt);
    assert((Si > 0)&&(Si <= NV));
    assert((Ti > 0)&&(Ti <= NV));
    tmpEdgeList[i].head   = Si-1;       //The S index 
    tmpEdgeList[i].tail   = Ti-1;    //The T index: One-based indexing
    tmpEdgeList[i].weight = fabs((double)Twt); //Make it positive and cast to Double
  }//End of outer for loop
  fclose(file); //Close the file
  time2 = omp_get_wtime(); 
  printf("Done reading from file: NE= %ld. Time= %lf\n", NE, time2-time1);
  
  //Remove duplicate entries:
 /* long NewEdges = removeEdges(NV, NE, tmpEdgeList);
  if (NewEdges < NE) {
    printf("Number of duplicate entries detected: %ld\n", NE-NewEdges);
    NE = NewEdges; //Only look at clean edges
  } else {
    printf("No dubplicates found.\n");
  }*/
  ///////////
  time1 = omp_get_wtime();
  long *edgeListPtr = (long *)  malloc((NV+1) * sizeof(long));
  assert(edgeListPtr != NULL);
  edge *edgeList = (edge *) malloc( 2*NE * sizeof(edge)); //Every edge stored twice
  assert( edgeList != NULL);
  time2 = omp_get_wtime();
  printf("Time for allocating memory for storing graph = %lf\n", time2 - time1);
#pragma omp parallel for
  for (long i=0; i <= NV; i++)
    edgeListPtr[i] = 0; //For first touch purposes
  
  //////Build the EdgeListPtr Array: Cumulative addition 
  time1 = omp_get_wtime();
#pragma omp parallel for
  for(long i=0; i<NE; i++) {
    __sync_fetch_and_add(&edgeListPtr[tmpEdgeList[i].head+1], 1); //Leave 0th position intact
    __sync_fetch_and_add(&edgeListPtr[tmpEdgeList[i].tail+1], 1);
  }
  for (long i=0; i<NV; i++) {
    edgeListPtr[i+1] += edgeListPtr[i]; //Prefix Sum:
  }
  //The last element of Cumulative will hold the total number of characters
  time2 = omp_get_wtime();
  printf("Done cumulative addition for edgeListPtrs:  %9.6lf sec.\n", time2 - time1);
  printf("Sanity Check: 2|E| = %ld, edgeListPtr[NV]= %ld\n", NE*2, edgeListPtr[NV]);
  printf("*********** (%ld)\n", NV);
  
  printf("About to build edgeList...\n");
  time1 = omp_get_wtime();
  //Keep track of how many edges have been added for a vertex:
  long  *added  = (long *)  malloc( NV  * sizeof(long)); assert( added != NULL);
#pragma omp parallel for
  for (long i = 0; i < NV; i++) 
    added[i] = 0;

  //Build the edgeList from edgeListTmp:
#pragma omp parallel for
  for(long i=0; i<NE; i++) {
    long head      = tmpEdgeList[i].head;
    long tail      = tmpEdgeList[i].tail;
    double weight  = tmpEdgeList[i].weight;
    
    long Where = edgeListPtr[head] + __sync_fetch_and_add(&added[head], 1);   
    edgeList[Where].head = head;
    edgeList[Where].tail = tail;
    edgeList[Where].weight = weight;
    //Now add the counter-edge:
    Where = edgeListPtr[tail] + __sync_fetch_and_add(&added[tail], 1);
    edgeList[Where].head = tail;
    edgeList[Where].tail = head;
    edgeList[Where].weight = weight;    
  }
  time2 = omp_get_wtime();
  printf("Time for building edgeList = %lf\n", time2 - time1);

  G->sVertices    = NV;
  G->numVertices  = NV;
  G->numEdges     = NE;
  G->edgeListPtrs = edgeListPtr;
  G->edgeList     = edgeList;
  
  //Clean up
  free(tmpEdgeList);
  free(added);

}//End of parse_Dimacs9FormatDirectedNewD()