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

long removeEdges(long NV, long NE, edge *edgeList) {
  printf("Within removeEdges()\n");
  long NGE = 0;
  long *head = (long *) malloc(NV * sizeof(long));     /* head of linked list points to an edge */   
  long *next = (long *) malloc(NE * sizeof(long));     /* ptr to next edge in linked list       */
  
/* Initialize linked lists */
  for (long i = 0; i < NV; i++) head[i] = -1;
  for (long i = 0; i < NE; i++) next[i] = -2;
  
  for (long i = 0; i < NE; i++) {
    long sv  = edgeList[i].head;
    long ev  = edgeList[i].tail;
    if (sv == ev) continue;    /* self edge */    
    long * ptr = head + sv;     /* start at head of list for this key */    
    while (1) {
      long edgeId = *ptr;      
      if (edgeId == -1) {         /* at the end of the list */
	edgeId = *ptr;             /* lock ptr               */	
	if (edgeId == -1) {       /* if still end of list   */	  
	  long newId = NGE;
	  NGE++;     /* increment number of good edges */
	  //edgeList[i].id = newId;                 /* set id of edge                 */	  	  
	  next[i] = -1;                           /* insert edge in linked list     */
	  *ptr = i;
	  break;
	}	
	*ptr = edgeId;	
      } else 
	if (edgeList[edgeId].tail == ev) break;     /* duplicate edge */
	else  ptr = next + edgeId;
    }
  }  
  /* Move good edges to front of edgeList                    */
  /* While edge i is a bad edge, swap with last edge in list */
  for (long i = 0; i < NGE; i++) {
    while (next[i] == -2) {
      long k = NE - 1;
      NE--;
      edgeList[i] = edgeList[k];
      next[i] = next[k];
    }
  }
  printf("About to free memory\n");
  free(head);
  free(next);
  printf("Exiting removeEdges()\n");
  return NGE;
}//End of removeEdges()

/* Since graph is undirected, sort each edge head --> tail AND tail --> head */
void SortEdgesUndirected(long NV, long NE, edge *list1, edge *list2, long *ptrs) {
  for (long i = 0; i < NV + 2; i++) 
    ptrs[i] = 0;
  ptrs += 2;

  /* Histogram key values */
  for (long i = 0; i < NE; i++) {
    ptrs[list1[i].head]++;
    ptrs[list1[i].tail]++;
  }
  /* Compute start index of each bucket */
  for (long i = 1; i < NV; i++) 
    ptrs[i] += ptrs[i-1];
  ptrs--;

  /* Move edges into its bucket's segment */
  for (long i = 0; i < NE; i++) {
    long head   = list1[i].head;
    long index          = ptrs[head]++;
    //list2[index].id     = list1[i].id;
    list2[index].head   = list1[i].head;
    list2[index].tail   = list1[i].tail;
    list2[index].weight = list1[i].weight;

    long tail   = list1[i].tail;
    index               = ptrs[tail]++;
    //list2[index].id     = list1[i].id;
    list2[index].head   = list1[i].tail;
    list2[index].tail   = list1[i].head;
    list2[index].weight = list1[i].weight;
  } 
}//End of SortEdgesUndirected2()

/* Sort each node's neighbors by tail from smallest to largest. */
void SortNodeEdgesByIndex(long NV, edge *list1, edge *list2, long *ptrs) {  
  for (long i = 0; i < NV; i++) {
    edge *edges1 = list1 + ptrs[i];
    edge *edges2 = list2 + ptrs[i];
    long size    = ptrs[i+1] - ptrs[i];

    /* Merge Sort */
    for (long skip = 2; skip < 2 * size; skip *= 2) {
      for (long sect = 0; sect < size; sect += skip)  {
	long j = sect;
	long l = sect;
	long half_skip = skip / 2;
	long k = sect + half_skip;
	
	long j_limit = (j + half_skip < size) ? j + half_skip : size;
	long k_limit = (k + half_skip < size) ? k + half_skip : size;
	
	while ((j < j_limit) && (k < k_limit)) {
	  if   (edges1[j].tail < edges1[k].tail) {edges2[l] = edges1[j]; j++; l++;}
	  else                                   {edges2[l] = edges1[k]; k++; l++;}
	}	
	while (j < j_limit) {edges2[l] = edges1[j]; j++; l++;}
	while (k < k_limit) {edges2[l] = edges1[k]; k++; l++;}
      }
      edge *tmp = edges1;
      edges1 = edges2;
      edges2 = tmp;
    }
    // result is in list2, so move to list1
    if (edges1 == list2 + ptrs[i])
      for (long j = ptrs[i]; j < ptrs[i+1]; j++) list1[j] = list2[j];
  } 
}//End of SortNodeEdgesByIndex2()

void parse_PajekFormat(graph * G, char *fileName) {
  printf("Parsing a Pajek File...\n");
  int nthreads;
#pragma omp parallel
  {
    nthreads = omp_get_num_threads();    
  }
  printf("parse_Pajek: Number of threads: %d\n", nthreads);
  
  double time1, time2;
  FILE *file = fopen(fileName, "r");
  if (file == NULL) {
    printf("Cannot open the input file: %s\n",fileName);
    exit(1);
  }
  //Parse the first line:
  char line[1024];
  fgets(line, 1024, file);  
  char  LS1[25], LS2[25];
  long NV = 0, NE=0;
  if (sscanf(line, "%s %s", LS1, LS2) != 2) {
    printf("parse_Pajek(): bad file format - 01");
    exit(1);
  }  
  //printf("(%s) --- (%s) \n", LS1, LS2);
  if ( strcmp(LS1,"*Vertices")!= 0 ) {
    printf("Error: The first line should start with *Vertices word \n");
    exit(1);
  }
  NV = atol(LS2);
  printf("|V|= %ld \n", NV);
  /* Ignore all the vertex lines */
  for (long i=0; i <= NV; i++) {
    fgets(line, 1024, file);
  }
  if (sscanf(line, "%s", LS1) != 1) {
    printf("parse_Pajek(): bad file format - 02");
    exit(1);
  }
  if ( strcmp(LS1,"*Edges")!= 0 ) {
    printf("Error: The next line should start with *Edges word \n");
    exit(1);
  }
  printf("Line read: %s\n",line);
  /*---------------------------------------------------------------------*/
  /* Read edge list                                                      */
  /* (i , j, value ) 1-based index                                       */
  /*---------------------------------------------------------------------*/  
  edge *edgeListTmp; //Read the data in a temporary list
  long Si, Ti;
  double weight = 1;
  long edgeEstimate = NV * 2; //25% density -- not valid for large graphs
  edgeListTmp = (edge *) malloc( edgeEstimate * sizeof(edge));
  assert(edgeListTmp != 0);

  printf("Parsing edges -- no weights\n");
  
  //while (fscanf(file, "%ld %ld %lf", &Si, &Ti, &weight) != EOF) {
  while (fscanf(file, "%ld %ld", &Si, &Ti) != EOF) {
    if(NE >= edgeEstimate) {
      printf("Temporary buffer is not enough. \n");
      exit(1);
    }
    //printf("%ld -- %ld\n", Si, Ti);
    Si--; Ti--;            // One-based indexing
    assert((Si >= 0)&&(Si < NV));
    assert((Ti >= 0)&&(Ti < NV));
    if(Si%99999 == 0)
      printf("%ld -- %ld\n", Si, Ti);
    if (Si == Ti) //Ignore self-loops
      continue; 
    //weight = fabs(weight); //Make it positive    : Leave it as is
    weight = 1.0; //Make it positive    : Leave it as is
    edgeListTmp[NE].head = Si;       //The S index
    edgeListTmp[NE].tail = Ti;    //The T index
    edgeListTmp[NE].weight = weight; //The value
    NE++;
  }
  fclose(file); //Close the file
  printf("Done reading from file.\n");
  printf("|V|= %ld, |E|= %ld \n", NV, NE);
  
  //Remove duplicate entries:
  long NewEdges = removeEdges(NV, NE, edgeListTmp);
  if (NewEdges < NE) {
    printf("Number of duplicate entries detected: %ld\n", NE-NewEdges);
    NE = NewEdges; //Only look at clean edges
  }
  
  //Allocate for Edge Pointer and keep track of degree for each vertex
  long  *edgeListPtr = (long *)  malloc((NV+1) * sizeof(long));
#pragma omp parallel for
  for (long i=0; i <= NV; i++)
    edgeListPtr[i] = 0; //For first touch purposes

#pragma omp parallel for
  for(long i=0; i<NE; i++) {
    __sync_fetch_and_add(&edgeListPtr[edgeListTmp[i].head + 1], 1); //Plus one to take care of the zeroth location
    __sync_fetch_and_add(&edgeListPtr[edgeListTmp[i].tail + 1], 1);
  }
  
  //////Build the EdgeListPtr Array: Cumulative addition 
  time1 = omp_get_wtime();
  for (long i=0; i<NV; i++) {
    edgeListPtr[i+1] += edgeListPtr[i]; //Prefix Sum:
  }
  //The last element of Cumulative will hold the total number of characters
  time2 = omp_get_wtime();
  printf("Done cumulative addition for edgeListPtrs:  %9.6lf sec.\n", time2 - time1);
  printf("Sanity Check: 2|E| = %ld, edgeListPtr[NV]= %ld\n", NE*2, edgeListPtr[NV]);
  
  /*---------------------------------------------------------------------*/
  /* Allocate memory for G & Build it                                    */
  /*---------------------------------------------------------------------*/    
  printf("About to allocate memory for graph data structures\n");
  time1 = omp_get_wtime();
  edge *edgeList = (edge *) malloc ((2*NE) * sizeof(edge)); //Every edge stored twice
  assert(edgeList != 0);
  //Keep track of how many edges have been added for a vertex:
  long  *added = (long *)  malloc (NV * sizeof(long)); assert (added != 0);
#pragma omp parallel for
  for (long i = 0; i < NV; i++) 
    added[i] = 0;
  time2 = omp_get_wtime();  
  printf("Time for allocating memory for edgeList = %lf\n", time2 - time1);
  
  time1 = omp_get_wtime();
  
  printf("About to build edgeList...\n");
  //Build the edgeList from edgeListTmp:
#pragma omp parallel for
  for(long i=0; i<NE; i++) {
    long head  = edgeListTmp[i].head;
    long tail  = edgeListTmp[i].tail;
    double weight      = edgeListTmp[i].weight;
    
    long Where = edgeListPtr[head] + __sync_fetch_and_add(&added[head], 1);   
    edgeList[Where].head = head; 
    edgeList[Where].tail = tail;
    edgeList[Where].weight = weight;
    //added[head]++;
    //Now add the counter-edge:
    Where = edgeListPtr[tail] + __sync_fetch_and_add(&added[tail], 1);
    edgeList[Where].head = tail;
    edgeList[Where].tail = head;
    edgeList[Where].weight = weight;
    //added[tail]++;
  }
  time2 = omp_get_wtime();
  printf("Time for building edgeList = %lf\n", time2 - time1);
  
  G->sVertices    = NV;
  G->numVertices  = NV;
  G->numEdges     = NE;
  G->edgeListPtrs = edgeListPtr;
  G->edgeList     = edgeList;
  
  free(edgeListTmp);
  free(added);
}

/*-------------------------------------------------------*
 * This function reads a Pajek file and builds the graph
 *-------------------------------------------------------*/
//Assume every edge is stored twice:
void parse_PajekFormatUndirected(graph * G, char *fileName) {
  printf("Parsing a Pajek File *** Undirected ***...\n");
  int nthreads;
#pragma omp parallel
  {
    nthreads = omp_get_num_threads();    
  }
  printf("parse_Pajek_undirected: Number of threads: %d\n", nthreads);
  
  double time1, time2;
  FILE *file = fopen(fileName, "r");
  if (file == NULL) {
    printf("Cannot open the input file: %s\n",fileName);
    exit(1);
  }
  //Parse the first line:
  char line[1024];
  fgets(line, 1024, file);  
  char  LS1[25], LS2[25];
  long NV = 0, NE=0;
  if (sscanf(line, "%s %s", LS1, LS2) != 2) {
    printf("parse_Pajek(): bad file format - 01");
    exit(1);
  }  
  //printf("(%s) --- (%s) \n", LS1, LS2);
  if ( strcmp(LS1,"*Vertices")!= 0 ) {
    printf("Error: The first line should start with *Vertices word \n");
    exit(1);
  }
  NV = atol(LS2);
  printf("|V|= %ld \n", NV);
  /* Ignore all the vertex lines */
  for (long i=0; i <= NV; i++) {
    fgets(line, 1024, file);
  }
  printf("Done parsing through vertex lines\n");
  if (sscanf(line, "%s", LS1) != 1) {
    printf("parse_Pajek(): bad file format - 02");
    exit(1);
  }
  if ( strcmp(LS1,"*Edges")!= 0 ) {
    printf("Error: The next line should start with *Edges word \n");
    exit(1);
  }
  printf("About to read edges -- no weights\n");
  /*---------------------------------------------------------------------*/
  /* Read edge list                                                      */
  /* (i , j, value ) 1-based index                                       */
  /*---------------------------------------------------------------------*/
  
  edge *edgeListTmp; //Read the data in a temporary list
  long Si, Ti;
  double weight = 1;
  long edgeEstimate = NV * NV / 8; //12.5% density -- not valid for dense graphs
  edgeListTmp = (edge *) malloc( edgeEstimate * sizeof(edge));

  //while (fscanf(file, "%ld %ld %lf", &Si, &Ti, &weight) != EOF) {
  while (fscanf(file, "%ld %ld ", &Si, &Ti) != EOF) {
    Si--; Ti--;            // One-based indexing
    assert((Si >= 0)&&(Si < NV));
    assert((Ti >= 0)&&(Ti < NV));
    if (Si == Ti) //Ignore self-loops 
      continue; 
    //weight = fabs(weight); //Make it positive    : Leave it as is
    weight = 1.0; //Make it positive    : Leave it as is
    edgeListTmp[NE].head = Si;       //The S index
    edgeListTmp[NE].tail = Ti;    //The T index
    edgeListTmp[NE].weight = weight; //The value
    NE++;
  }
  fclose(file); //Close the file
  printf("Done reading from file.\n");
  printf("|V|= %ld, |E|= %ld \n", NV, NE);
  
  //Remove duplicate entries:
  /*
  long NewEdges = removeEdges(NV, NE, edgeListTmp);
  if (NewEdges < NE) {
    printf("Number of duplicate entries detected: %ld\n", NE-NewEdges);
    NE = NewEdges; //Only look at clean edges
  } else
    printf("No duplicates were found\n");
  */

  //Allocate for Edge Pointer and keep track of degree for each vertex
  long  *edgeListPtr = (long *) malloc((NV+1) * sizeof(long));
  assert(edgeListPtr != 0);

#pragma omp parallel for
  for (long i=0; i <= NV; i++) {
    edgeListPtr[i] = 0; //For first touch purposes
  }
#pragma omp parallel for
  for(long i=0; i<NE; i++) {
    __sync_fetch_and_add(&edgeListPtr[edgeListTmp[i].head + 1], 1); //Plus one to take care of the zeroth location
    //__sync_fetch_and_add(&edgeListPtr[edgeListTmp[i].tail + 1], 1); //No need
  }
  
  //////Build the EdgeListPtr Array: Cumulative addition 
  time1 = omp_get_wtime();
  for (long i=0; i<NV; i++) {
    edgeListPtr[i+1] += edgeListPtr[i]; //Prefix Sum:
  }
  //The last element of Cumulative will hold the total number of characters
  time2 = omp_get_wtime();
  printf("Done cumulative addition for edgeListPtrs:  %9.6lf sec.\n", time2 - time1);
  printf("Sanity Check: |E| = %ld, edgeListPtr[NV]= %ld\n", NE, edgeListPtr[NV]);
  /*---------------------------------------------------------------------*/
  /* Allocate memory for G & Build it                                    */
  /*---------------------------------------------------------------------*/    
  time1 = omp_get_wtime();
  printf("Size of edge: %ld  and size of NE*edge= %ld\n",  sizeof(edge), NE*sizeof(edge) );
  edge *edgeList = (edge *) malloc( NE * sizeof(edge)); //Every edge stored twice
  assert(edgeList != 0);
  //Keep track of how many edges have been added for a vertex:
  long  *added    = (long *)  malloc( NV  * sizeof(long)); assert (added != 0);
 #pragma omp parallel for
  for (long i = 0; i < NV; i++) 
    added[i] = 0;
  
  time2 = omp_get_wtime();
  printf("Time for allocating memory for edgeList = %lf\n", time2 - time1);
  
  time1 = omp_get_wtime();

  printf("About to build edgeList...\n");
  //Build the edgeList from edgeListTmp:
#pragma omp parallel for
  for(long i=0; i<NE; i++) {
    long head  = edgeListTmp[i].head;
    long tail  = edgeListTmp[i].tail;
    double weight      = edgeListTmp[i].weight;
    
    long Where = edgeListPtr[head] + __sync_fetch_and_add(&added[head], 1);   
    edgeList[Where].head = head; 
    edgeList[Where].tail = tail;
    edgeList[Where].weight = weight;
  }
  time2 = omp_get_wtime();
  printf("Time for building edgeList = %lf\n", time2 - time1);
  
  G->sVertices    = NV;
  G->numVertices  = NV;
  G->numEdges     = NE/2; //Each edge had been presented twice
  G->edgeListPtrs = edgeListPtr;
  G->edgeList     = edgeList;
  
  free(edgeListTmp);
  free(added);
}



