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
#include "utilityStringTokenizer.hpp"

void parse_SNAP(graph * G, char *fileName) {
  printf("Parsing a SNAP formatted file as a general graph...\n");
  printf("WARNING: Assumes that the graph is directed -- an edge is stored only once.\n");
  printf("       : Graph will be stored as undirected, each edge appears twice.\n");
  int nthreads = 0;
#pragma omp parallel
  {
    nthreads = omp_get_num_threads();
  }
  printf("parse_Dimacs9FormatDirectedNewD: Number of threads: %d\n ", nthreads);
  
  long   NV=0,  NE=0;
  string oneLine, myDelimiter(" "), myDelimiter2("\t"), oneWord; //Delimiter is a blank space
  //ifstream fin;
  char comment;

  double time1, time2;
  /*
  FILE *file = fopen(fileName, "r");
  if (file == NULL) {
    printf("Cannot open the input file: %s\n",fileName);
    exit(1);
  }  
  */
  ifstream fin;
  fin.open(fileName);
  if(!fin) {
    cerr<<"Within Function: loadSNAPFileFormat() \n";
    cerr<<"Could not open the file.. \n";
    exit(1);
  }
  
  do { //Parse the comment lines for problem size
    getline(fin, oneLine);
    cout<<"Read line: "<<oneLine<<endl;
    comment = oneLine[0];
    if (comment == '#') { //Check if this line has problem sizes
      StringTokenizer* ST = new StringTokenizer(oneLine, myDelimiter);
      if ( ST->HasMoreTokens() )
	oneWord = ST->GetNextToken(); //Ignore #
      if ( ST->HasMoreTokens() )
	oneWord = ST->GetNextToken(); //Ignore #
      if(oneWord == "Nodes:") {	
	//# Nodes: 65608366 Edges: 1806067135	
	NV  = atol( ST->GetNextToken().c_str() ); //Number of Vertices
	oneWord = ST->GetNextToken(); //Ignore Edges:
	NE  = atol( ST->GetNextToken().c_str() ); //Number of Edges   
      }
      delete ST;
    }
  } while ( comment == '#');
      
  printf("|V|= %ld, |E|= %ld \n", NV, NE);  
  printf("Weight of 1 will be assigned to each edge.\n");
  cout << oneLine <<endl;
  /*---------------------------------------------------------------------*/
  /* Read edge list: a U V W                                             */
  /*---------------------------------------------------------------------*/  
  edge *tmpEdgeList = (edge *) malloc( NE * sizeof(edge)); //Every edge stored ONCE
  assert( tmpEdgeList != NULL);
  long Si, Ti;

  map<long, long> clusterLocalMap; //Renumber vertices contiguously from zero
  map<long, long>::iterator storedAlready;
  long numUniqueVertices = 0;
    
  //Parse the first edge already read from the file and stored in oneLine
  long i=0;
  do {
    StringTokenizer* ST = new StringTokenizer(oneLine, myDelimiter2);
    if ( ST->HasMoreTokens() )
      Si  = atol( ST->GetNextToken().c_str() );   
    if ( ST->HasMoreTokens() )
      Ti  = atol( ST->GetNextToken().c_str() ); 
    delete ST;  
    
    storedAlready = clusterLocalMap.find(Si); //Check if it already exists
    if( storedAlready != clusterLocalMap.end() ) {	//Already exists
      Si = storedAlready->second; //Renumber the cluster id
    } else {
      clusterLocalMap[Si] = numUniqueVertices; //Does not exist, add to the map
      Si = numUniqueVertices; //Renumber the vertex id
      numUniqueVertices++; //Increment the number
    }
    
    storedAlready = clusterLocalMap.find(Ti); //Check if it already exists
    if( storedAlready != clusterLocalMap.end() ) {	//Already exists
      Ti = storedAlready->second; //Renumber the cluster id
    } else {
      clusterLocalMap[Ti] = numUniqueVertices; //Does not exist, add to the map
      Ti = numUniqueVertices; //Renumber the vertex id
      numUniqueVertices++; //Increment the number
    }
    tmpEdgeList[i].head   = Si;  //The S index 
    tmpEdgeList[i].tail   = Ti;  //The T index: One-based indexing
    tmpEdgeList[i].weight = 1;     //default weight of one 
    //cout<<" Adding edge ("<<Si<<", "<<Ti<<")\n";
    i++;
    //Read-in the next line
    getline(fin, oneLine);
    if ((i % 999999) == 1) {
      cout <<"Reading Line: "<<i<<endl;
    }
  } while ( !fin.eof() );//End of while

  fin.close(); //Close the file
  time2 = omp_get_wtime(); 
  printf("Done reading from file: NE= %ld. Time= %lf\n", NE, time2-time1);
  
  
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

}