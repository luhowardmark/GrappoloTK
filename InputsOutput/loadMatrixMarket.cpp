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


/*-------------------------------------------------------*
 * This function reads a MATRIX MARKET file and builds the graph
 *-------------------------------------------------------*/
void parse_MatrixMarket(graph * G, char *fileName) {
  printf("Parsing a Matrix Market File...\n");
  int nthreads;
#pragma omp parallel
  {
    nthreads = omp_get_num_threads();    
  }
  printf("parse_MatrixMarket: Number of threads: %d\n", nthreads);

  double time1, time2;
  FILE *file = fopen(fileName, "r");
  if (file == NULL) {
    printf("Cannot open the input file: %s\n",fileName);
    exit(1);
  }
  
  /* -----      Read File in Matrix Market Format     ------ */
  //Parse the first line:
  char line[1024];
  fgets(line, 1024, file);
  char  LS1[25], LS2[25], LS3[25], LS4[25], LS5[25];    
  if (sscanf(line, "%s %s %s %s %s", LS1, LS2, LS3, LS4, LS5) != 5) {
    printf("parse_MatrixMarket(): bad file format - 01");
    exit(1);
  }
  printf("%s %s %s %s %s\n", LS1, LS2, LS3, LS4, LS5);
  if ( strcmp(LS1,"%%MatrixMarket") != 0 ) {
    printf("Error: The first line should start with %%MatrixMarket word \n");
    exit(1);
  }
  if ( !( strcmp(LS2,"matrix")==0 || strcmp(LS2,"Matrix")==0 || strcmp(LS2,"MATRIX")==0 ) ) {
    printf("Error: The Object should be matrix or Matrix or MATRIX \n");
    exit(1);
  }
  if ( !( strcmp(LS3,"coordinate")==0 || strcmp(LS3,"Coordinate")==0 || strcmp(LS3,"COORDINATE")==0) ) {
    printf("Error: The Object should be coordinate or Coordinate or COORDINATE \n");
    exit(1);
  }
  int isComplex = 0;    
  if ( strcmp(LS4,"complex")==0 || strcmp(LS4,"Complex")==0 || strcmp(LS4,"COMPLEX")==0 ) {
    isComplex = 1;
    printf("Warning: Will only read the real part. \n");
  }
  int isPattern = 0;
  if ( strcmp(LS4,"pattern")==0 || strcmp(LS4,"Pattern")==0 || strcmp(LS4,"PATTERN")==0 ) {
    isPattern = 1;
    printf("Note: Matrix type is Pattern. Will set all weights to 1.\n");
    //exit(1);
  }
  int isSymmetric = 0, isGeneral = 0;
  if ( strcmp(LS5,"general")==0 || strcmp(LS5,"General")==0 || strcmp(LS5,"GENERAL")==0 )
    isGeneral = 1;
  else {
    if ( strcmp(LS5,"symmetric")==0 || strcmp(LS5,"Symmetric")==0 || strcmp(LS5,"SYMMETRIC")==0 ) {
      isSymmetric = 1;
      printf("Note: Matrix type is Symmetric: Converting it into General type. \n");
    }
  }	
  if ( (isGeneral==0) && (isSymmetric==0) ) 	  {
    printf("Warning: Matrix type should be General or Symmetric. \n");
    exit(1);
  }
  
  /* Parse all comments starting with '%' symbol */
  do {
    fgets(line, 1024, file);
  } while ( line[0] == '%' );
  
  /* Read the matrix parameters */
  long NS=0, NT=0, NV = 0;
  long NE=0;
  if (sscanf(line, "%ld %ld %ld",&NS, &NT, &NE ) != 3) {
    printf("parse_MatrixMarket(): bad file format - 02");
    exit(1);
  }
  NV = NS + NT;
  printf("|S|= %ld, |T|= %ld, |E|= %ld \n", NS, NT, NE);

  /*---------------------------------------------------------------------*/
  /* Read edge list                                                      */
  /* S vertices: 0 to NS-1                                               */
  /* T vertices: NS to NS+NT-1                                           */
  /*---------------------------------------------------------------------*/
  //Allocate for Edge Pointer and keep track of degree for each vertex
  long  *edgeListPtr = (long *)  malloc((NV+1) * sizeof(long));
#pragma omp parallel for
  for (long i=0; i <= NV; i++)
    edgeListPtr[i] = 0; //For first touch purposes
  
  edge *edgeListTmp; //Read the data in a temporary list
  long newNNZ = 0;    //New edges because of symmetric matrices
  long Si, Ti;
  double weight = 1;
  if( isSymmetric == 1 ) {
    printf("Matrix is of type: Symmetric Real or Complex\n");
    printf("Weights will be converted to positive numbers.\n");
    edgeListTmp = (edge *) malloc(2 * NE * sizeof(edge));
    for (long i = 0; i < NE; i++) {
      if (isPattern == 1)
	fscanf(file, "%ld %ld", &Si, &Ti);
      else
	fscanf(file, "%ld %ld %lf", &Si, &Ti, &weight);
      Si--; Ti--;            // One-based indexing
      assert((Si >= 0)&&(Si < NV));
      assert((Ti >= 0)&&(Ti < NV));
      weight = fabs(weight); //Make it positive  : Leave it as is
      if ( Si == Ti ) {
	edgeListTmp[i].head = Si;       //The S index
	edgeListTmp[i].tail = NS+Ti;    //The T index 
	edgeListTmp[i].weight = weight; //The value
	edgeListPtr[Si+1]++;
	edgeListPtr[NS+Ti+1]++;
      }
      else { //an off diagonal element: Also store the upper part
	//LOWER PART:
	edgeListTmp[i].head = Si;       //The S index 
	edgeListTmp[i].tail = NS+Ti;    //The T index 
	edgeListTmp[i].weight = weight; //The value
	edgeListPtr[Si+1]++;
	edgeListPtr[NS+Ti+1]++;
	//UPPER PART:
	edgeListTmp[NE+newNNZ].head = Ti;       //The S index
	edgeListTmp[NE+newNNZ].tail = NS+Si;    //The T index
	edgeListTmp[NE+newNNZ].weight = weight; //The value
	newNNZ++; //Increment the number of edges
	edgeListPtr[Ti+1]++;
	edgeListPtr[NS+Si+1]++;
      }
    }
  } //End of Symmetric
  /////// General Real or Complex ///////
  else {
    printf("Matrix is of type: Unsymmetric Real or Complex\n");
    printf("Weights will be converted to positive numbers.\n");
    edgeListTmp = (edge *) malloc( NE * sizeof(edge));
    for (long i = 0; i < NE; i++) {
      if (isPattern == 1)
	fscanf(file, "%ld %ld", &Si, &Ti);
      else
	fscanf(file, "%ld %ld %lf", &Si, &Ti, &weight);
      //printf("(%d, %d) %lf\n",Si, Ti, weight);
      Si--; Ti--;            // One-based indexing
      assert((Si >= 0)&&(Si < NV));
      assert((Ti >= 0)&&(Ti < NV));
      weight = fabs(weight); //Make it positive    : Leave it as is
      edgeListTmp[i].head = Si;       //The S index
      edgeListTmp[i].tail = NS+Ti;    //The T index
      edgeListTmp[i].weight = weight; //The value
      edgeListPtr[Si+1]++;
      edgeListPtr[NS+Ti+1]++;
    }
  } //End of Real or Complex
  fclose(file); //Close the file
  printf("Done reading from file.\n");

  if( isSymmetric ) {
    printf("Modified the number of edges from %ld ",NE);
    NE += newNNZ; //#NNZ might change
    printf("to %ld \n",NE);
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
  time1 = omp_get_wtime();
  edge *edgeList = (edge *) malloc( 2*NE * sizeof(edge)); //Every edge stored twice
  assert(edgeList != 0);
  //Keep track of how many edges have been added for a vertex:
  long  *added    = (long *)  malloc( NV  * sizeof(long)); assert(added != 0);
#pragma omp parallel for
  for (long i = 0; i < NV; i++) 
    added[i] = 0;
  time2 = omp_get_wtime();
  printf("Time for allocating memory for marks and edgeList = %lf\n", time2 - time1);
  
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
    //Now add the counter-edge:
    Where = edgeListPtr[tail] + __sync_fetch_and_add(&added[tail], 1);
    edgeList[Where].head = tail;
    edgeList[Where].tail = head;
    edgeList[Where].weight = weight;
  }
  time2 = omp_get_wtime();
  printf("Time for building edgeList = %lf\n", time2 - time1);
  
  G->sVertices    = NS;
  G->numVertices  = NV;
  G->numEdges     = NE;
  G->edgeListPtrs = edgeListPtr;
  G->edgeList     = edgeList;
  
  free(edgeListTmp);
  free(added);
}


/*-------------------------------------------------------*
 * This function reads a MATRIX MARKET file and build the graph
 * graph is nonbipartite: each diagonal entry is a vertex, and 
 * each non-diagonal entry becomes an edge. Assume structural and 
 * numerical symmetry.
 *-------------------------------------------------------*/
void parse_MatrixMarket_Sym_AsGraph(graph * G, char *fileName) {
  printf("Parsing a Matrix Market File as a general graph...\n");
  int nthreads = 0;
#pragma omp parallel
  {
    nthreads = omp_get_num_threads();
  }
  printf("parse_MatrixMarket: Number of threads: %d\n ", nthreads);

  double time1, time2;
  FILE *file = fopen(fileName, "r");
  if (file == NULL) {
    printf("Cannot open the input file: %s\n",fileName);
    exit(1);
  }
  /* -----      Read File in Matrix Market Format     ------ */
  //Parse the first line:
  char line[1024];
  fgets(line, 1024, file);
  char  LS1[25], LS2[25], LS3[25], LS4[25], LS5[25];    
  if (sscanf(line, "%s %s %s %s %s", LS1, LS2, LS3, LS4, LS5) != 5) {
    printf("parse_MatrixMarket(): bad file format - 01");
    exit(1);
  }
  printf("%s %s %s %s %s\n", LS1, LS2, LS3, LS4, LS5);
  if ( strcmp(LS1,"%%MatrixMarket") != 0 ) {
    printf("Error: The first line should start with %%MatrixMarket word \n");
    exit(1);
  }
  if ( !( strcmp(LS2,"matrix")==0 || strcmp(LS2,"Matrix")==0 || strcmp(LS2,"MATRIX")==0 ) ) {
    printf("Error: The Object should be matrix or Matrix or MATRIX \n");
    exit(1);
  }
  if ( !( strcmp(LS3,"coordinate")==0 || strcmp(LS3,"Coordinate")==0 || strcmp(LS3,"COORDINATE")==0) ) {
    printf("Error: The Object should be coordinate or Coordinate or COORDINATE \n");
    exit(1);
  }
  int isComplex = 0;    
  if ( strcmp(LS4,"complex")==0 || strcmp(LS4,"Complex")==0 || strcmp(LS4,"COMPLEX")==0 ) {
    isComplex = 1;
    printf("Warning: Will only read the real part. \n");
  }
  int isPattern = 0;
  if ( strcmp(LS4,"pattern")==0 || strcmp(LS4,"Pattern")==0 || strcmp(LS4,"PATTERN")==0 ) {
    isPattern = 1;
    printf("Note: Matrix type is Pattern. Will set all weights to 1.\n");
    //exit(1);
  }
  int isSymmetric = 0, isGeneral = 0;
  if ( strcmp(LS5,"general")==0 || strcmp(LS5,"General")==0 || strcmp(LS5,"GENERAL")==0 )
    isGeneral = 1;
  else {
    if ( strcmp(LS5,"symmetric")==0 || strcmp(LS5,"Symmetric")==0 || strcmp(LS5,"SYMMETRIC")==0 ) {
      isSymmetric = 1;
      printf("Note: Matrix type is Symmetric: Converting it into General type. \n");
    }
  }	
  if ( isSymmetric==0 ) 	  {
    printf("Warning: Matrix type should be Symmetric for this routine. \n");
    exit(1);
  }
  
  /* Parse all comments starting with '%' symbol */
  do {
    fgets(line, 1024, file);
  } while ( line[0] == '%' );
  
  /* Read the matrix parameters */
  long NS=0, NT=0, NV = 0;
  long NE=0;
  if (sscanf(line, "%ld %ld %ld",&NS, &NT, &NE ) != 3) {
    printf("parse_MatrixMarket(): bad file format - 02");
    exit(1);
  }
  NV = NS;
  printf("|S|= %ld, |T|= %ld, |E|= %ld \n", NS, NT, NE);

  /*---------------------------------------------------------------------*/
  /* Read edge list                                                      */
  /* S vertices: 0 to NS-1                                               */
  /* T vertices: NS to NS+NT-1                                           */
  /*---------------------------------------------------------------------*/
  //Allocate for Edge Pointer and keep track of degree for each vertex
  long *edgeListPtr = (long *)  malloc((NV+1) * sizeof(long));
#pragma omp parallel for
  for (long i=0; i <= NV; i++)
    edgeListPtr[i] = 0; //For first touch purposes
  
  edge *edgeListTmp; //Read the data in a temporary list
  long newNNZ = 0;    //New edges because of symmetric matrices
  long Si, Ti;
  double weight = 1;
  printf("Matrix is of type: Symmetric Real or Complex\n");
  printf("Weights will be converted to positive numbers.\n");
  edgeListTmp = (edge *) malloc(2 * NE * sizeof(edge));
  for (long i = 0; i < NE; i++) {
    if (isPattern == 1)
      fscanf(file, "%ld %ld", &Si, &Ti);
    else
      fscanf(file, "%ld %ld %lf", &Si, &Ti, &weight);
    Si--; Ti--;            // One-based indexing
    assert((Si >= 0)&&(Si < NV));
    assert((Ti >= 0)&&(Ti < NV));
    weight = fabs(weight); //Make it positive  : Leave it as is
    if ( Si == Ti ) {
      //Do nothing...
    }
    else { //an off diagonal element: store the edge
      //LOWER PART:
      edgeListTmp[newNNZ].head = Si;       //The S index 
      edgeListTmp[newNNZ].tail = Ti;       //The T index 
      edgeListTmp[newNNZ].weight = weight; //The value
      edgeListPtr[Si+1]++;
      edgeListPtr[Ti+1]++;
      newNNZ++;
    }//End of Else
  }//End of for loop
  fclose(file); //Close the file
  //newNNZ = newNNZ / 2;
  printf("Done reading from file.\n");
  printf("Modified the number of edges from %ld ", NE);
  NE = newNNZ; //#NNZ might change
  printf("to %ld \n", NE);

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
  time1 = omp_get_wtime();
  edge *edgeList = (edge *) malloc( 2*NE * sizeof(edge)); //Every edge stored twice
  assert(edgeList != 0);
  //Keep track of how many edges have been added for a vertex:
  long  *Counter = (long *) malloc (NV  * sizeof(long)); assert(Counter != 0);
#pragma omp parallel for
  for (long i = 0; i < NV; i++) {
    Counter[i] = 0;
  }
  time2 = omp_get_wtime();
  printf("Time for allocating memory for edgeList = %lf\n", time2 - time1);
  printf("About to build edgeList...\n");

  time1 = omp_get_wtime();
  //Build the edgeList from edgeListTmp:
#pragma omp parallel for
  for(long i=0; i<NE; i++) {
    long head     = edgeListTmp[i].head;
    long tail     = edgeListTmp[i].tail;
    double weight = edgeListTmp[i].weight;

    long Where    = edgeListPtr[head] + __sync_fetch_and_add(&Counter[head], 1);
    edgeList[Where].head = head;
    edgeList[Where].tail = tail;
    edgeList[Where].weight = weight;
    //Now add the edge the other way:
    Where                  = edgeListPtr[tail] + __sync_fetch_and_add(&Counter[tail], 1);
    edgeList[Where].head   = tail;
    edgeList[Where].tail   = head;
    edgeList[Where].weight = weight;
  }
  time2 = omp_get_wtime();
  printf("Time for building edgeList = %lf\n", time2 - time1);

  G->sVertices    = NV;
  G->numVertices  = NV;
  G->numEdges     = NE;
  G->edgeListPtrs = edgeListPtr;
  G->edgeList     = edgeList;

  free(edgeListTmp);
  free(Counter);

}//End of parse_MatrixMarket_Sym_AsGraph()
