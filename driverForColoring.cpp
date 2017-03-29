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
#include "coloring.h"
#include "input_output.h"
#include "basic_util.h"
int main(int argc, char** argv) {

   //Parse Input parameters:
  clustering_parameters opts;
  if (!opts.parse(argc, argv)) {
    return -1;
  }
  int nT = 1; //Default is one thread
#pragma omp parallel
  {
      nT = omp_get_num_threads();
  }
  if (nT <= 1) {
	printf("The number of threads should be greater than one.\n");
	return 0;
  }
  graph* G = (graph *) malloc (sizeof(graph));

  int fType = opts.ftype; //File type
  char *inFile = (char*) opts.inFile;

	if(fType == 1)
     parse_MatrixMarket_Sym_AsGraph(G, inFile);
  else if(fType == 2)
     parse_Dimacs9FormatDirectedNewD(G, inFile);
  else if(fType == 3)
     parse_PajekFormat(G, inFile);
  else if(fType == 4)
     parse_PajekFormatUndirected(G, inFile);
  else if(fType == 5)
     loadMetisFileFormat(G, inFile); 
  else if(fType == 6)
    parse_UndirectedEdgeList(G, inFile);
  else if(fType == 7)
		parse_DirectedEdgeList(G, inFile);
  else if(fType == 8)
    parse_SNAP(G, inFile);
  else if(fType == 9)
    parse_EdgeListBinaryNew(G,inFile);
  else {
    cout<<"Not a valid file type"<<endl;
    exit(1);
  }

  
  // Use -c and -y for iteration and number of hash
  //parse_EdgeListBinaryNew(G,inFile); 
  displayGraphCharacteristics(G);

  printf("Iteration: %d , NHash :%d\n", opts.coloring,opts.syncType);
  //Build vector for storing colors
  int numColors = 0;
  int *colors = (int *) malloc (G->numVertices * sizeof(int)); assert (colors != 0);

	#pragma omp parallel for
  for (long i=0; i<G->numVertices; i++) {
	colors[i] = -1;
  }
  double tmpTime;
  //numColors = algoDistanceOneVertexColoringOpt(G, colors, nT, &tmpTime);
	numColors = algoColoringMultiHashMaxMin(G, colors, 1, &tmpTime, opts.syncType, opts.coloring);
  
  printf("Time to color: %lf\n", tmpTime);
  //return 0;

  /* Step : Set up output parameters */
 /* char buf1[256];
  sprintf(buf1,"%s_colorFreq.txt",argv[1]);
  FILE* out1 = fopen(buf1,"w");
  char buf2[256];
  sprintf(buf2,"%s_colorTiming.txt",argv[1]);
  FILE* out2 = fopen(buf2,"w");
 */
 
 printf("What is this %s\n",argv[1]);
 /*
  int curThread = 2; //Start with two threads
  while (curThread <= nT) {
	printf("\n\n***************************************\n");
	printf("Starting run with %d threads.\n", curThread);
	printf("***************************************\n");
	fprintf(out2, "*******************************\n");
	fprintf(out2, "Number of threads: %d\n",curThread);
	fprintf(out2, "*******************************\n");
				
#pragma omp parallel for
	for (long i=0; i<G->numVertices; i++) {
		colors[i] = -1;
	}
	numColors = algoDistanceOneVertexColoringOpt(G, colors, curThread, &tmpTime);
	fprintf(out2, "Time to color : %lf\n",tmpTime);
		
	curThread = curThread*2;
  }//End of while()
*/	
  //Count the frequency of colors:
  long *colorFreq = (long *) malloc (numColors * sizeof(long)); assert(colorFreq != 0);
#pragma omp parallel for
  for(long i = 0; i < numColors; i++) {
	colorFreq[i] = 0;
  }
	
  #pragma omp parallel for
  for(long i = 0; i < G->numVertices; i++) {
	__sync_fetch_and_add(&colorFreq[colors[i]],1);
  }
	
  for(long i=0; i < numColors; i++) {
 	printf("%ld \t %ld\n", i, colorFreq[i]);
  }
	
  /* Step 5: Clean up */
  free(G->edgeListPtrs);
  free(G->edgeList);
  free(G);
	
  free(colorFreq);
  free(colors);
  
  return 0;

}
