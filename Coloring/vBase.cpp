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

#include "coloringUtils.h"
#include "defs.h"
#include "coloring.h"

/* The redistritbuted coloring step, no balance */
int vBaseRedistribution(graph* G, int* vtxColor, int ncolors, int type)
{
#ifdef PRINT_DETAILED_STATS_
  printf("Vertex base redistribution\n");
#endif
	
  double time1=0, time2=0, totalTime=0;
  //Get the iterators for the graph:
  long NVer    = G->numVertices;
  long NEdge   = G->numEdges;  
  long *verPtr = G->edgeListPtrs;   //Vertex Pointer: pointers to endV
  edge *verInd = G->edgeList;       //Vertex Index: destination id of an edge (src -> dest)

#ifdef PRINT_DETAILED_STATS_
  printf("Vertices: %ld  Edges: %ld\n", NVer, NEdge);
#endif

  //Build a vector of random numbers
  double *randValues = (double*) malloc (NVer * sizeof(double));
  assert(randValues != 0);
  generateRandomNumbers(randValues, NVer);

	long *Q    = (long *) malloc (NVer * sizeof(long)); assert(Q != 0);
  long *Qtmp = (long *) malloc (NVer * sizeof(long)); assert(Qtmp != 0);
  long *Qswap;    
  if( (Q == NULL) || (Qtmp == NULL) ) {
    printf("Not enough memory to allocate for the two queues \n");
    exit(1);
  }
	
	
	// initialize the color to baseColor
	int *baseColors = (int *) malloc (NVer * sizeof(int)); assert (baseColors != 0);
	#pragma omp parallel for
	for(long i = 0; i<NVer;i++)
		baseColors[i]=vtxColor[i];

	// Put uncolor vertices in the queue
	long QTail=0;    //Tail of the queue 
  long QtmpTail=0; //Tail of the queue (implicitly will represent the size)
  long realMaxDegree = 0;
	
	#pragma omp parallel for
  for (long i=0; i<NVer; i++) {
      Q[i]= i;     //Natural order
      Qtmp[i]= -1; //Empty queue
  }
  QTail = NVer;	//Queue all vertices

	
	// Cal real Maximum degree, no used
	#pragma omp parallel for reduction(max: realMaxDegree)
	for (long i = 0; i < NVer; i++) {
		long adj1, adj2, de;
		adj1 = verPtr[i];
		adj2 = verPtr[i+1];
		de = adj2-adj1;
		if ( de > realMaxDegree)
			realMaxDegree = de;
	}
	  
	/////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////// START THE WHILE LOOP ///////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////
  long nConflicts = 0; //Number of conflicts 
  int nLoops = 0;     //Number of rounds of conflict resolution

	// Holder for frequency, could use realMaxDegree here
	ColorVector freq(ncolors,0);
	BitVector overSize(ncolors,false);
	long avg = (long)ceil((double)NVer/(double)ncolors);

	// calculate the frequency 
	computeBinSizes(freq,baseColors,NVer,ncolors);
	
	// Find the overSize bucket (can do some Optimization here)
	#pragma omp parallel for
	for(size_t ci = 0U; ci <ncolors; ci++)
		if(freq[ci]>avg)
			overSize[ci]= true;

	/* Begining of Redistribution */
	std::cout << "VR start "<< std::endl;


	// Coloring Main Loop
	do{
		time1 = omp_get_wtime();
		#pragma omp parallel for
    for (long Qi=0; Qi<QTail; Qi++) {
      long v = Q[Qi]; //Q.pop_front();
			int maxColor = 0;
			
			if( overSize[baseColors[v]] == false)
				continue;
			if( (vtxColor[v] != -1) && (freq[vtxColor[v]] <= avg))
				continue;
			
			BitVector mark(MaxDegree, false);
			maxColor = distanceOneMarkArray(mark,G,v,vtxColor);
			
			int myColor = -1;
			int permissable = 0;
			
			if(type == 0){	// First Fit
				for (myColor=0; myColor<=ncolors; myColor++) {
					if ( (mark[myColor] == false) && (freq[myColor]<avg) && (overSize[myColor]!= true))
						break;
				}
			}
			else if(type == 1){ // Least use
				for(int ci = 0; ci<ncolors;ci++){
					if(mark[ci] != true && freq[ci]<avg && overSize[ci]!=true){
						if(myColor==-1||freq[myColor]>freq[ci]){
							myColor = ci;
						}
					}
				}
			}
			
			if(vtxColor[v]==-1 && (myColor==-1 || myColor ==ncolors) )
				myColor=baseColors[v];
			
			if(myColor != ncolors && myColor !=-1){
				#pragma omp atomic update
				freq[myColor]++;
				if(vtxColor[v] != -1){
					#pragma omp atomic update
					freq[vtxColor[v]]--;
				}
				vtxColor[v] = myColor;
			}
		}	// End of vertex wise redistribution

		time2 = omp_get_wtime();
		
		#pragma omp parallel for
		for (long Qi=0; Qi<QTail; Qi++) {
			long v = Q[Qi]; //Q.pop_front();
			distanceOneConfResolution(G, v, vtxColor, randValues, &QtmpTail, Qtmp, freq, 1);
		} //End of outer for loop: for each vertex
  
		time2  = omp_get_wtime() - time2;
		totalTime += time2;    
		nConflicts += QtmpTail;
		nLoops++;

#ifdef PRINT_DETAILED_STATS_
    printf("Num conflicts      : %ld \n", QtmpTail);
    printf("Time for detection : %lf sec\n", time2);
#endif

    //Swap the two queues:
    Qswap = Q;
    Q = Qtmp; //Q now points to the second vector
    Qtmp = Qswap;
    QTail = QtmpTail; //Number of elements
    QtmpTail = 0; //Symbolic emptying of the second queue    
  } while (QTail > 0);

	//Sanity check;
	distanceOneChecked(G,NVer,vtxColor);
}

