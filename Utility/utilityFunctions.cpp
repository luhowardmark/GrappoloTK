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
#include "RngStream.h"
#include <algorithm>

using namespace std;

//WARNING: Assume that colorSize is populated with the frequency for each color
double computeGiniCoefficient(long *colorSize, int numColors) {

    //Step 1: Sort the color size vector -- use STL sort function
    double time1 = omp_get_wtime();
    sort(colorSize, colorSize+numColors);
    double time2 = omp_get_wtime();
    printf("Time for sorting: %g secs\n", time2-time1);
    //Step 2: Compute Gini coefficient
    double numFunc=0.0, denFunc=0.0;
    for (long i=0; i < numColors; i++) {
        numFunc = numFunc + ((i+1)*colorSize[i]);
        denFunc = denFunc + colorSize[i];
    }
    printf("Numerator = %g  Denominator = %g\n", numFunc, denFunc);
    //printf("Negative component = %g\n", ((double)(numColors+1)/(double)numColors));
    double giniCoeff = ((2*numFunc)/(numColors*denFunc)) - ((double)(numColors+1)/(double)numColors);
    
    //printf("Numerator = %g\n", ((2*numFunc)/(numColors*denFunc)));
    //printf("*** Gini coeff = %g", giniCoeff);
    
    return giniCoeff; //Return the Gini coefficient
}//End of computeGiniCoefficient()


void generateRandomNumbers(double *RandVec, long size) {
	int nT;
#pragma omp parallel
	{
		nT = omp_get_num_threads();
	}
#ifdef PRINT_DETAILED_STATS_
	printf("Within generateRandomNumbers() -- Number of threads: %d\n", nT);
#endif	
    //Initialize parallel pseudo-random number generator
	unsigned long seed[6] = {1, 2, 3, 4, 5, 6};
	RngStream::SetPackageSeed(seed);
	RngStream RngArray[nT]; //array of RngStream Objects
	
	long block = size / nT;
#ifdef PRINT_DETAILED_STATS_
	cout<<"Each thread will add "<<block<<" edges\n";
#endif
	//Each thread will generate m/nT edges each
	double start = omp_get_wtime();
#pragma omp parallel
	{
		int myRank = omp_get_thread_num();
#pragma omp for schedule(static)
		for (long i=0; i<size; i++) {
			RandVec[i] =  RngArray[myRank].RandU01();
		}
	}//End of parallel region	
} //End of generateRandomNumbers()

void displayGraph(graph *G) {
  long    NV        = G->numVertices;  
  long    NE        = G->numEdges;
  long    *vtxPtr   = G->edgeListPtrs;
  edge    *vtxInd   = G->edgeList;  
  printf("***********************************");
  printf("|V|= %ld, |E|= %ld \n", NV, NE);
  printf("***********************************");
  for (long i = 0; i < NV; i++) {
    long adj1 = vtxPtr[i];
    long adj2 = vtxPtr[i+1];
    printf("\nVtx: %ld [%ld]: ",i+1,adj2-adj1);
    for(long j=adj1; j<adj2; j++) {      
      printf("%ld (%g), ", vtxInd[j].tail+1, vtxInd[j].weight);
    }
  }
  printf("\n***********************************\n");
}

void duplicateGivenGraph(graph *Gin, graph *Gout) {
	long    NV        = Gin->numVertices;  
	long    NE        = Gin->numEdges;
	long    *vtxPtr   = Gin->edgeListPtrs;
	edge    *vtxInd   = Gin->edgeList;  
#ifdef PRINT_DETAILED_STATS_
	printf("|V|= %ld, |E|= %ld \n", NV, NE);
#endif		
	double time1 = omp_get_wtime();
	long *edgeListPtr = (long *)  malloc((NV+1) * sizeof(long));
	assert(edgeListPtr != NULL);		
#pragma omp parallel for
	for (long i=0; i<=NV; i++) {
		edgeListPtr[i] = vtxPtr[i]; //Prefix Sum
	}
	
	//WARNING: There is a bug in edge counting when self-loops exist
	edge *edgeList = (edge *) malloc( 2*NE * sizeof(edge)); 
	assert( edgeList != NULL);	
#pragma omp parallel for
	for (long i=0; i<NV; i++) {
		for (long j=vtxPtr[i]; j<vtxPtr[i+1]; j++) {
			edgeList[j].head = vtxInd[j].head;
			edgeList[j].tail = vtxInd[j].tail;
			edgeList[j].weight = vtxInd[j].weight;
		}		
	}

	//The last element of Cumulative will hold the total number of characters
	double time2 = omp_get_wtime();
#ifdef PRINT_DETAILED_STATS_
	printf("Done duplicating the graph:  %9.6lf sec.\n", time2 - time1);
#endif	
	Gout->sVertices    = NV;
	Gout->numVertices  = NV;
	Gout->numEdges     = NE;
	Gout->edgeListPtrs = edgeListPtr;
	Gout->edgeList     = edgeList;	
} //End of duplicateGivenGraph()


void displayGraphEdgeList(graph *G) {
	long    NV        = G->numVertices;  
	long    NE        = G->numEdges;
	long    *vtxPtr   = G->edgeListPtrs;
	edge    *vtxInd   = G->edgeList;  
	printf("***********************************");
	printf("|V|= %ld, |E|= %ld \n", NV, NE);
	for (long i = 0; i < NV; i++) {
		long adj1 = vtxPtr[i];
		long adj2 = vtxPtr[i+1];
		for(long j=adj1; j<adj2; j++) {      
			printf("%ld %ld %g\n", i+1, vtxInd[j].tail+1, vtxInd[j].weight);
		}
	}
	printf("\n***********************************\n");
}


void displayGraphEdgeList(graph *G, FILE* out) {
	long    NV        = G->numVertices;  
	long    NE        = G->numEdges;
	long    *vtxPtr   = G->edgeListPtrs;
	edge    *vtxInd   = G->edgeList;  
	printf("********PRINT OUTPUT********************");
	fprintf(out,"p sp %ld %ld \n", NV, NE/2);
	for (long i = 0; i < NV; i++) {
		long adj1 = vtxPtr[i];
		long adj2 = vtxPtr[i+1];
		for(long j=adj1; j<adj2; j++) {    
			if( i+1 < vtxInd[j].tail+1)
			{
				fprintf(out,"a %ld %ld %g\n", i+1, vtxInd[j].tail+1, vtxInd[j].weight);
			}
		}
	}
}

void writeEdgeListToFile(graph *G, FILE* out) {
	long    NV        = G->numVertices;
	long    *vtxPtr   = G->edgeListPtrs;
	edge    *vtxInd   = G->edgeList;
	for (long i = 0; i < NV; i++) {
		long adj1 = vtxPtr[i];
		long adj2 = vtxPtr[i+1];
		for(long j=adj1; j<adj2; j++) {    
	//		if( i < vtxInd[j].tail) {
				fprintf(out,"%ld %ld\n", i, vtxInd[j].tail);
	//		}
		}
	}
}

void displayGraphCharacteristics(graph *G) {
  printf("Within displayGraphCharacteristics()\n");        
  long    sum = 0, sum_sq = 0;
  double  average, avg_sq, variance, std_dev;
  long    maxDegree = 0;
  long    isolated  = 0;
  long    degreeOne = 0;
  long    NS        = G->sVertices;    
  long    NV        = G->numVertices;
  long    NT        = NV - NS;
  long    NE        = G->numEdges;
  long    *vtxPtr   = G->edgeListPtrs;
  long    tNV       = NV; //Number of vertices    
     
  if ( (NS == 0)||(NS == NV) ) {  //Nonbiparite graph	
    for (long i = 0; i < NV; i++) {
      long degree = vtxPtr[i+1] - vtxPtr[i];
      sum_sq += degree*degree;
      sum    += degree;
      if (degree > maxDegree)
	maxDegree = degree;
      if ( degree == 0 )
	isolated++;
      if ( degree == 1 )
	degreeOne++;
    }	
    average  = (double) sum / tNV;
    avg_sq   = (double) sum_sq / tNV;
    variance = avg_sq - (average*average);
    std_dev  = sqrt(variance);
    
    printf("*******************************************\n");
    printf("General Graph: Characteristics :\n");
    printf("*******************************************\n");
    printf("Number of vertices   :  %ld\n", NV);
    printf("Number of edges      :  %ld\n", NE);
    printf("Maximum out-degree is:  %ld\n", maxDegree);
    printf("Average out-degree is:  %lf\n",average);
    printf("Expected value of X^2:  %lf\n",avg_sq);
    printf("Variance is          :  %lf\n",variance);
    printf("Standard deviation   :  %lf\n",std_dev);
    printf("Isolated vertices    :  %ld (%3.2lf%%)\n", isolated, ((double)isolated/tNV)*100);
    printf("Degree-one vertices  :  %ld (%3.2lf%%)\n", degreeOne, ((double)degreeOne/tNV)*100);
    printf("Density              :  %lf%%\n",((double)NE/(NV*NV))*100);
    printf("*******************************************\n");
    
  }//End of nonbipartite graph
  else { //Bipartite graph
    
    //Compute characterisitcs from S side:
    for (long i = 0; i < NS; i++) {
      long degree = vtxPtr[i+1] - vtxPtr[i];
      sum_sq += degree*degree;
      sum    += degree;
      if (degree > maxDegree)
	maxDegree = degree;
      if ( degree == 0 )
	isolated++;
      if ( degree == 1 )
	degreeOne++;
    }    
    average  = (double) sum / NS;
    avg_sq   = (double) sum_sq / NS;
    variance = avg_sq - (average*average);
    std_dev  = sqrt(variance);
    
    printf("*******************************************\n");
    printf("Bipartite Graph: Characteristics of S:\n");
    printf("*******************************************\n");
    printf("Number of S vertices :  %ld\n", NS);
    printf("Number of T vertices :  %ld\n", NT);
    printf("Number of edges      :  %ld\n", NE);
    printf("Maximum out-degree is:  %ld\n", maxDegree);
    printf("Average out-degree is:  %lf\n",average);
    printf("Expected value of X^2:  %lf\n",avg_sq);
    printf("Variance is          :  %lf\n",variance);
    printf("Standard deviation   :  %lf\n",std_dev);
    printf("Isolated (S)vertices :  %ld (%3.2lf%%)\n", isolated, ((double)isolated/NS)*100);
    printf("Degree-one vertices  :  %ld (%3.2lf%%)\n", degreeOne, ((double)degreeOne/tNV)*100);
    printf("Density              :  %lf%%\n",((double)NE/(NS*NS))*100);
    printf("*******************************************\n");
    
    sum = 0; 
    sum_sq = 0;
    maxDegree = 0;
    isolated  = 0;
    //Compute characterisitcs from T side:
    for (long i = NS; i < NV; i++) {
      long degree = vtxPtr[i+1] - vtxPtr[i];
      sum_sq += degree*degree;
      sum    += degree;
      if (degree > maxDegree)
	maxDegree = degree;
      if ( degree == 0 )
	isolated++;
      if ( degree == 1 )
	degreeOne++;
    }
    
    average  = (double) sum / NT;
    avg_sq   = (double) sum_sq / NT;
    variance = avg_sq - (average*average);
    std_dev  = sqrt(variance);
    
    printf("Bipartite Graph: Characteristics of T:\n");	
    printf("*******************************************\n");
    printf("Number of T vertices :  %ld\n", NT);
    printf("Number of S vertices :  %ld\n", NS);
    printf("Number of edges      :  %ld\n", NE);
    printf("Maximum out-degree is:  %ld\n", maxDegree);
    printf("Average out-degree is:  %lf\n",average);
    printf("Expected value of X^2:  %lf\n",avg_sq);
    printf("Variance is          :  %lf\n",variance);
    printf("Standard deviation   :  %lf\n",std_dev);
    printf("Isolated (T)vertices :  %ld (%3.2lf%%)\n", isolated, ((double)isolated/NT)*100);
    printf("Degree-one vertices  :  %ld (%3.2lf%%)\n", degreeOne, ((double)degreeOne/tNV)*100);
    printf("Density              :  %lf%%\n",((double)NE/(NT*NT))*100);
    printf("*******************************************\n");    
  }//End of bipartite graph    
}


//Convert a directed graph into an undirected graph:
//Parse through the directed graph and add edges in both directions
graph * convertDirected2Undirected(graph *G) {
  printf("Within convertDirected2Undirected()\n");
  int nthreads;
#pragma omp parallel
  {
    nthreads = omp_get_num_threads();
  }
    
  double time1=0, time2=0, totalTime=0;
  //Get the iterators for the graph:
  long NVer     = G->numVertices;
  long NEdge    = G->numEdges;       //Returns the correct number of edges (not twice)
  long *verPtr  = G->edgeListPtrs;   //Vertex Pointer: pointers to endV
  edge *verInd  = G->edgeList;       //Vertex Index: destination id of an edge (src -> dest)
  printf("N= %ld  NE=%ld\n", NVer, NEdge);
  
  long *degrees = (long *) malloc ((NVer+1) * sizeof(long));
  assert(degrees != NULL);
  
  //Count the edges from source --> sink (for sink > source)
#pragma omp parallel for
  for (long v=0; v < NVer; v++ ) {
    long adj1 = verPtr[v];
    long adj2 = verPtr[v+1];
    for(long k = adj1; k < adj2; k++ ) {
      long w = verInd[k].tail;
      if (w < v)
	continue;
      __sync_fetch_and_add(&degrees[v+1], 1); //Increment by one to make space for zero
      __sync_fetch_and_add(&degrees[w+1], 1);
    }//End of for(k)
  }//End of for(i)
  
  //Build the pointer array:
  long m=0;
  for (long i=1; i<=NVer; i++) {
    m += degrees[i]; //Accumulate the number of edges
    degrees[i] += degrees[i-1];
  }
  //Sanity check:
  if(degrees[NVer] != 2*m) {
    printf("Number of edges added is not correct (%ld, %ld)\n", degrees[NVer], 2*m);
    exit(1);
  }
  printf("Done building pointer array\n");
  
  //Build CSR for Undirected graph:
  long* counter = (long *) malloc (NVer * sizeof(long));
  assert(counter != NULL);
#pragma omp parallel for
  for (long i=0; i<NVer; i++)
    counter[i] = 0;
  
  //Allocate memory for Edge list:
  edge *eList = (edge *) malloc ((2*m) * sizeof (edge));
  assert(eList != NULL);
  
#pragma omp parallel for
  for (long v=0; v < NVer; v++ ) {
    long adj1 = verPtr[v];
    long adj2 = verPtr[v+1];
    for(long k = adj1; k < adj2; k++ ) {
      long w = verInd[k].tail;
      double weight = verInd[k].weight;
      if (w < v)
	continue;
      //Add edge v --> w
      long location = degrees[v] + __sync_fetch_and_add(&counter[v], 1);
      if (location >= 2*m) {
	printf("location is out of bound: %ld \n", location);
	exit(1);
      }
      eList[location].head   = v;
      eList[location].tail   = w;
      eList[location].weight = weight;
      
      //Add edge w --> v
      location = degrees[w] + __sync_fetch_and_add(&counter[w], 1);
      if (location >= 2*m) {
	printf("location is out of bound: %ld \n", location);
	exit(1);
      }
      eList[location].head   = w;
      eList[location].tail   = v;
      eList[location].weight = weight;      
    }//End of for(k)
  }//End of for(v)
  
  //Clean up:
  free(counter);

  //Build and return a graph data structure
  graph * Gnew = (graph *) malloc (sizeof(graph));
  Gnew->numVertices  = NVer;
  Gnew->sVertices    = NVer;
  Gnew->numEdges     = m;
  Gnew->edgeListPtrs = degrees;
  Gnew->edgeList     = eList;
  
  return Gnew;
  
}//End of convertDirected2Undirected()

