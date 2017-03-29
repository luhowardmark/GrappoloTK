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
#include "stdlib.h"
#include "time.h"
//Return the number of colors used (zero is a valid color)
//Algorithm: Adaptation of Luby-Jones-Plusman
//Source: http://on-demand.gputechconf.com/gtc/2012/presentations/S0332-Efficient-Graph-Matching-and-Coloring-on-GPUs.pdf
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////

void generateRandomNumbers2(double* randValues, long NVer)
{
	for(int v = 0; v<NVer; v++)
	{
		randValues[v] = (double)rand();
	}
}

int algoColoringMultiHashMaxMin(graph *G, int *vtxColor, int nThreads, double *totTime, int nHash, int nItrs)
{
	srand (time(NULL));
#ifdef PRINT_DETAILED_STATS_
    printf("Within algoColoringMultiHashMaxMin()\n");
#endif
    
    if (nThreads < 1)
        omp_set_num_threads(1); //default to one thread
    else
        omp_set_num_threads(nThreads);
    int nT;
    
#pragma omp parallel
    {
        nT = omp_get_num_threads();
    }
    
#ifdef PRINT_DETAILED_STATS_
    printf("Actual number of threads: %d (requested: %d)\n", nT, nThreads);
#endif
    
    assert(nItrs > 0);
    assert(nHash > 0);
    // NOTE: Cheating for now:
    // nHash = 1;
    double time1=0, time2=0, totalTime=0;
    //Get the iterators for the graph:
    long NVer    = G->numVertices;
    long NEdge   = G->numEdges;
    long *verPtr = G->edgeListPtrs;   //Vertex Pointer: pointers to endV
    edge *verInd = G->edgeList;       //Vertex Index: destination id of an edge (src -> dest)
    
    int maxColor = 2 * nHash * nItrs;
	int totalColored = 0;
#ifdef PRINT_DETAILED_STATS_
    printf("Vertices: %ld  Edges: %ld   Max color: %d\n", NVer, NEdge, maxColor);
#endif
    
    //Build a vector of random numbers:
    //Note: Cheating a little bit now -- need to fix this with a hash function
		// Cheating Big Time now with multiple randValues
		/*
		double *randValues = (double*) malloc (NVer * sizeof(double));
    assert(randValues != 0);*/
    
	double ** randValuesPtr = (double**) malloc(nHash* sizeof(double*));
	for(int i = 0; i < nHash; i++)
	{
		randValuesPtr[i] = (double*) malloc (NVer * sizeof(double));
		generateRandomNumbers2(randValuesPtr[i], NVer);
	}
	
	// Use this to simulate the distributed enviroment
	int *vtxColor_2 = (int *) malloc (G->numVertices * sizeof(int));
    	
    //Color all the vertices to a maximum number (means that the vertex did not get colored)
	#pragma omp parallel for
    for (long v=0; v<NVer; v++) {
        vtxColor[v] = maxColor; //Set the color to maximum
		vtxColor_2[v] = maxColor;
    }
    
    //Loop through the iterations:
    for (int itr=0; itr<nItrs; itr++) {
        //Iterate for the number of hashes
		#pragma omp parallel for
        for (long v=0; v<NVer; v++) {
			if(vtxColor[v] != maxColor)
                    continue;
			for (int ihash=0; ihash<nHash; ihash++) {
				//Iterate over all the vertices:
                //Check if this vertex has already been colored
                //Vertex v has not been colored. Check to see if it is a local max or a local min
                long adj1 = verPtr[v];
                long adj2 = verPtr[v+1];
                if(vtxColor_2[v] != maxColor)
					continue;
				//Browse the adjacency set of vertex v
                bool isMax = true, isMin = true;
                for(long k = adj1; k < adj2; k++ ) {
                    if ( v == verInd[k].tail ) //Self-loops
                        continue;
                    if(vtxColor[verInd[k].tail] < maxColor)
                        continue; //It has already been colored -- ignore this neighbor
                    if ( randValuesPtr[ihash][v] <= randValuesPtr[ihash][verInd[k].tail] ) {
                        isMax = false;
                    }
                    if ( randValuesPtr[ihash][v] >= randValuesPtr[ihash][verInd[k].tail] ) {
                        isMin = false;
                    }
                    //Corner case: if all neighbors have been colored,
                    //both isMax and isMin will be true
                }//End of for(k)
				if (isMax == true) {
                    vtxColor_2[v] = (itr * 2 * nHash) + 2*ihash;
                } else if (isMin == true) {
                    vtxColor_2[v] = (itr * 2 * nHash) + 2*ihash + 1;
                }
            }//End of for(ihash)
        } //End of for(v)
		
		/*//Update hash_function if we have multiple function to use, but wanted in the same iteration
		for(int i = 0; i < nHash; i++)
		{
			generateRandomNumbers2(randValuesPtr[i], NVer);
		}*/
		
		//Update v_color, cheaing here, it has a better way to update
		int iterFreq = 0;
		#pragma omp parallel for
		for (long v=0; v<NVer; v++)
		{
			if(vtxColor[v] != vtxColor_2[v])
			{
				__sync_fetch_and_add(&iterFreq,1);
				vtxColor[v] = vtxColor_2[v];	
			}
		}
		totalColored += iterFreq;
		printf("Iteration : %d, Num. V Colored in this iteration: %d (%lf), Number of v colored total: %d (%lf)\n", itr, iterFreq, (double)iterFreq/NVer, totalColored, (double)totalColored/NVer);
		
    } //End of for(itr)
    
#ifdef PRINT_DETAILED_STATS_
    printf("***********************************************\n");
    printf("Total number of colors used: %d \n", maxColor);
    printf("Total Time                 : %lf sec\n", totalTime);
    printf("***********************************************\n");
#endif
    *totTime = totalTime;
    
    //Cleanup:
    for(int i = 0; i < nHash; i++)
	{
	//	for(int j =0; j < NVer; j++)
	//		printf("%lf",randValuesPtr[i][j]);
		free(randValuesPtr[i]);
	}
	free(randValuesPtr);
	free(vtxColor_2);
	
    return maxColor; //Return the number of colors used
}


