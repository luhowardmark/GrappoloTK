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
#include "input_output.h"
#include "basic_comm.h"
#include "basic_util.h"
#include "utilityClusteringFunctions.h"
#include "color_comm.h"
#include "sync_comm.h"

using namespace std;

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
    
    // File Loading
    double time1, time2;
    graph* G = (graph *) malloc (sizeof(graph));
    
    /* Step 2: Parse the graph in Matrix Market format */
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
    
    
    
    displayGraphCharacteristics(G);
    int threadsOpt = 0;
    if(opts.threadsOpt)
        threadsOpt =1;
    threadsOpt =1;
    
    /* Vertex Following option */
    if( opts.VF ) {
        printf("Vertex following is enabled.\n");
        time1 = omp_get_wtime();
        long numVtxToFix = 0; //Default zero
        long *C = (long *) malloc (G->numVertices * sizeof(long)); assert(C != 0);
        numVtxToFix = vertexFollowing(G,C); //Find vertices that follow other vertices
        if( numVtxToFix > 0) {  //Need to fix things: build a new graph
            printf("Graph will be modified -- %ld vertices need to be fixed.\n", numVtxToFix);
            graph *Gnew = (graph *) malloc (sizeof(graph));
            long numClusters = renumberClustersContiguously(C, G->numVertices);
            buildNewGraphVF(G, Gnew, C, numClusters);
            //Get rid of the old graph and store the new graph
            free(G->edgeListPtrs);
            free(G->edgeList);
            free(G);
            G = Gnew;
        }
        free(C); //Free up memory
        printf("Graph after modifications:\n");
        displayGraphCharacteristics(G);
    }//End of if( VF == 1 )
    
	   
    // Datastructures to store clustering information
    long NV = G->numVertices;
    long *C_orig = (long *) malloc (NV * sizeof(long)); assert(C_orig != 0);
    
    //Call the clustering algorithm:
    if ( opts.strongScaling ) { //Strong scaling enabled
        /* Unsupprted right now !!!!! Need to be fix !!!! */
        //Retain the original copy of the graph:
        /*	graph* G_original = (graph *) malloc (sizeof(graph)); //The original version of the graph
         time1 = omp_get_wtime();
         duplicateGivenGraph(G, G_original);
         time2 = omp_get_wtime();
         printf("Time to duplicate : %lf\n", time2-time1);
         
         //Run the algorithm in powers of two for the maximum number of threads available
         int curThread = 2; //Start with two threads
         while (curThread <= nT) {
         printf("\n\n***************************************\n");
         printf("Starting run with %d threads.\n", curThread);
         printf("***************************************\n");
         //Call the clustering algorithm:
         #pragma omp parallel for
         for (long i=0; i<G->numVertices; i++) {
         C_orig[i] = -1;
         }
         runMultiPhaseLouvainAlgorithm(G, C_orig, coloring, replaceMap,opts.minGraphSize, opts.threshold, opts.C_thresh, curThread,threadsOpt);
         //Increment thread and revert back to original graph
         if (curThread < nT) {
         //Skip copying at the very end
         //Old graph is already destroyed in the above function
         G = (graph *) malloc (sizeof(graph)); //Allocate new space
         duplicateGivenGraph(G_original, G); //Copy the original graph to G
         }
         curThread = curThread*2; //Increment by powers of two
         }//End of while() */
    }
    
    else { //No strong scaling -- run once with max threads
#pragma omp parallel for
        for (long i=0; i<NV; i++) {
            C_orig[i] = -1;
        }
        
        //runMultiPhaseLouvainAlgorithm(G, C_orig, coloring, replaceMap, opts.minGraphSize, opts.threshold, opts.C_thresh, nT,threadsOpt);
        // Change to each sub function that belong to the folder
        if(opts.coloring != 0){
            runMultiPhaseColoring(G, C_orig, opts.coloring, opts.minGraphSize, opts.threshold, opts.C_thresh, nT,threadsOpt);
        }else if(opts.syncType != 0){
            runMultiPhaseSyncType(G, C_orig, opts.syncType, opts.minGraphSize, opts.threshold, opts.C_thresh, nT,threadsOpt);
        }else{
            runMultiPhaseBasic(G, C_orig, opts.basicOpt, opts.minGraphSize, opts.threshold, opts.C_thresh, nT,threadsOpt);
        }
        
    }
    
    //Check if cluster ids need to be written to a file:
    if( opts.output ) {
        char outFile[256];
        sprintf(outFile,"%s_clustInfo", opts.inFile);
        printf("Cluster information will be stored in file: %s\n", outFile);
        FILE* out = fopen(outFile,"w");
        for(long i = 0; i<NV;i++) {
            fprintf(out,"%ld\n",C_orig[i]);
        }		
        fclose(out);
    }
    
    //Cleanup:
    if(C_orig != 0) free(C_orig);
    //Do not free G here -- it will be done in another routine.
    
    return 0;
}//End of main()
