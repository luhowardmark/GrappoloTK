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
void writeGraphMetisSimpleFormat(graph* G, char *filename) {
  //Get the iterators for the graph:
  long NVer     = G->numVertices;
  long NEdge    = G->numEdges;       //Returns the correct number of edges (not twice)
  long *verPtr  = G->edgeListPtrs;   //Vertex Pointer: pointers to endV
  edge *verInd = G->edgeList;       //Vertex Index: destination id of an edge (src -> dest)
  printf("NVer= %ld --  NE=%ld\n", NVer, NEdge);

  printf("Writing graph in Metis format - each edge represented twice -- no weights; 1-based indices\n");
  printf("Graph will be stored in file: %s\n", filename);
  
   
  FILE *fout;
  fout = fopen(filename, "w");
  if (!fout) {
    printf("Could not open the file \n");
    exit(1);
  }
  //First Line: #Vertices #Edges
  fprintf(fout, "%ld %ld\n", NVer, NEdge);
  //Write the edges:
  for (long v=0; v<NVer; v++) {
    long adj1 = verPtr[v];
    long adj2 = verPtr[v+1];
    for(long k = adj1; k < adj2; k++ ) {
      fprintf(fout, "%ld ", (verInd[k].tail+1) );
    }
    fprintf(fout, "\n");
  }
  fclose(fout);
  printf("Graph has been stored in file: %s\n",filename);
}//End of writeGraphPajekFormat()