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

void parse_EdgeListBinaryNew(graph * G, char *fileName) {
  printf("Parsing a file in binary format...\n");
  printf("WARNING: Assumes that the graph is undirected -- every edge is stored twice.\n");
  int nthreads = 0;

  #pragma omp parallel
  {
    nthreads = omp_get_num_threads();
  }
  
  double time1, time2;

  std::ifstream ifs;  
  ifs.open(fileName, std::ifstream::in | std::ifstream::binary);
  if (!ifs) {
    std::cerr << "Error opening binary format file: " << fileName << std::endl;
    exit(EXIT_FAILURE);
  }

  long NV, NE, weighted;
  //Parse line-1: #Vertices #Edges
  ifs.read(reinterpret_cast<char*>(&NV), sizeof(NV));
  ifs.read(reinterpret_cast<char*>(&NE), sizeof(NE));
  ifs.read(reinterpret_cast<char*>(&weighted), sizeof(weighted));

  long* verPtrRaw = (long*) malloc( (NV+1)*sizeof(long)); assert(verPtrRaw != 0);
  edge* edgeListRaw = (edge*) malloc(2*NE*sizeof(edge)); assert(edgeListRaw != 0);

  ifs.read(reinterpret_cast<char*>(verPtrRaw), sizeof(long) * (NV+1));
  ifs.read(reinterpret_cast<char*>(edgeListRaw), sizeof(edge) * (2*NE));
 
  ifs.close(); //Close the file

  G->sVertices    = NV;
  G->numVertices  = NV;
  G->numEdges     = NE;
  G->edgeListPtrs = verPtrRaw;
  G->edgeList     = edgeListRaw;
  
  //Clean up

  //displayGraph(G);
}//End of parse_Dimacs9FormatDirectedNewD()
