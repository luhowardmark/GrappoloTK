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

#ifndef _DEFS_H
#define _DEFS_H

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <omp.h>
#include <iostream>
#include <time.h>
#include <fstream>
#include <map>
#include <vector>
#include <unistd.h> //For getopts()

#define MilanRealMax HUGE_VAL       // +INFINITY
#define MilanRealMin -MilanRealMax  // -INFINITY

#define PRINT_DETAILED_STATS_
//#define PRINT_TERSE_STATS_

typedef struct comm
{
  long size;
  double degree;
}Comm;

typedef struct
{
    long cid;       //community ID
    double Counter; //Weight relative to that community
} mapElement;

typedef struct /* the edge data structure */
{
  long head;
  long tail;
  double weight;
} edge;

typedef struct /* the graph data structure */
{
  long numVertices;        /* Number of columns                                */
  long sVertices;          /* Number of rows: Bipartite graph: number of S vertices; T = N - S */
  long numEdges;           /* Each edge stored twice, but counted once        */
  long * edgeListPtrs;     /* start vertex of edge, sorted, primary key        */
  edge * edgeList;         /* end   vertex of edge, sorted, secondary key      */
} graph;

struct clustering_parameters 
{
  const char *inFile; //Input file
  int ftype;  //File type

  bool strongScaling; //Enable strong scaling
  bool output; //Printout the clustering data
  bool VF; //Vertex following turned on
  int coloring; // Type of coloring
  int syncType; // Type of synchronization method
  int basicOpt; //If map data structure is replaced with a vector
  bool threadsOpt;
  double C_thresh; //Threshold with coloring on
  long minGraphSize; //Min |V| to enable coloring
  double threshold; //Value of threshold
       
  clustering_parameters();
  void usage();    
  
  //Define in parseInputParameter.cpp
  bool parse(int argc, char *argv[]);
};

#endif
