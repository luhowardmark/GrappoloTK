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

#ifndef __ColorUtils__
#define __ColorUtils__
#include <cstdint>

#include <vector>
#include <omp.h>
#include <cassert>
#include <cstdint>

#include <algorithm>
#include <deque>
#include <iostream>
#include <sstream>
#include <vector>
#include <fstream>
#include "RngStream.h"
//##include <timer.h>
#include <ctime>
#include <cstdlib>
#include <omp.h>
#include "defs.h"



/*struct MoveInfo
{
	ColorElem source;
	ColorElem target;
	GraphElem numVertices;
	GraphElem startPost;
	
};*/

typedef std::vector<bool> BitVector;
typedef std::vector<long> ColorVector;
typedef int ColorElem;
#define MaxDegree 4096
//using namespace std;

int distanceOneMarkArray(BitVector &mark, graph *G, long v, int *vtxColor);
void computeBinSizes(ColorVector &binSizes, int* colors, long nv, int numColors);
void distanceOneConfResolution(graph* G, long v, int* vtxColor, double* randValues, long* QtmpTail, long* Qtmp, ColorVector& freq, int type);
void distanceOneChecked(graph* G, long nv ,int* colors);
void buildColorsIndex(int* colors, const int numColors, const long nv, ColorVector& colorPtr,  ColorVector& colorIndex, ColorVector& binSizes);

/******* UtiliyFunctions *****
void computeBinSizes(ColorVector &binSizes, const ColorVector &colors, const GraphElem nv, const ColorElem numColors);
ColorElem getDegree(const GraphElem ci, const Graph &g);
void computeBinSizesWeighted(ColorVector &binSizes, const ColorVector &colors, const GraphElem nv, const ColorElem numColors, const Graph &g);
void outPut(const ColorVector &colors, std::string output, const ColorVector& freq, const ColorElem ncolors);

void buildColorsIndex(const ColorVector& colors, const ColorElem numColors, const ColorElem nv, ColorVector& colorPtr,  ColorVector& colorIndex, ColorVector& binSizes);

//bool findConflicts(const Graph &g, const ColorElem targetColor, const ColorVector &colors, const GraphElem vertex);

void distanceOneConfResolution(const Graph &g, const GraphElem &sv,  ColorVector &colors, ColorVector &freq, const RandVec& randVec,  ColorQueue& qtmp, GraphElem& qtmpPos, int type);

void distanceOneConfResolutionWeighted(const Graph &g, const GraphElem &sv,  ColorVector &colors, ColorVector &freq, const RandVec& randVec,  ColorQueue& qtmp, GraphElem& qtmpPos, int type);


ColorElem distanceOneMarkArray(BitVector &mark, const Graph &g, const GraphElem &sv, const ColorVector &colors);

void distanceOneChecked(const Graph &g, const GraphElem nv ,const ColorVector &colors);
void generateRandomNumbers(std::vector<double> &randVec);

/******* Coloring Functions ******

/* Basic coloring (unbalanced) in initialColoring.cpp 
ColorElem initColoring(const Graph &g, ColorVector &colors, std::string input);
/* Basic coloiring (ab-inital) in initialColoringLU.cpp 
ColorElem initColoringLU(const Graph &g, ColorVector &colors, std::string input);

/* Vertex base redistribution in vBase.cpp
 * type: 0) FF, 1) LU			
ColorElem vBaseRedistribution(const Graph &g, ColorVector &baseColors, std::string input, ColorElem ncolors, int type);
ColorElem TrueSerialvBaseRedistribution(const Graph &g, ColorVector &baseColors, std::string input, ColorElem ncolors, int type);
ColorElem wBaseRedistribution(const Graph &g, ColorVector &baseColors, std::string input, ColorElem ncolors, int type);


ColorElem mBaseRedistribution(const Graph &g, ColorVector &baseColors, std::string input, ColorElem ncolors, int type);

/* Color base redistribution in cBase.cpp
 * type: 0) FF, 1) LU			
ColorElem cBaseRedistribution(const Graph &g, ColorVector &baseColors, std::string input, ColorElem ncolors, int type);

ColorElem reColor(const Graph &g, ColorVector &baseColors, std::string input, ColorElem ncolors,double factor);
ColorElem schRedistribution(const Graph &g, ColorVector &colors, std::string input, ColorElem ncolors);*/
#endif
