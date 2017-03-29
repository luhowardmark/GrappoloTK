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
//extern GraphElem MaxDegree;
/* The redistritbuted coloring step, no balance */
ColorElem reColor(const Graph &g, ColorVector &baseColors, std::string input, ColorElem ncolors,double factor)
{
	double t1,t2;
	t1 = mytimer();
	const GraphElem nv = g.getNumVertices();
	RandVec randVec(nv);
	generateRandomNumbers(randVec);

	srand(time(NULL));
	std::string outb;
	outb = input+".ReC"+std::to_string(factor);

	// Rebuild indirection for coloring
	ColorVector colorPtr, colorIndex, freq;
	colorPtr.resize(ncolors+1);
	colorIndex.resize(nv);
	freq.resize(ncolors);
	buildColorsIndex(baseColors,ncolors,nv,colorPtr,colorIndex,freq);
	GraphElem avg = ceil(nv/ncolors);

	// initialize the color -1, prepare for recolor
	ColorVector newColors(nv);
	ColorVector newFreq(MaxDegree,0);
	#pragma omp parallel for default(none), shared(newColors,baseColors), schedule(static)
	for(GraphElem i = 0L; i<nv;i++)
		newColors[i]= -1;

	// reverse the vertex from highest color to lowest color
	ColorQueue q(nv), qtmp;
	GraphElem qtmpPos = 0L;
	#pragma omp parallel for default(none), shared(q,colorIndex),schedule(static)
	for (GraphElem i = 0; i < nv; i++)
		q[i] = colorIndex[nv-1-i];

	// Conflicts check statistics
	ColorElem nconflicts=0;
	int nloops = 0;
	GraphElem realMaxDegree = -1;

	/* Cal real Maximum degree, no used
	#pragma omp parallel for default(none), shared(g), reduction(max: realMaxDegree), schedule(static)
        for (GraphElem i = 0L; i < nv; i++) {
                GraphElem e0, e1, er;
                g.getEdgeRangeForVertex(i, e0, e1);
                er = e1 - e0;
                if (er > realMaxDegree)
                        realMaxDegree = er;
        }
        static_assert(sizeof(int) == sizeof(int32_t), "int should be 32-bit in size");
        assert((realMaxDegree < INT32_MAX) && (realMaxDegree > 0L));*/


	/* Begining of Redistribution */
	std::cout << "ReColor start "<< std::endl;


	// Coloring Main Loop
	do {
                size_t qsz = q.size();
                qtmp.resize(qsz);
                double mst=0, cst=0;
                #pragma omp parallel default(none), shared(ncolors,newFreq,q, qtmp, qtmpPos, randVec, g, newColors, qsz, nloops, nconflicts, std::cerr, avg,factor)//, MaxDegree)
                {
			// Travel unprocessed overfilled vertices
			#pragma omp for firstprivate(nloops, qsz), schedule(guided)
                        for (GraphElem qi = 0L; qi < qsz; qi++) {
                                GraphElem v = q[qi];
				ColorElem maxColor = -1;
                        	BitVector mark(MaxDegree, false);
				// Mark the used color
				maxColor = distanceOneMarkArray(mark,g,v,newColors);
				ColorElem myColor=-1;
				
				//Pick the target
				for(myColor =0; myColor <MaxDegree; myColor++){
					if(mark[myColor] != true && newFreq[myColor]<(avg)){
						break;
					}
				}	
			
				if(myColor == MaxDegree){
					std::cerr<< "Increase too much color, please check the code" << std::endl;
					exit(1);
				}
                                newColors[v] = myColor;
				#pragma omp atomic update
				newFreq[myColor] ++;
			}// End of Vertex wise coloring (for)
		
			//Conflicts resloution step
			#pragma omp for firstprivate(qsz), schedule(guided)
			for (GraphElem qi = 0L; qi < qsz; qi++) {
				GraphElem v = q[qi];
				distanceOneConfResolution(g,v,newColors,newFreq,randVec,qtmp,qtmpPos,1);
			} //End of identify all conflicts (re-set conflicts to -1)
			#pragma omp single
			{
				q.resize(qtmpPos);
			}
			#pragma omp for schedule(static)
			for (GraphElem qi = 0; qi<qtmpPos;qi++)
				q[qi] = qtmp[qi];
		} //End of parallel step

		nconflicts += qtmpPos;
		nloops++;
		qtmpPos = 0;
	}while(!q.empty()); // End of the Coloring main loop
	t2 = mytimer();
	std::cout << outb << " Recoloring Time: " << t2-t1<<std::endl;

	//Sanity check;
	distanceOneChecked(g,nv,newColors);

	// Out put
        ColorElem newNcolors = -1;
        double variance;
        #pragma omp parallel for default(none), shared(newColors), reduction(max: newNcolors), schedule(static)
        for (size_t ci = 0U; ci < nv; ci++)
        {
                if (newColors[ci] > newNcolors)
                        newNcolors = newColors[ci];
        }
	ColorVector newFreq2(MaxDegree,0);
  computeBinSizes(newFreq2,newColors,nv,newNcolors);

	outPut(newColors,outb,newFreq2,newNcolors);
	/*	
	std::cout<<"Total Number of Colors used: " << newNcolors<<std::endl;
	for(int ci = 0; ci <nv; ci++)
	{
		std::cout<< ci << " : " <<newColors[ci]<<std::endl;
	}*/
	baseColors = newColors;
}

