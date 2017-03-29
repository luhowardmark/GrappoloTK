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
ColorElem cBaseRedistribution(graph* G, int* vtxColor, int ncolors, int type)
{
	printf("Color Base Redistribution")
	
	double time1=0, time2=0, totalTime=0;
	long NVer    = G->numVertices;
  long NEdge   = G->numEdges;  
  long *verPtr = G->edgeListPtrs;   //Vertex Pointer: pointers to endV
  edge *verInd = G->edgeList;       //Vertex Index: destination id of an edge (src -> dest)
	
	
	// Rebuild indirection for coloring
	ColorVector colorPtr, colorIndex, freq;
	colorPtr.resize(ncolors+1);
	colorIndex.resize(NVer);
	freq.resize(ncolors);
	
	buildColorsIndex(vtxColor, ncolors, NVer, colorPtr, colorIndex, freq);
	
	BitVector overSize(ncolors,false);
	long avg = ceil(nv/ncolors);
       
	// Find the overSize bucket (can do some Optimization here)
	#pragma omp parallel for
	for(int ci = 0U; ci <ncolors; ci++){
		if(freq[ci]>avg)
			overSize[ci]= true;
	}

	
	/* Begining of Redistribution */
	std::cout <<"AVG:"<<avg<< " CR start "<< std::endl;

		
	// Color Base Redist. 
	#pragma omp parallel default(none), shared(overSize, colorPtr,colorIndex,type,freq,avg,ncolors,g,colors,std::cerr,std::cout)
	{
		// Travel all colors
		for(ColorElem CI = 0; CI<ncolors && overSize[CI] ==true ;CI++){
			long cadj1	= colorPtr[CI];
			long cadj2 = colorPtr[CI+1];
			
			// Move the vetex in same bin together
			#pragma omp for schedule(guided)
			for(long ki=cadj1; ki<cadj2; ki++){
				long v = colorIndex[ki];
				
				if(freq[CI] <= avg)
					continue;
				
				//std::cout<<"Move: " <<v<<std::endl;
				// Build mark array for movement canadi.
				BitVector mark(ncolors,false);
				distanceOneMarkArray(mark,g,v,colors);
				if(colors[v]!=-1)
					mark[colors[v]] = true;	
					
				//Pick target
				ColorElem myColor = -1;
				if(type == 0){	//First Fit
					for(myColor = 0; myColor <ncolors; myColor++)
						if(mark[myColor] != true && freq[myColor]<avg)
							break;
				}else if(type == 1){ //Least Used
					for(ColorElem ci = 0; ci<ncolors; ci++)
						if(mark[ci]!=true && freq[ci] <avg)
							if(myColor == -1 || freq[myColor]>freq[ci])
								myColor = ci;
				}
				
				//Update the color
				if(myColor != -1 && myColor < ncolors){
					#pragma omp atomic update
					freq[myColor]++;
					#pragma omp atomic update
					freq[colors[v]]--;
					colors[v] = myColor;
				}
			}// End of singl color
		}//End of all colors.
	} //End of parallel step
	t2 = mytimer();
	std::cout << outb << " CRE Time: " << t2-t1<<std::endl;

	//Sanity check;
	distanceOneChecked(g,nv,colors);

	// Out put
	std::cout<<"Total Number of Colors used: " << ncolors<<std::endl;
	for(int ci = 0; ci <nv; ci++)
	{
		std::cout<< ci << " : " <<colors[ci]<<std::endl;
	}
	baseColors = colors;
}

