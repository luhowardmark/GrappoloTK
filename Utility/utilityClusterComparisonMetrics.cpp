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
#include "utilityClusteringFunctions.h"

using namespace std;

/* Merkin distance computed as described in "Comparing Clustering -- An 
 Axiomatic View" by Marina Meila. Proceedings of the 22nd International 
 Conference on Machine Learning, Bonn, Germany. 2005.  
 */

//Assume that clusters have been numbered in a contiguous manner
double computeMerkinMetric(long* C1, long N1, long* C2, long N2) {
    assert(N1>0 && N2>0);
    assert((C1 != 0) && (C2 != 0));
    //Compute number of communities in each set:
    //Assume zero is a valid community id
    long nC1=-1;
    for(long i = 0; i < N1; i++) {
        if (C1[i] > nC1) {
            nC1 = C1[i];
        }
    }
    assert(nC1>0);
    
    //STEP 1: Create a CSR-like datastructure for communities in C1
    long * commPtr1 = (long *) malloc ((nC1+1) * sizeof(long)); assert(commPtr1 != 0);
    long * commIndex1 = (long *) malloc (N1 * sizeof(long)); assert(commIndex1 != 0);
    long * commAdded1 = (long *) malloc (nC1 * sizeof(long)); assert(commAdded1 != 0);
    
    // Initialization
#pragma omp parallel for
    for(long i = 0; i < nC1; i++) {
        commPtr1[i] = 0;
        commAdded1[i] = 0;
    }
    commPtr1[nC1] = 0;
    // Count the size of each community
#pragma omp parallel for
    for(long i = 0; i < N1; i++) {
        __sync_fetch_and_add(&commPtr1[(long)C1[i]+1],1);
    }
    //Prefix sum:
    for(long i=0; i<nC1; i++) {
        commPtr1[i+1] += commPtr1[i];
    }
    //Group vertices with the same color in particular order
#pragma omp parallel for
    for (long i=0; i<N1; i++) {
        long tc = (long)C1[i];
        long Where = commPtr1[tc] + __sync_fetch_and_add(&(commAdded1[tc]), 1);
        commIndex1[Where] = i; //The vertex id
    }
    
    //Compare all pairs of vertices in each community from C1 to those in C2:
    long nDisagree = 0;
    #pragma omp parallel for
    for(long ci = 0; ci < nC1; ci++) {
        long adj1 = commPtr1[ci];
        long adj2 = commPtr1[ci+1];
        for(long i=adj1; i<adj2; i++) {
            for(long j=i+1; j<adj2; j++) {
                //Check if the two vertices belong to the same community in C2
                if(C2[commIndex1[i]] != C2[commIndex1[j]])
                    __sync_fetch_and_add(&nDisagree,1); //Increment the counter
            }//End of for(j)
        }//End of for(i)
    }//End of for(ci)
    
    double dM = (2 * nDisagree) / (N1 * N2);
    
    //Cleanup:
    free(commPtr1); free(commIndex1); free(commAdded1);
    
    return dM;
    
} //End of computeMerkinDistance()



//Assume that clusters have been numbered in a contiguous manner
double computeVanDongenMetric(long* C1, long N1, long* C2, long N2) {
    cout << "Function computeVanDongenMetric() has not been implemented.\n";
}

  
