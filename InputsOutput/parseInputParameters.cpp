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

using namespace std;

clustering_parameters::clustering_parameters()
: ftype(7), strongScaling(false), output(false), VF(false), coloring(0), syncType(0),
threadsOpt(false), basicOpt(0), C_thresh(0.01), minGraphSize(100000), threshold(0.000001)
{}

void clustering_parameters::usage() {
    cout << "***************************************************************************************"<< endl;
    cout << "Basic usage: Driver <Options> FileName\n";
    cout << "***************************************************************************************"<< endl;
    cout << "Input Options: \n";
    cout << "***************************************************************************************"<< endl;
    cout << "File-type  : -f <1-8>   -- default=7" << endl;
    cout << "File-Type  : (1) Matrix-Market  (2) DIMACS#9 (3) Pajek (each edge once) (4) Pajek (twice) \n";
    cout << "           : (5) Metis (DIMACS#10) (6) Simple edge list twice (7) Binary format (8) SNAP\n";
    cout << "--------------------------------------------------------------------------------------" << endl;
    cout << "Strong scaling : -s   [default=false]							" << endl;
    cout << "VF             : -v   [default=false]							" << endl;
    cout << "Output         : -o   [default=false]							" << endl;
    cout << "Coloring       : -c   [default=0]   							" << endl;
    cout << "BasicOpt       : -b   [default=0]  (0) basic (1) replaceMap    " << endl;
    cout << "syncType       : -y   [default=0]  (1) FullSync (2) NeighborSync (3) EarlyTerm (4) 1+3   " << endl;
    cout << "--------------------------------------------------------------------------------------" << endl;
    cout << "Min-size       : -m <value> -- default=100000" << endl;
    cout << "C-threshold    : -d <value> -- default=0.01" << endl;
    cout << "Threshold      : -t <value> -- default=0.000001" << endl;
    cout << "***************************************************************************************"<< endl;
}//end of usage()

bool clustering_parameters::parse(int argc, char *argv[]) {
    static const char *opt_string = "c:b:y:svof:t:d:m:";
    int opt = getopt(argc, argv, opt_string);
    while (opt != -1) {
        switch (opt) {
                
            case 'c': coloring = atol(optarg); break;
            case 'y': syncType = atol(optarg); break;
            case 'b': basicOpt = atol(optarg); break;
            case 's': strongScaling = true; break;
            case 'v': VF = true; break;
            case 'o': output = true; break;
                
            case 'f': ftype = atoi(optarg);
                if((ftype >10)||(ftype<0)) {
                    cout << "ftype must be an integer between 1 to 8" << endl;
                    return false;
                }
                break;
                
            case 't': threshold = atof(optarg);
                if (threshold < 0.0) {
                    cout << "Threshold must be non-negative" << endl;
                    return false;
                }
                break;
                
            case 'd': C_thresh = atof(optarg);
                if (C_thresh < 0.0) {
                    cout << "Threshold must be non-negative" << endl;
                    return false;
                }
                break;
                
            case 'm': minGraphSize = atol(optarg);
                if(minGraphSize <0) {
                    cout << "minGraphSize must be non-negative" << endl;
                    return false;
                }
                break;
                
            default:
                cerr << "unknown argument" << endl;
                return false;
        }
        opt = getopt(argc, argv, opt_string);
    }
    
    if (argc - optind != 1) {
        cout << "Problem name not specified.  Exiting." << endl;
        usage();
        return false;
    } else {
        inFile = argv[optind];
    }
    
#ifdef PRINT_DETAILED_STATS_
    cout << "********************************************"<< endl;
    cout << "Input Parameters: \n";
    cout << "********************************************"<< endl;
    cout << "Input File: " << inFile << endl;
    cout << "File type : " << ftype  << endl;
    cout << "Threshold  : " << threshold << endl;
    cout << "C-threshold: " << C_thresh << endl;
    cout << "Min-size   : " << minGraphSize << endl;
    cout << "basicOpt   : " << basicOpt << endl;
    cout << "SyncType   : " << syncType << endl;
    cout << "--------------------------------------------" << endl;
    if (coloring)
        cout << "Coloring   : TRUE" << endl;
    else
        cout << "Coloring   : FALSE" << endl;
    if (basicOpt)
        cout << "Replace map : TRUE" << endl;
    else
        cout << "Replace map : FALSE" << endl;
    if (strongScaling)
        cout << "Strong scaling : TRUE" << endl; 
    else
        cout << "Strong scaling : FALSE" << endl; 
    if(VF)
        cout << "VF         : TRUE" << endl;
    else
        cout << "VF         : FLASE" << endl;
    if(output)
        cout << "Output     : TRUE"  << endl;
    else
        cout << "Output     : FALSE"  << endl;    
    cout << "********************************************"<< endl;    
#endif
    
    return true;
}
