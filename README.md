Grappolo-TK: The Grappolo graph toolkit 

This graph toolkit contains multithreaded C++ implementations for graph community detection and balanced graph coloring. 

Community detection: Community detection is a widely-used operation in graph theory. The goal is to partition the vertex set of an input graph into tightly-knit groups or “communities", such that vertices that are assigned to a community have a higher density of edges among them than to the vertices in the rest of the network.  Grappolo (Lu et al. ParCo 2015) is a parallelization of the Louvain heuristic (Blondel et al. 2008), which is a widely used sequential implementation for community detection. Grappolo deploys a combination of heuristics - these include: coloring and vertex following for preprocessing, and minimum labelling during iterative processing. In particular, coloring helps obtain a partial ordering of vertices that aids in improving both the parallel performance and the serial convergence rate of the heuristic. 

Balanced graph coloring: Distance-1 coloring is a color assignment to the vertices of a graph such that no two adjacent vertices (i.e., connected by an edge) are assigned the same color. Traditional methods for coloring try to minimize the number of colors. However, in the context of many parallel processing applications, it also becomes important to obtain a balanced distribution of the color sizes. In this package, we provide various heuristics for obtaining a balanced coloring. These heuristics are described in (Lu et al. IPDPS'15, Lu et al. TPDS'17). The default variant of balanced coloring that is supported is "Vertex First-Fit (VFF)".  The balanced coloring code is multithreaded.


PAPER CITATIONS

If you use the Grappolo-TK implementations, please cite the following papers:

(Community Detection)

H. Lu, M. Halappanavar, A. Kalyanaraman. Parallel heuristics for scalable community detection. Parallel Computing, vol. 47, pp. 19-37, 2015, DOI: 10.1016/j.parco.2015.03.003.

(Balanced Graph Coloring)

H. Lu, M. Halappanavar, D. Chavarria-Miranda, A. Gebremedhin, A. Panyala, A. Kalyanaraman. Algorithms for Balanced Colorings with Applications in Parallel Computing. IEEE Transactions on Parallel and Distributed Systems, vol. 28, no. 5, pp. 1240-1256, May 1 2017, DOI: 10.1109/TPDS.2016,2620142. 
