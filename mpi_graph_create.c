#include <mpi.h>

int MPIX_GRAPH_CREATE(MPI_Comm comm_old,
                      int nnodes, const int index[], const int edges[],
                      int reorder, MPI_Comm *comm_graph) {
    int retVal = MPI_SUCCESS;

    int myRank;
    MPI_Comm_rank(comm_old, &myRank);
    int base = myRank ? index[myRank-1] : 0;

    // first, look for others in my edges
    int outdegree = myRank ? index[myRank]-index[myRank-1] : index[myRank];
    int destinations[outdegree];
    int destweights[outdegree];
    for (int i=0;i<outdegree;++i) {
        destinations[i] = edges[base+i];
        destweights[i] = 0;
    }
    
    // Direction of edges is an arbitrary choice. 
    // For neighbourhood collectives to be legal,
    // the adjacency matrix must be symmetric and
    // this must produce the same answer as above
    // but it takes longer and it is more complex
    
    // next, look for me in others' edges
    int indegree = 0, otherRank = 0;
    int sources[index[nnodes-1]];
    int sourceweights[index[nnodes-1]];
    for (int i=0;i<index[nnodes-1];++i) {
        if (edges[i]!=myRank) continue;
        ++indegree;
        while (index[otherRank]-1<i) ++otherRank;
        sources[indegree] = otherRank;
        sourceweights[indegree] = 0;
    }

    // optional: create a custom info key & value
    // for later identification by this extension
    MPI_Info info;
    MPI_Info_create(&info);
    MPI_Info_set(info, "MPIX_TOPOL_TYPE", "GRAPH");

    retVal = MPI_Dist_graph_create_adjacent(comm_old,
                 indegree, sources, sourceweights,
                 outdegree, destinations, destweights,
                 info, reorder, comm_graph);

    return retVal;
}
