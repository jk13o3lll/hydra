// g++ -D NETWORK1 test.cpp -o test.exe

#include <stdio.h>
#include <stdlib.h>
#include <queue>
#include "data.h"

// Notes 1:
// Pressure loss by Darcy-Weisbach equation (dP = r * Q ^ n), and form loop equations
// Net node flowrate by mass conservation, and form node equations

// Notes 2:
// For bc is flowrate to node (bcType = 0), E = (N - 1) + L
// (number of unknown = E, number of node equation = N - 1, number of loop equation = L)
// for bc is pressure at node (bcType = 1), E = (N - N0) + (L + (N0 - 1))
// (number of unknown = E, number of node equation = N - N0, number of loop equation = L + N0 - 1)
// (p.s. N0 - 1 loops equation by connecting one source node to others alternately)
// E: number of flow (= number of edge)
// L: number of independent loop (not conjugated to each other)
// N: number of nodes (including source nodes, number of source nodes is N0)

// 2D index to 1D: r = row, c = column, cols = number of columns
#define RC2I(r,c,cols) ((r)*(cols)+(c))

struct Adjacency {  // id of flow on that edge
    int id, state;  // -1 for direction opposite to defined flow direction,
};                  // 0 for no connnection, +1 for direction indentical to defined flow direction.

// compute residual, R = A(x) - b
void computeR(const double *incNodes, const double *conNodes,   // node eqs
              const double *incLoops, const double *conLoops,   // loop eqs
              const double *x,                                  // vector of flows
              const int nNeq, const int nLeq,                   // dims
              double *R){                                       // residual (output)
    // ...
}

// compute jacobian of R
void computeJ(const double *incNodes, const double *conNodes,   // node eqs
              const double *incLoops, const double *conLoops,   // loop eqs
              const double *x,                                  // vector of flows
              const int nNeq, const int nLeq,                   // dims
              double *J){                                       // jacobian (output)
    // ...
}

// Calculate dx by BiCG or BiCGSTAB instead of calculating inverve Jacobian
void bicg(const double *A, double *x, const double *b, const int n){ // A should be square (can be nonsymmetric)
    // ...
}

// main function
int main(int argc, char *argv[]){
    const int nE = sizeof(edgeList) / sizeof(Edge);         // number of edges
    const int nN0 = sizeof(srcList) / sizeof(SourceNode);   // number of source node
    int nN, nL, nNeq, nLeq;                                 // number of nodes, loops, node eqns, loop eqns
    Adjacency *adj;                 // adjacency matrix for pipes
    double *incNodes, *conNodes;    // incidence matrix for node eqs, constant at the same side with unknowns
    double *incLoops, *conLoops;    // incidence matrix for loop eqs, constant at the same side with unknowns
    double *inlets;                 // inlet flow on every nodes
    int root, *parent, *depth;      // spanning tree info
    bool *hasNeq;                   // mark nodes to construct node eq
    bool *visited;                  // record visited nodes for BFS
    bool *visitedEdges;             // record visited edges while BFS
    std::queue<int> buffer;         // BFS buffer
    int i, j, k, ii, jj, kk, tmp;

    // TODO: add function to modify network (replace some pipe with equivalent pipe)
    // 1. parallel pipes (> 2 pipe between 2 nodes)
    // 2. serial pipes (node deg = 2)
    // 3. multiple source flow on one node
    // analyse modified network, and use result of modified network to derive result of original network

    // get dimensions
    for(i = 0; i < nE; ++i){
        if(edgeList[i].a > nN)  nN = edgeList[i].a;
        if(edgeList[i].b > nN)  nN = edgeList[i].b;
    }
    ++nN;
    nL = nE - nN + 1;
    if(bcType == 0){
        nNeq = nN - 1;          // ignore one dependent source
        nLeq = nE - nNeq;
    }
    else if(bcType == 1){
        nNeq = nN - nN0;        // ignore all source nodes
        nLeq = nL + nN0 - 1;    // include loops for source nodes
    }
    printf("nE = %d, nN = %d, nN0 = %d, nL = %d\n", nE, nN, nN0, nL);
    printf("nNeq = %d, nLeq = %d\n", nNeq, nLeq);

    // allocate memory
    adj = (Adjacency*) malloc(nN * nN * sizeof(Adjacency));
    incNodes = (double*) malloc(nNeq * nE * sizeof(double));
    conNodes = (double*) malloc(nNeq * sizeof(double));
    incLoops = (double*) malloc(nLeq * nE * sizeof(double));
    conLoops = (double*) malloc(nLeq * sizeof(double));
    inlets = (double*) malloc(nN * sizeof(double));
    parent = (int*) malloc(nN * sizeof(int));
    depth = (int*) malloc(nN * sizeof(int));
    visited = (bool*) malloc(nN * sizeof(bool));
    visitedEdges = (bool*) malloc(nE * sizeof(bool));

    // construct adjacency matrix
    for(i = 0; i < nE; ++i)
        for(j = i; j < nE; ++j)
            adj[RC2I(i, j, nN)].state = adj[RC2I(j, i, nN)].state = 0;
    for(i = 0; i < nE; ++i){
        ii = RC2I(edgeList[i].a, edgeList[i].b, nN);
        jj = RC2I(edgeList[i].b, edgeList[i].a, nN);
        adj[jj].id = (adj[ii].id = i);
        adj[jj].state = -(adj[ii].state = 1);
    }

    // // construct incidence matrix (node equations)
    // for(i = 0; i < nN; ++i) hasNeq[i] = true, inlets[i] = 0.;
    // for(i = 0; i < nN0; ++i) inlets[srcList[i].id] = srcList[i].x;
    // if(bcType != 0){ // bc is presure
    //     for(i = 0; i < nN0; ++i) hasNeq[srcList[i].id] = false; // exclude all source nodes
    //     for(tmp = 0, i = 0; i < nN; ++i)
    //         if(hasNeq[i]){
    //             for(j = 0; j < nN; ++j){
    //                 ii = RC2I(i, j, nN);
    //                 if(adj[ii].state != 0)
    //                     incNodes[RC2I(tmp, adj[ii].id, nE)] = i == edgeList[adj[ii].id].b? 1. : -1.; // inlet as positive
    //             }
    //             conNodes[tmp] = 0.;
    //             ++tmp;
    //         }
    // }
    // else{   // bcType == 0, bc is flowrate
    //     hasNeq[srcList[nN0 - 1].id] = false; // only remove one
    //     for(tmp = 0, i = 0; i < nN; ++i)
    //         if(hasNeq[i]){
    //             for(j = 0; j < nN; ++j){
    //                 ii = RC2I(i, j, nN);
    //                 if(adj[ii].state != 0)
    //                     incNodes[RC2I(tmp, adj[ii].id, nE)] = i == edgeList[adj[ii].id].b? 1. : -1.; // inlet as positive
    //             }
    //             conNodes[tmp] = inlets[i];
    //             ++tmp;
    //         }
    // }

    // construct incidence matrix (loop equations)
    // find node with max degree as root
    ii = -1; // max degree
    for(i = 0; i < nN; ++i){
        for(tmp = 0, j = 0; j < nN; ++j)
            if(adj[RC2I(i, j, nN)].state != 0)
                ++tmp;
        if(tmp > ii) ii = tmp, root = i;
    }
    // BFS and get spanning tree
    for(i = 0; i < nN; ++i) visited[i] = false;
    for(i = 0; i < nE; ++i) visitedEdges[i] = false;
    parent[root] = -1, depth[root] = 0, visited[root] = true, buffer.push(root);
    while(!buffer.empty()){
        i = buffer.front(); buffer.pop();
        for(j = 0; j < nN; ++j){
            ii = RC2I(i, j, nN);
            if(adj[ii].state != 0 && !visited[j]){
                parent[j] = i, depth[j] = depth[i] + 1, visited[j] = true, buffer.push(i);
                visitedEdges[adj[ii].id] = true;
            }
        }
    }
    // for(i = 0; i < nN; ++i)
    //     printf("parent of %d: %d\n", i, parent[i]);

    // form loop by connecting edges not in tree
    for(tmp = 0, i = 0; i < nE; ++i)
        if(!visitedEdges[i]){
            for(ii = edgeList[i].a, jj = edgeList[i].b; ii != jj; )
                if(depth[ii] > depth[jj]){
                    kk = RC2I(parent[ii], ii, nN); // since using dir of a->b, dir is parent to child
                    k = adj[kk].id;
                    incLoops[RC2I(tmp, k, nE)] = adj[kk].state == 1? edgeList[k].r : -edgeList[k].r;
                    ii = parent[ii];
                }
                else{   // depth[ii] <= depth[jj]
                    kk = RC2I(jj, parent[jj], nN);
                    k = adj[kk].id;
                    incLoops[RC2I(tmp, k, nE)] = adj[kk].state == 1? edgeList[k].r : -edgeList[k].r;
                    jj = parent[jj];
                }
            incLoops[RC2I(tmp, i, nE)] = edgeList[i].r;
            conLoops[tmp] = 0.; // pressure loss == 0
            ++tmp;
        }
    // form source loop by connecting src node
    if(bcType == 1){
        for(i = 1; i < nN0; ++i){
            for(ii = srcList[0].id, jj = srcList[i].id; ii != jj; )
                if(depth[ii] > depth[jj]){
                    kk = RC2I(parent[ii], ii, nN); // since using dir of 0->jj0, dir is parent to child
                    k = adj[kk].id;
                    incLoops[RC2I(tmp, k, nE)] = adj[kk].state == 1? edgeList[k].r : -edgeList[k].r;
                    ii = parent[ii];
                }
                else{   // depth[ii] <= depth[jj]
                    kk = RC2I(jj, parent[jj], nN);
                    k = adj[kk].id;
                    incLoops[RC2I(tmp, k, nE)] = adj[kk].state == 1? edgeList[k].r : -edgeList[k].r;
                    jj = parent[jj];
                }
            conLoops[tmp] = srcList[i].x - srcList[0].x; // dir is 0->jj0
            ++tmp;
        }
    }

    // print eqns
    // for(i = 0; i < nNeq; ++i){
    //     printf("Node eq %d: ", i);
    //     for(j = 0; j < nE; ++j)
    //         printf("%lf ", incNodes[RC2I(i, j, nE)]);
    //     printf(", const = %lf\n", conNodes[i]);
    // }
    for(i = 0; i < nLeq; ++i){
        printf("Loop eq %d: ", i);
        for(j = 0; j < nE; ++j)
            printf("%lf ", incLoops[RC2I(i, j, nE)]);
        printf(", const = %lf\n", conLoops[i]);
    }

    // release memory not necessary in the following
    free(adj);
    free(inlets);
    free(parent);
    free(depth);
    free(visited);
    free(visitedEdges);

    // allocate memory required

    // Newton's method

    // release memory
    system("pause");
    free(incNodes);
    free(conNodes);
    free(incLoops);
    free(conLoops);
    system("pause");

    // print result (flowrate on every pipe)

    return 0;
}
