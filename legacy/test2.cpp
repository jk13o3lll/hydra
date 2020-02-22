#include <stdio.h>
#include <stdlib.h>
#include <queue>
#define RC2I(r,c,cols) ((r)*(cols)+(c))

struct Edge {
    int a, b;
    double r;
};
struct SourceNode {
    int id;
    double x;
};
struct Adjacency {
    int id, state;
};

#define bcType 0
#define n 2.0
#define nN 4
#define nN0 2
#define nE 5
#define nL 2
#define nNeq 3
#define nLeq 2
const Edge edgeList[nE] = { 
    {0, 1, 5.}, {0, 2, 1.}, {1, 2, 1.}, {1, 3, 1.}, {2, 3, 4.}  
};
const SourceNode srcList[nN0] = {
    {0, 10.}, {3, -10.}
};

int main(int argc, char *argv[]){
    Adjacency adj[nN * nN];     // adjacency matrix for pipes
    double incNodes[nNeq * nE]; // incidence matrix for node eqs, 
    double conNodes[nNeq];      // constant at the same side with unknowns
    double incLoops[nLeq * nE]; // incidence matrix for loop eqs
    double conLoops[nLeq];      // constant at the same side with unknowns
    double inlets[nN];                  // inlet flow on every nodes
    bool hasNeq[nN];                    // mark nodes to construct node eq
    int root, parent[nN], depth[nN];    // spanning tree info
    bool visited[nN];                   // record visited nodes for BFS
    bool visitedEdges[nE];              // record visited edges while BFS
    std::queue<int> buffer;             // BFS buffer
    int i, j, k, ii, jj, kk, tmp;

    printf("nE = %d, nN = %d, nN0 = %d, nL = %d\n", nE, nN, nN0, nL);
    printf("nNeq = %d, nLeq = %d\n", nNeq, nLeq);

    // construct adjacency matrix
    for(i = 0; i < nN; ++i)
        for(j = i; j < nN; ++j)
            adj[RC2I(i, j, nN)].state = adj[RC2I(j, i, nN)].state = 0;
    for(i = 0; i < nE; ++i){
        ii = RC2I(edgeList[i].a, edgeList[i].b, nN);
        jj = RC2I(edgeList[i].b, edgeList[i].a, nN);
        adj[jj].id = (adj[ii].id = i);
        adj[jj].state = -(adj[ii].state = 1);
    }

    // construct incidence matrix (node equations)
    for(i = 0; i < nN; ++i) hasNeq[i] = true, inlets[i] = 0.;
    for(i = 0; i < nN0; ++i) inlets[srcList[i].id] = srcList[i].x;

    if(bcType != 0){ // bc is presure
        for(i = 0; i < nN0; ++i) hasNeq[srcList[i].id] = false; // exclude all source nodes
        for(tmp = 0, i = 0; i < nN; ++i)
            if(hasNeq[i]){
                for(j = 0; j < nN; ++j){
                    ii = RC2I(i, j, nN);
                    if(adj[ii].state != 0)
                        incNodes[RC2I(tmp, adj[ii].id, nE)] = i == edgeList[adj[ii].id].b? 1. : -1.; // inlet as positive
                }
                conNodes[tmp] = 0.;
                ++tmp;
            }
    }
    else{   // bcType == 0, bc is flowrate
        hasNeq[srcList[nN0 - 1].id] = false; // only remove one
        for(tmp = 0, i = 0; i < nN; ++i)
            if(hasNeq[i]){
                for(j = 0; j < nN; ++j){
                    ii = RC2I(i, j, nN);
                    if(adj[ii].state != 0)
                        incNodes[RC2I(tmp, adj[ii].id, nE)] = i == edgeList[adj[ii].id].b? 1. : -1.; // inlet as positive
                }
                conNodes[tmp] = inlets[i];
                ++tmp;
            }
    }

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
            // conLoops[tmp] = srcList[i].x - srcList[0].x;
            conLoops[tmp] = srcList[0].x - srcList[i].x; // dir is 0->jj0
            ++tmp;
        }
    }

    // print eqns
    for(i = 0; i < nNeq; ++i){
        printf("Node eq %d:", i);
        for(j = 0; j < nE; ++j)
            printf(" %lf", incNodes[RC2I(i, j, nE)]);
        printf(", const = %lf\n", conNodes[i]);
    }
    for(i = 0; i < nLeq; ++i){
        printf("Loop eq %d:", i);
        for(j = 0; j < nE; ++j)
            printf(" %lf", incLoops[RC2I(i, j, nE)]);
        printf(", const = %lf\n", conLoops[i]);
    }

    // Newton's method



    // print result (flowrate on every pipe)

    return 0;
}
