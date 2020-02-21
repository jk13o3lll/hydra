#define NETWORK1
// #define NETWORK2
// #define NETWORK3
// #define NETWORK4

#include <stdio.h>
#include <stdlib.h>
#include <queue>

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

// basic data type
struct Edge {
    int a, b;
    double r;   // r is constant related to pipe diameter, pipe length, and friction factors
};
struct SourceNode {
    int id;
    double x;   // flow or pressure
};
struct Adjacency {
    int id;     // id of flow on that edge
    int state;  // -1 for direction opposite to defined flow direction,
                // 0 for no connnection,
                // +1 for direction indentical to defined flow direction.
};


#ifdef NETWORK1
// network 1 (boundary condition by flowrate)
const int bcType = 0; // boundary condition is flowrate to node
const double n = 2.;    // constant in pressure loss eqn, 1.5 < n < 2.0
const Edge edgeList[] = { 
    {0, 1, 5.}, // a, b, r (from a to b, with pipe radius = r)
    {0, 2, 1.}, // all node should be contain in single graph
    {1, 2, 1.}, // (should not have multple separate graphs)
    {1, 3, 1.}, // id of node should start from zero and continuous
    {2, 3, 4.}  // use original direction as assumed direction to
};              // generate adjacency matrix (adj store +1 or -1) and incidence matrix
const SourceNode srcList[] = {  // direction for node eqation is 
    {0, 10.},                   // flow into node is positive
    {3, -10.}                   // flow out from node is negative
};                              // p.s. one node should be dependent due to mass conservation
#endif // NETWORK1

#ifdef NETWORK2
// network 2 (boundary condition by flowrate)
const int bcType = 0; // boundary condition is flowrate to node
const double n = 2.;
const Edge edgeList[] = { 
    {0, 1, 1.},
    {1, 6, 5.},
    {1, 2, 1.},
    {0, 3, 5.},
    {3, 4, 3.},
    {2, 4, 1.},
    {2, 5, 1.},
    {5, 6, 3.},
    {3, 8, 3.},
    {4, 7, 1.},
    {5, 7, 1.},
    {7, 8, 2.},
    {8, 9, 3.},
    {6, 9, 3.},
};
const SourceNode srcList[] = {
    {0, 10.},
    {9, -10.}
};
#endif // NETWORK2

#ifdef NETWORK3
// network 3 (boundary condition by pressure)
const int bcType = 1; // boundary condition is pressure at node
const double n = 2.;
const Edge edgeList[] = {
    {0, 1, 1.},
    {1, 6, 5.},
    {1, 2, 1.},
    {0, 3, 5.},
    {3, 4, 3.},
    {2, 4, 1.},
    {2, 5, 1.},
    {5, 6, 3.},
    {3, 14, 1.},
    {4, 7, 1.},
    {5, 7, 1.},
    {7, 8, 2.},
    {8, 9, 3.},
    {6, 9, 3.},
    {8, 10, 1.},
    {14, 11, 2.},
    {13, 16, 5.},
    {14, 15, 5.},
    {10, 11, 3.},
    {11, 12, 2.},
    {9, 16, 1.},
    {12, 15, 2.},
    {10, 13, 2.},
    {13, 12, 1.},
    {16, 17, 1.},
    {15, 17, 2.},
    {18, 0, 1.}, // edge connect to source
    {9, 19, 1.}, // edge connect to source
    {17, 20, 1.} // edge connect to source
};
const SourceNode srcList[] = {  // direction for loop equation is
    {18, 10.},                  // same as loop direction is positive
    {19, 0.},                   // opposite to loop direction is negative
    {20, 0.}    // source node should be exclude when construct node eq
};
#endif // NETWORK3

#ifdef NETWORK4
// network 4 (boundary condition by pressure)
const int bcType = 1; // boundary condition is pressure at node
const double n = 2.;
const Edge edgeList[] = {
    {0, 1, 1.},
    {1, 6, 5.},
    {1, 2, 1.},
    {0, 3, 5.},
    {3, 4, 3.},
    {2, 4, 1.},
    {2, 5, 1.},
    {5, 6, 3.},
    {3, 14, 1.},
    {4, 7, 1.},
    {5, 7, 1.},
    {7, 8, 2.},
    {8, 9, 3.},
    {6, 9, 3.},
    {8, 10, 1.},
    {14, 11, 2.},
    {13, 16, 5.},
    {14, 15, 5.},
    {10, 11, 3.},
    {11, 12, 2.},
    {9, 16, 1.},
    {12, 15, 2.},
    {10, 13, 2.},
    {13, 12, 1.},
    {16, 17, 1.},
    {15, 17, 2.},
    {18, 0, 1.}, // edge connect to source
    {9, 19, 1.}, // edge connect to source
    {17, 20, 1.} // edge connect to source
};
const SourceNode srcList[] = {
    {18, 100.},
    {19, 0.},
    {20, 0.}
};
#endif // NETWORK4

// get basic dimensions
const int nE = sizeof(edgeList) / sizeof(Edge); // number of edges
const int nN0 = sizeof(srcList) / sizeof(SourceNode);
int nN = 0, nL = 0; // number of nodes (include sources), number of loops, obtained later

// To calculate matrix inverse (square matrix) by BiCG or BiCGSTAB
void inv(double *A, double *Ainv, const int nRows, const int nCols){
    return;
}


// main function
int main(int argc, char *argv[]){
    int i, j;
    Adjacency *adj;     // adjacency matrix for pipes, value is +1 or -1
                        // (+1 for originally defined direction, -1 for opposite)
    int nNeq, nLeq;     // number node eqns, number loop eqns
    double *incNodes;   // incidence matrix for node equations
    double *incLoops;   // incidence matrix for loop equations

    // ======================================================================

    // get dimensions
    for(i = 0; i < nE; ++i){
        if(edgeList[i].a > nN)  nN = edgeList[i].a;
        if(edgeList[i].b > nN)  nN = edgeList[i].b;
    }
    ++nN;
    nL = nE - nN + 1;

    // allocate memory
    adj = (Adjacency*) malloc(nN * nN * sizeof(Adjacency));
    if(bcType == 0){
        nNeq = nN - 1;          // ignore one dependent source
        nLeq = nE - nNeq;
    }
    else if(bcType == 1){
        nNeq = nN - nN0;        // ignore all source nodes
        nLeq = nL + nN0 - 1;    // include loops for source nodes
    }
    incNodes = (double*) malloc(nNeq * nE * sizeof(double));
    incLoops = (double*) malloc(nLeq * nE * sizeof(double));

    // --- test ---
    printf("nE = %d, nN = %d, nN0 = %d, nL = %d\n", nE, nN, nN0, nL); // dimensions
    printf("nNeq = %d, nLeq = %d\n", nNeq, nLeq); // num of eqs
    // construct adjacency matrix

    // construct incidence matrix (node equations)

    // construct incidence matrix (loop equations)

    // release memory not necessary in the following
    free(adj);

    // ======================================================================

    // allocate memory required

    // Newton's method

    // print result

    // release all memory

    return 0;
}

// // test spaning tree
// int i, j, a, b, tmp;
// int maxDegree;
// int adj[N][N], adj2[N][N];      // adjacency matrix (undirected) (weight is radius)
// int edge[E][3];                 // node1, node2, weight
// int root, parent[N], depth[N];  // spanning tree info
// bool visit[N];                  // for BFS
// std::queue<int> buffer;         // BFS buffer
// int loop[E][N][N];              // TODO: change to more efficient data structure later
// int edgeflow[E][2];             // assume edge flow from node a to node b, where a < b
// int A[E][N];                    // edgeflow-node incidence matrix
// int M[L][E];                    // loop-edgeflow incidence matrix
// // adj to store edge id for faster searching

// // init arrays
// // construct and copy adj
// // ...

// // find node with max degree as root
// maxDegree = -1;
// for(i = 0; i < N; ++i){
//     for(tmp = 0, j = 0; j < N; ++j)
//         if(adj[i][j] > 0)
//             ++tmp;
//     if(tmp > maxDegree)
//         maxDegree = tmp, root = i;
// }

// // bfs
// parent[root] = -1, depth[root] = 0, visit[root] = true, buffer.push(i);
// while(!buffer.empty()){
//     tmp = buffer.front(); buffer.pop();
//     for(i = 0; i < N; ++i)
//         if(adj[tmp][i] > 0 && !visit[i]){
//             parent[i] = tmp, depth[i] = depth[tmp] + 1, visit[i] = true, buffer.push(i);
//             adj2[tmp][i] = adj2[i][tmp] = 0; // edges left in adj2 will be used later
//         }
// }

// // find indepednet loops
// for(tmp = 0, i = 0; i < N; ++i)
//     for(j = i + 1; j < N; ++j)
//         if(adj2[i][j] > 0){
//             loop[tmp][i][j] = -1, loop[tmp][j][i] = 1;
//             // form the loop, watch out direciton of edge
//             for(a = i, b = j; a != b; )
//                 if(depth[a] > depth[b])
//                     loop[tmp][a][parent[a]] = 1, loop[tmp][parent[a]][a] = -1, a = parent[a];
//                 else // depth[a] <= depth[b], loop dir is opposite to a
//                     loop[tmp][b][parent[b]] = -1, loop[tmp][parent[b]][b] = 1, b = parent[b];
//             ++tmp;
//         }

// // construct equation

// // solve equations by Newton's (use bicg to compute inverse for gradient)