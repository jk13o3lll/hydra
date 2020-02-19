#include <stdio.h>
#include <stdlib.h>
#include <queue>

#define N 10    // number of nodes
#define E 100   // number of edges
#define L 0     // number of loops

int main(int argc, char *argv[]){
    // test spaning tree
    int i, j, a, b, tmp;
    int maxDegree;
    int adj[N][N], adj2[N][N];      // adjacency matrix (undirected) (weight is radius)
    int edge[E][3];                 // node1, node2, weight
    int root, parent[N], depth[N];  // spanning tree info
    bool visit[N];                  // for BFS
    std::queue<int> buffer;         // BFS buffer
    int loop[E][N][N];              // TODO: change to more efficient data structure later
    int edgeflow[E][2];             // assume edge flow from node a to node b, where a < b
    int A[E][N];                    // edgeflow-node incidence matrix
    int M[L][E];                    // loop-edgeflow incidence matrix
    // adj to store edge id for faster searching

    // init arrays
    // construct and copy adj
    // ...

    // find node with max degree as root
    maxDegree = -1;
    for(i = 0; i < N; ++i){
        for(tmp = 0, j = 0; j < N; ++j)
            if(adj[i][j] > 0)
                ++tmp;
        if(tmp > maxDegree)
            maxDegree = tmp, root = i;
    }

    // bfs
    parent[root] = -1, depth[root] = 0, visit[root] = true, buffer.push(i);
    while(!buffer.empty()){
        tmp = buffer.front(); buffer.pop();
        for(i = 0; i < N; ++i)
            if(adj[tmp][i] > 0 && !visit[i]){
                parent[i] = tmp, depth[i] = depth[tmp] + 1, visit[i] = true, buffer.push(i);
                adj2[tmp][i] = adj2[i][tmp] = 0; // edges left in adj2 will be used later
            }
    }

    // find indepednet loops
    for(tmp = 0, i = 0; i < N; ++i)
        for(j = i + 1; j < N; ++j)
            if(adj2[i][j] > 0){
                loop[tmp][i][j] = -1, loop[tmp][j][i] = 1;
                // form the loop, watch out direciton of edge
                for(a = i, b = j; a != b; )
                    if(depth[a] > depth[b])
                        loop[tmp][a][parent[a]] = 1, loop[tmp][parent[a]][a] = -1, a = parent[a];
                    else // depth[a] <= depth[b], loop dir is opposite to a
                        loop[tmp][b][parent[b]] = -1, loop[tmp][parent[b]][b] = 1, b = parent[b];
                ++tmp;
            }
    
    // construct equation

    // solve equations by Newton's (use bicg to compute inverse for gradient)
            

    return 0;
}