#include <stdio.h>
#include <stdlib.h>
#include <queue>

#define RC2I(r,c,cols) ((r)*(cols)+(c))

struct Edge {
    int a, b;
    double r;
};
struct Source {
    int node;
    double bc;
};
struct Adjacency {
    int edge;
    int dir;
};

bool loadData(const char *filename, Edge *&edgeList, Source *&srcList, int &bcType, double &n, int &nN, int &nN0, int &nE, int &nL, int &nNeq, int &nLeq){
    int i;
    FILE *fp;
    int a, b;
    double r;

    // load data
    fp = fopen(filename, "r");
    if(fp == NULL) return false;
    fscanf(fp, "%d%d%d%d%lf", &nN, &nE, &nN0, &bcType, &n);
    if(nN <= 0 || nE <= 0 || nN <= 1 || (bcType != 0 && bcType != 1) || (n < 1.5 || n > 2.1))
        return false;
    edgeList = (Edge*)malloc(nE * sizeof(Edge));
    srcList = (Source*)malloc(nN0 * sizeof(Source));
    for(i = 0; i < nE; ++i)
        fscanf(fp, "%d%d%lf", &edgeList[i].a, &edgeList[i].b, &edgeList[i].r);
    for(i = 0; i < nN0; ++i)
        fscanf(fp, "%d%lf", &srcList[i].node, &srcList[i].bc);
    fclose(fp);

    // get other dimensions
    nL = nE - nN + 1;
    if(bcType == 0){      nNeq = nN - 1;    nLeq = nE - nNeq; }
    else if(bcType == 1){ nNeq = nN - nN0;  nLeq = nL + nN0 - 1; }
    return true;
}

// bool toEquivalentNetwork(vector<int> original[nN][nN], int eqivalent[nN][nN]){
    // generate new edge list based on Q'
    // (necessary) 1. flowrate sources to one node (just merge)
    // (necessary) 2. merge multiple edges between two node (compute equivalent r)
    // (optional) 3. form adj, check degree, if deg==2, travel two endpoints until stop and replace by one edge
    // generate k-list nad o-list
//     return true;
// }

// bool fromEquivalentNetwork(vector<int> original[nN][nN], int eqivalent[nN][nN], double x[nE]){
    // normal Q = Q'
    // parallel Q = (r_eq / r)^(1/n) Q' = kQ'
    // serial Q1 = Q'; Q2 = Q'

    // k[] = [+1, -1, k0, 1, k1, ...] include dir also
    // o[] = [0, 1, 2, 2, 3, 4, ...]
    // Q[i] = k[i] * Q'[o[i]]
//     return true;
// }

// this function is not suitable for parallelize
bool getEquations(const Edge *edgeList, const Source *srcList, const int bcType, const double n, const int nN, const int nN0, const int nE, const  int nL, const int nNeq, const int nLeq, double *&incNode, double *&conNode, double *&incLoop, double *&conLoop){
    Adjacency *adj;
    // for node eq
    bool *hasNeq;
    double *inlet;
    // for loop eq
    int root;
    bool *visited, *visitedE;
    int *parent, *depth;
    std::queue<int> buffer;
    // tmp
    int i, j, k, ii, jj, kk, tmp;

    // construct adjacency matrix
    adj = (Adjacency*)malloc(nN * nN * sizeof(Adjacency));
    for(i = 0; i < nN; ++i)
        for(j = i; j < nN; ++j)
            adj[RC2I(i, j, nN)].dir = adj[RC2I(j, i, nN)].dir = 0;
    for(i = 0; i < nE; ++i){
        ii = RC2I(edgeList[i].a, edgeList[i].b, nN);
        jj = RC2I(edgeList[i].b, edgeList[i].a, nN);
        adj[jj].edge = (adj[ii].edge = i);
        adj[jj].dir = -(adj[ii].dir = 1);
    }

    // construct incidence and constants matrix for node eq
    incNode = (double*)malloc(nNeq * nE * sizeof(double));
    conNode = (double*)malloc(nNeq * sizeof(double));
    for(i = 0; i < nNeq; ++i){
        conNode[i] = 0.;
        for(j = 0; j < nE; ++j)
            incNode[RC2I(i, j, nE)] = 0.;
    }
    hasNeq = (bool*)malloc(nN * sizeof(bool));
    for(i = 0; i < nN; ++i) hasNeq[i] = true;
    if(bcType == 0){
        inlet = (double*)malloc(nN * sizeof(double));
        for(i = 0; i < nN; ++i)  inlet[i] = 0.;
        for(i = 0; i < nN0; ++i) inlet[srcList[i].node] = srcList[i].bc;
        hasNeq[srcList[nN0 - 1].node] = false;
        for(tmp = 0, i = 0; i < nN; ++i)
            if(hasNeq[i]){
                for(j = 0; j < nN; ++j){
                    ii = RC2I(i, j, nN);
                    if(adj[ii].dir != 0)
                        incNode[RC2I(tmp, adj[ii].edge, nE)] = i == edgeList[adj[ii].edge].b? 1. : -1.;
                }
                conNode[tmp++] = inlet[i];
            }
        free(inlet);
    }
    else if(bcType == 1){
        for(i = 0; i < nN0; ++i) hasNeq[srcList[i].node] = false;
        for(tmp = 0, i = 0; i < nN; ++i)
            if(hasNeq[i]){
                for(j = 0; j < nN; ++j){
                    ii = RC2I(i, j, nN);
                    if(adj[ii].dir !=0)
                        incNode[RC2I(tmp, adj[ii].edge, nE)] = i == edgeList[adj[ii].edge].b? 1. : -1.;
                }
                conNode[tmp++] = 0.;
            }
    }
    free(hasNeq);

    // construct incidence and constants marrix for loop eq
    // use max degree node as root
    for(ii = 0, i = 0; i < nN; ++i){
        for(tmp = 0, j = 0; j < nN; ++j)
            if(adj[RC2I(i, j, nN)].dir != 0)
                ++tmp;
        if(tmp > ii){ ii = tmp; root = i; }
    }
    // BFS and get spanning tree
    visited = (bool*)malloc(nN * sizeof(bool));
    visitedE = (bool*)malloc(nE * sizeof(bool));
    parent = (int*)malloc(nN * sizeof(int));
    depth = (int*)malloc(nN * sizeof(int));
    for(i = 0; i < nN; ++i) visited[i] = false;
    for(i = 0; i < nE; ++i) visitedE[i] = false;
    parent[root] = -1; depth[root] = 0; visited[root] = true; buffer.push(root);
    while(!buffer.empty()){
        i = buffer.front(); buffer.pop();
        for(j = 0; j < nN; ++j){
            ii = RC2I(i, j, nN);
            if(adj[ii].dir != 0 && !visited[j]){
                parent[j] = i; depth[j] = depth[i] + 1; visited[j] = true; buffer.push(j);
                visitedE[adj[ii].edge] = true;
            }
        }
    }
    free(visited);
    // for(i = 0; i < nN; ++i)
    //     printf("parent of node %d is node %d\n", i, parent[i]);
    // system("pause");
    // for(tmp = 0, i = 0; i < nE; ++i)
    //     if(!visitedE[i])
    //         ++tmp;
    // printf("%d edges are not in BFS spanning tree\n", tmp);

    // form loop by connecting edges not in tree
    incLoop = (double*)malloc(nLeq * nE * sizeof(double));
    conLoop = (double*)malloc(nLeq * sizeof(double));
    for(i = 0; i < nLeq; ++i){
        conLoop[i] = 0.;
        for(j = 0; j < nE; ++j)
            incLoop[RC2I(i, j, nE)] = 0.;
    }
    for(tmp = 0, i = 0; i < nE; ++i)
        if(!visitedE[i]){
            for(ii = edgeList[i].a, jj = edgeList[i].b; ii != jj; )
                if(depth[ii] > depth[jj]){
                    kk = RC2I(parent[ii], ii, nN); // assume a->b as loop dir, dir is parent to child
                    k = adj[kk].edge;
                    incLoop[RC2I(tmp, k, nE)] = adj[kk].dir == 1? edgeList[k].r : -edgeList[k].r;
                    ii = parent[ii];
                }
                else{ // depth[ii] <= depth[jj]
                    kk = RC2I(jj, parent[jj], nN);
                    k = adj[kk].edge;
                    incLoop[RC2I(tmp, k, nE)] = adj[kk].dir == 1? edgeList[k].r : -edgeList[k].r;
                    jj = parent[jj];
                }
            incLoop[RC2I(tmp, i, nE)] = edgeList[i].r;
            conLoop[tmp++] = 0.;
        }
    free(visitedE);
    // for source loop by coonecting src node
    if(bcType == 1){
        for(i = 1; i < nN0; ++i){
            for(ii = srcList[0].node, jj = srcList[i].node; ii != jj; )
                if(depth[ii] > depth[jj]){
                    kk = RC2I(parent[ii], ii, nN);
                    k = adj[kk].edge;
                    incLoop[RC2I(tmp, k, nE)] = adj[kk].dir == 1? edgeList[k].r : -edgeList[k].r;
                    ii = parent[ii];
                }
                else{ // depth[ii] <= depth[jj]
                    kk = RC2I(jj, parent[jj], nN);
                    k = adj[kk].edge;
                    incLoop[RC2I(tmp, k, nE)] = adj[kk].dir == 1? edgeList[k].r : -edgeList[k].r;
                    jj = parent[jj];
                }
            conLoop[tmp++] = srcList[0].bc - srcList[i].bc;
        }
    }
    free(parent);
    free(depth);
    // printf("%d loop eqns are founded\n", tmp);

    // release memory
    free(adj);
    return true;
}

bool newton(const double *incNode, const double *conNode, const double *incLoop, const double *conLoop, const double n, double *&x, const int nNeq, const int nLeq, const int nE){

    // set x to initial x (rand or zero)
    // set dx = inf

    // for(nIter = 0; norm(dx) > tol && nIter < nMaxIter; ++nIter){
    //     // compute F
    //     // compute J
    //     // bicg(F, dx, J)
    //     // x += dx
    // }

    // return nIter < nMaxIter? true : false;
}

int main(int argc, char *argv[]){
    Edge *edgeList;
    Source *srcList;
    int bcType;
    double n;
    int nN, nN0, nE, nL, nNeq, nLeq;
    double *incNode, *conNode, *incLoop, *conLoop; // constants at the same side with unknowns
    double *x;
    bool ret;
    int i, j;

    // load data
    if(argc != 2){ puts("Wrong input arguments"); return 1; }
    ret = loadData(argv[1], edgeList, srcList, bcType, n, nN, nN0, nE, nL, nNeq, nLeq);
    if(!ret){ puts("Fail to load data"); return 1; }
    printf("nE = %d, nN = %d, nN0 = %d, nL = %d\n", nE, nN, nN0, nL);
    printf("nNeq = %d, nLeq = %d\n", nNeq, nLeq);
    printf("n = %lf\n", n);
    printf("edgeList =\n");
    for(i = 0; i < nE; ++i)
        printf("%d %d %lf\n", edgeList[i].a, edgeList[i].b, edgeList[i].r);
    printf("srcList =\n");
    for(i = 0; i < nN0; ++i)
        printf("%d %lf\n", srcList[i].node, srcList[i].bc);

    // transform into equivalent network
    // toEquivalentNetwork();

    // get equations
    ret = getEquations(edgeList, srcList, bcType, n, nN, nN0, nE, nL, nNeq, nLeq, incNode, conNode, incLoop, conLoop);
    free(edgeList);
    free(srcList);
    if(!ret){ puts("Fail to get equations"); return 1; }

    // show equations
    for(i = 0; i < nNeq; ++i){
        printf("Node eq %d:", i);
        for(j = 0; j < nE; ++j)
            printf(" %.2lf", incNode[RC2I(i, j, nE)]);
        printf(", con = %.2lf\n", conNode[i]);
    }
    for(i = 0; i < nLeq; ++i){
        printf("Loop eq %d:", i);
        for(j = 0; j < nE; ++j)
            printf(" %.2lf", incLoop[RC2I(i, j, nE)]);
        printf(", con = %.2lf\n", conLoop[i]);
    }

    // Newton's method
    // ret = newton();
    free(incNode);
    free(conNode);
    free(incLoop);
    free(conLoop);
    if(!ret){ puts("Exceed maximum iteration"); return 1; }

    // tranform from equivalent network to original network
    // fromEquivalentNetwork();

    // save result
    // fp = fopen()
    // for  fprintf(fp, x)
    free(x);

    return 0;
}