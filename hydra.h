#include "solver.h"
#include <queue>

extern "C" {
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
}

void loadData(const char *filename, Edge *&edgeList, Source *&srcList, int &bcType, double &n, int &nN, int &nN0, int &nE, int &nL, int &nNeq, int &nLeq){
    FILE *fp;
    int i, j, k, a, b;
    double r;

    fp = fopen(filename, "r");
    // load parameters
    if(fp == NULL){ fprintf(stderr, "Failed to load data\n"); return; }
    fscanf(fp, "%d %d %d %d %lf", &nN, &nE, &nN0, &bcType, &n);
    if(nN <= 0 || nE <= 0 || nN <= 1 ||
        (bcType != 0 && bcType != 1) ||
        (n < 1.5 || n > 2.1)){
        fprintf(stderr, "Wrong parameters in data\n");
        return;
    }
    // allocate
    edgeList = (Edge*) malloc(nE * sizeof(Edge));
    srcList = (Source*) malloc(nN0 * sizeof(Source));
    // load pipe data
    for(i = 0; i < nE; ++i)
        fscanf(fp, "%d %d %lf", &edgeList[i].a, &edgeList[i].b, &edgeList[i].r);
    for(i = 0; i < nN0; ++i)
        fscanf(fp, "%d %lf", &srcList[i].node, &srcList[i].bc);
    fclose(fp);
    // retrieve other dimensions
    nL = nE - nN + 1;
    if(bcType == 0)         nNeq = nN - 1,   nLeq = nE - nNeq;
    else if(bcType == 1)    nNeq = nN - nN0, nLeq = nL + nN0 - 1;
}

// void toEqivalentNetwork(){}
// void toOriginalNetwork(){}

void getEquations(Edge *edgeList, Source *srcList, int bcType, double n, int nN, int nN0, int nE, int nL, int nNeq, int nLeq, double *&incNode, double *&conNode, double *&incLoop, double *&conLoop){
    int i, j, k, ii, jj, kk, tmp;
    Adjacency *adj;
    // for node eq
    bool *hasNeq;  // whether that node has node eq (sometimes src nodes are removed)
    double *inlet; // inlet flow for every node (from srcList)
    // for loop eq
    int root;       // root of spanning tree
    bool *visited;  // record visited nodes while building spanning tree
    bool *visitedE; // record vistied edge (for using non-visited edges later)
    int *parent;    // parent of every node in spanning tree
    int *detph;     // depth of every node in spanning tree
    std::queue<int> buffer; // buffer for BFS

    // convert original network to equivalent one
    // toEqivalentNetwork(); // ex. merge parallel edge, merge sources, ...

    // construct adjacency matrix
    adj = (Adjacency*) malloc(nN * nN * sizeof(Adjacency));
    for(i = 0; i < nN; ++i)
        for(j = i; j < nN; ++j)
            adj[i*nN+j].dir = adj[j*nN+i].dir = 0;
    for(i = 0; i < nE; ++i){
        ii = edgeList[i].a * nN + edgeList[i].b;
        jj = edgeList[i].b * nN + edgeList[i].a;
        adj[jj].edge = (adj[ii].edge = i);
        adj[jj].dir = -(adj[ii].dir = 1);
    }
    // // test
    // for(i = 0; i < nN; ++i){
    //     for(j = 0; j < nN; ++j)
    //         printf("%2d ", adj[i*nN+j].dir);
    //     putchar('\n');
    // }

    // construct incidence and constants matrix for node eq
    incNode = (double*) malloc(nNeq * nE * sizeof(double));
    conNode = (double*) malloc(nNeq * sizeof(double));
    memset(conNode, 0, nNeq * sizeof(double));
    memset(incNode, 0, nNeq * nE * sizeof(double));
    if(bcType == 0){ // bc is flowrate
        inlet = (double*) malloc(nN * sizeof(double));
        memset(inlet, 0, nN * sizeof(double));
        for(i = 0; i < nN0; ++i) inlet[srcList[i].node] = srcList[i].bc;
        for(tmp = 0, i = 0; i < nN; ++i)    // one node eq is dependent,
            if(i != srcList[nN0 - 1].node){ // exclude one arbitrary src node
                for(j = 0; j < nN; ++j)
                    if(adj[(ii = i*nN+j)].dir != 0)
                        incNode[tmp*nE+adj[ii].edge] =
                            i == edgeList[adj[ii].edge].b? 1.0 : -1.0;
                conNode[tmp++] = inlet[i];
            }
        free(inlet);
    }
    else if(bcType == 1){ // bc is pressure
        hasNeq = (bool*) malloc(nN * sizeof(bool));
        memset(hasNeq, 255, nN * sizeof(bool));
        for(i = 0; i < nN0; ++i) hasNeq[srcList[i].node] = false; // exclude all src
        for(tmp = 0, i = 0; i < nN; ++i)
            if(hasNeq[i]){
                for(j = 0; j < nN; ++j)
                    if(adj[(ii = i*nN+j)].dir != 0)
                        incNode[tmp*nE+adj[ii].edge] = 
                            i == edgeList[adj[ii].edge].b? 1.0 : -1.0;
                conNode[tmp++] = 0.0;
            }
        free(hasNeq);
    }
    // // test
    // for(i = 0; i < nNeq; ++i){
    //     printf("Node eq %d:", i);
    //     for(j = 0; j < nE; ++j)
    //         printf(" %.2lf", incNode[i*nE+j]);
    //     printf(", c = %.2lf\n", conNode[i]);
    // }

    // use BFS and get spanning tree (for finding independent loops)
    

    // construct incidence and constants matrix for loop eq
    incLoop = (double*) malloc(nLeq * nE * sizeof(double));
    conLoop = (double*) malloc(nLeq * sizeof(double));


    // release
    free(adj);
}

void computeF(){}
void computeJ(){}
void solve(){
    // solve by newton's method
        // has multile attempts

        // set x to initial x (rand or zero)
        // set dx = inf

        // for(nIter = 0; norm(dx) > tol && nIter < nMaxIter; ++nIter){
        //     // compute F
        //     // compute J
        //     // bicg(F, dx, J)
        //     // x += dx
        // }

        // return nIter < nMaxIter? true : false;






    // toOriginalNetwork(); // convert back to original network
}