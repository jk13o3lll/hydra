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
    int *depth;     // depth of every node in spanning tree
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
    memset(incNode, 0, nNeq * nE * sizeof(double));
    memset(conNode, 0, nNeq * sizeof(double));
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
    for(ii = 0, i = 0; i < nN; ++i){ // use node with max degree as root
        for(tmp = 0, j = 0; j < nN; ++j)
            if(adj[i*nN+j].dir != 0)
                ++tmp;
        if(tmp > ii) ii = tmp, root = i;
    }
    // printf("root is node %d (deg = %d)\n", root, ii);
    visited = (bool*) malloc(nN * sizeof(bool));
    visitedE = (bool*) malloc(nE * sizeof(bool));
    parent = (int*) malloc(nN * sizeof(int));
    depth = (int*) malloc(nN * sizeof(int));
    memset(visited, 0, nN * sizeof(bool));
    memset(visitedE, 0, nE * sizeof(bool));
    parent[root] = -1, depth[root] = 0, visited[root] = true, buffer.push(root);
    while(!buffer.empty()){
        i = buffer.front(), buffer.pop();
        for(j = 0; j < nN; ++j)
            if(adj[(ii = i*nN+j)].dir != 0 && !visited[j]){
                parent[j] = i, depth[j] = depth[i] + 1, visited[j] = true, buffer.push(j);
                visitedE[adj[ii].edge] = true;
            }
    }
    free(visited);
    // test
    for(i = 0; i < nN; ++i)
        printf("parent of node %d is node %d\n", i, parent[i]);
    for(tmp = 0, i = 0; i < nE; ++i)
        if(!visitedE[i])
            ++tmp;
    if(bcType == 0)
        printf("%d edges are not in BFS spanning tree (nLeq = %d)\n", tmp, nLeq);
    else if(bcType == 1)
        printf("%d edges are not in BFS spanning tree (nLeq - %d = %d)\n", tmp, nN0-1, nLeq-nN0+1);

    // construct incidence and constants matrix for loop eq
    incLoop = (double*) malloc(nLeq * nE * sizeof(double));
    conLoop = (double*) malloc(nLeq * sizeof(double));
    memset(incLoop, 0, nLeq * nE * sizeof(double));
    memset(conLoop, 0, nLeq * sizeof(double));
    for(tmp = 0, i = 0; i < nE; ++i)
        if(!visitedE[i]){ // start from non-visited edge
            for(ii = edgeList[i].a , jj = edgeList[i].b; ii != jj; ) // find lowest common ancient
                if(depth[ii] > depth[jj]){ // for ii, parent -> child is same as loop dir
                    kk = parent[ii] * nN + ii;
                    k = adj[kk].edge;
                    // whether direction of loop is same as defined flow direction
                    incLoop[tmp*nE+k] = adj[kk].dir == 1? edgeList[k].r : -edgeList[k].r;
                    ii = parent[ii];
                }
                else{ // depth[ii] <= depth[jj], for jj, parent -> child is opposite to loop dir
                    kk = jj * nN + parent[jj];
                    k = adj[kk].edge;
                    // whether direction of loop is same as defined flow direction
                    incLoop[tmp*nE+k] = adj[kk].dir == 1? edgeList[k].r : -edgeList[k].r;
                    jj = parent[ii];
                }
            // use edgeList[i].a -> edgeList[i].b as loop dir (loop dir = edge dir)
            incLoop[tmp*nE+i] = edgeList[i].r;
            conLoop[tmp++] = 0.0;
        }
    if(bcType == 1){ // additional loops of passing src (if bc is pressure)
        for(i = 1; i < nN0; ++i){ // connect the first src node to others
            for(ii = srcList[0].node, jj = srcList[i].node; ii != jj; )
                if(depth[ii] > depth[jj]){
                    kk = parent[ii] * nN + ii;
                    k = adj[kk].edge;
                    // whether direction of loop is same as defined flow direction
                    incLoop[tmp*nE+k] = adj[kk].dir == 1? edgeList[k].r : -edgeList[k].r;
                    ii = parent[ii];
                }
                else{ // depth[ii] <= depth[jj]
                    kk = jj * nN + parent[jj];
                    k = adj[kk].edge;
                    // whether direction of loop is same as defined flow direction
                    incLoop[tmp*nE+k] = adj[kk].dir == 1? edgeList[k].r : -edgeList[k].r;
                    jj = parent[ii];
                }
            // constant pressure change between source nodes
            // use srcList[0].node -> srcList[i].node as loop dir (rQ is presure drop)
            conLoop[tmp++] = srcList[0].bc - srcList[i].bc; // watch out dir
        }
    }
    free(visitedE); free(parent); free(depth); free(adj);
    // // test
    // for(i = 0; i < nLeq; ++i){
    //     printf("Loop eq %d:", i);
    //     for(j = 0; j < nE; ++j)
    //         printf(" %.2lf", incLoop[i*nE+j]);
    //     printf(", c = %.2lf\n", conLoop[i]);
    // }
}

// A(x) + c = 0
// R is vector, R = -F = -(A(x) + c)
void computeR(int nE, int nLeq, int nNeq, double *incLoop, double *conLoop, double *incNode, double *conNode, double *x, double n, double *tmp, double *res){
    int i, j, offset = nLeq;
    #pragma omp parallel for
    for(i = 0; i < nE; ++i) // tmp for store pow(x)
        tmp[i] = pow(x[i], n);
    for(i = 0; i < nLeq; ++i){
        res[i] = -conLoop[i];
        for(j = 0; j < nE; ++j)
            res[i] -= incLoop[i*nE+j] * tmp[j];
    }
    for(i = 0; i < nNeq; ++i){
        res[offset + i] = -conNode[i];
        for(j = 0; j < nE; ++j)
            res[offset + i] -= incNode[i*nE+j] * x[j];
    }
}
// J is matrix, J = Jacobian(F)
void computeJ(int nE, int nLeq, int nNeq, double *incLoop, double *incNode, double *x, double n, double *tmp, double *res){
    int i, j, offset = nLeq * nE;
    #pragma omp parallel for
    for(i = 0; i < nE; ++i)
        tmp[i] = pow(x[i], n - 1.0);
    for(i = 0; i < nLeq; ++i)
        for(j = 0; j < nE; ++j)
            res[i*nE+j] = incLoop[i*nE+j] * n * tmp[i];
    for(i = 0; i < nNeq; ++i)
        for(j = 0; j < nE; ++j)
            res[offset + i*nE+j] = incNode[i*nE+j];
}
void solve(int nE, int nLeq, int nNeq, double *incLoop, double *conLoop, double *incNode, double *conNode, double *&x, double n, int maxiter = 100000, int maxattempts = 10){
    int i, j, k;
    double *R, *J, *bufferx, *dx;
    // allocate
    R = (double*) malloc(nE * sizeof(double));
    J = (double*) malloc((nLeq + nNeq) * nE * sizeof(double));
    bufferx = (double*) malloc(nE * sizeof(double));
    dx = (double*) malloc(nE * sizeof(double));
    x = (double*) malloc(nE * sizeof(double));
    // solve by newton's method
    for(i = 0; i < 100; ++i){ // has multiple attempts
        printf("Attempt %d\n", i + 1);
        // zerox(nE, x);
        randx(nE, x, 20.0, -10.0);
        for(j = 0; j < maxiter; ++j){
            // get residual and jacobian
            computeR(nE, nLeq, nNeq, incLoop, conLoop, incNode, conNode, x, n, bufferx, R);
            computeJ(nE, nLeq, nNeq, incLoop, incNode, x, n, bufferx, J);
            // get dx and update
            // bicgstab(nE, J, R, dx);
            gaussian(nE, J, R, dx);
            xplussy(nE, x, 0.5, dx, x); // half step size
            // check
            if((k = allzero(nE, R, 1e-3)) > 0) break;
            // debug
            // printA("J = \n", nE, J);
            printx("R = ", nE, R);
            // printx("dx = ", nE, dx);
            // printx("x = ", nE, x);
            // getchar();
        }
        if(j == maxiter)
            puts("Caonnot converge");
        else if(k == 2)
            puts("Detect NaN");
        else{
            printf("Converge at iteration %d\n", j + 1);
            break;
        }
    }
    free(R); free(J); free(bufferx); free(dx);

    // toOriginalNetwork(); // convert back to original network
}