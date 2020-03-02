// TODO: add bicg also
// TODO: check and optimize

#include "solver2.h"
#include <queue>
#include <stack>

extern "C" {
    struct Edge {
        int a, b;
        long double r;
    };
    struct Source {
        int node;
        long double bc;
    };
    struct Adjacency {
        int edge;
        int dir;
    };
}

void loadData(const char *filename, Edge *&edgeList, Source *&srcList, int &bcType, long double &n, int &nN, int &nN0, int &nE, int &nL, int &nNeq, int &nLeq){
    FILE *fp;
    int i, j, k, a, b;
    long double r;

    fp = fopen(filename, "r");
    // load parameters
    if(fp == NULL){ fprintf(stderr, "Failed to load data\n"); return; }
    fscanf(fp, "%d %d %d %d %Lf", &nN, &nE, &nN0, &bcType, &n);
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
        fscanf(fp, "%d %d %Lf", &edgeList[i].a, &edgeList[i].b, &edgeList[i].r);
    for(i = 0; i < nN0; ++i)
        fscanf(fp, "%d %Lf", &srcList[i].node, &srcList[i].bc);
    fclose(fp);
    // retrieve other dimensions
    nL = nE - nN + 1;
    if(bcType == 0)         nNeq = nN - 1,   nLeq = nE - nNeq;
    else if(bcType == 1)    nNeq = nN - nN0, nLeq = nL + nN0 - 1;
}

// void toEqivalentNetwork(){}
// void toOriginalNetwork(){}

void getEquations(Edge *edgeList, Source *srcList, int bcType, long double n, int nN, int nN0, int nE, int nL, int nNeq, int nLeq, long double *&incNode, long double *&conNode, long double *&incLoop, long double *&conLoop){
    int i, j, k, ii, jj, kk, tmp;
    Adjacency *adj;
    // for node eq
    bool *hasNeq;  // whether that node has node eq (sometimes src nodes are removed)
    long double *inlet; // inlet flow for every node (from srcList)
    // for loop eq
    int root;       // root of spanning tree
    bool *visited;  // record visited nodes while building spanning tree
    bool *visitedE; // record vistied edge (for using non-visited edges later)
    int *parent;    // parent of every node in spanning tree
    int *depth;     // depth of every node in spanning tree
    std::queue<int> buffer; // buffer for BFS
    // std::stack<int> buffer; // buffer for DFS

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
    incNode = (long double*) malloc(nNeq * nE * sizeof(long double));
    conNode = (long double*) malloc(nNeq * sizeof(long double));
    memset(incNode, 0, nNeq * nE * sizeof(long double));
    memset(conNode, 0, nNeq * sizeof(long double));
    if(bcType == 0){ // bc is flowrate
        inlet = (long double*) malloc(nN * sizeof(long double));
        memset(inlet, 0, nN * sizeof(long double));
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
    //         printf(" %.2Lf", incNode[i*nE+j]);
    //     printf(", c = %.2Lf\n", conNode[i]);
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
        i = buffer.front(), buffer.pop(); // BFS
        // i = buffer.top(), buffer.pop(); // DFS
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
    incLoop = (long double*) malloc(nLeq * nE * sizeof(long double));
    conLoop = (long double*) malloc(nLeq * sizeof(long double));
    memset(incLoop, 0, nLeq * nE * sizeof(long double));
    memset(conLoop, 0, nLeq * sizeof(long double));
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
                    jj = parent[jj];
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
                    jj = parent[jj];
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
    //         printf(" %.2Lf", incLoop[i*nE+j]);
    //     printf(", c = %.2Lf\n", conLoop[i]);
    // }
}

// A(x) + c = 0
// R is vector, R = -F = -(A(x) + c)
void computeR(int nE, int nLeq, int nNeq, long double *incLoop, long double *conLoop, long double *incNode, long double *conNode, long double *x, long double n, long double *tmp, long double *res){
    int i, j, offset = nLeq;
    // #pragma omp parallel for
    for(i = 0; i < nE; ++i) // tmp for store pow(x)
        tmp[i] = x[i] * pow(fabs(x[i]), n-1.0);
        // tmp[i] = x[i] * fabs(x[i]); // n == 2
    for(i = 0; i < nLeq; ++i){
        res[i] = conLoop[i];
        for(j = 0; j < nE; ++j)
            res[i] += incLoop[i*nE+j] * tmp[j];
    }
    for(i = 0; i < nNeq; ++i){
        res[offset + i] = conNode[i];
        for(j = 0; j < nE; ++j)
            res[offset + i] += incNode[i*nE+j] * x[j];
    }
}
// J is matrix, J = Jacobian(F)
void computeJ(int nE, int nLeq, int nNeq, long double *incLoop, long double *incNode, long double *x, long double n, long double *tmp, long double *res){
    int i, j, offset = nLeq * nE;
    // #pragma omp parallel for
    for(i = 0; i < nE; ++i)
        tmp[i] = (n-1.0)*pow(fabs(x[i]), n-3.0) + pow(fabs(x[i]), n-1.0);
        // tmp[i] = fabs(x[i]); // n == 2
    for(i = 0; i < nLeq; ++i)
        for(j = 0; j < nE; ++j)
            res[i*nE+j] = incLoop[i*nE+j] * tmp[i];
    memcpy(res+offset, incNode, nNeq * nE * sizeof(long double));
}
void solve(int nE, int nLeq, int nNeq, long double *incLoop, long double *conLoop, long double *incNode, long double *conNode, long double *&x, long double n, int maxiter = 100000, int maxattempts = 1000){
    int i, j, k;
    long double *R, *J, *bufferx, *dx;
    srand(time(NULL));
    // allocate
    R = (long double*) malloc(nE * sizeof(long double));
    J = (long double*) malloc((nLeq + nNeq) * nE * sizeof(long double));
    bufferx = (long double*) malloc(nE * sizeof(long double));
    dx = (long double*) malloc(nE * sizeof(long double));
    x = (long double*) malloc(nE * sizeof(long double));

    // solve by newton's method
    for(i = 0; i < maxattempts; ++i){ // has multiple attempts
        printf("Attempt %d\n", i + 1);
        // randx(nE, x, 1.0, -0.5);
        randx(nE, x, 2.0, -1.0);
        // randx(nE, x, 0.2, -0.1);
        // randx(nE, x, 10.0);
        // randx(nE, x, 100.0, -50.0);
        // randx(nE, x, 0.02, 0.99);
        // randx(nE, x, 0.5, 0.1);
        for(j = 0; j < maxiter; ++j){
            // get residual and jacobian
            computeR(nE, nLeq, nNeq, incLoop, conLoop, incNode, conNode, x, n, bufferx, R);
            computeJ(nE, nLeq, nNeq, incLoop, incNode, x, n, bufferx, J);
            // get dx and update
            // bicgstab(nE, J, R, dx); // cannot converge for large matrix?
            // gaussian(nE, J, R, dx);
            gaussian2(nE, J, R, dx, bufferx);
            // xminussy(nE, x, dx, x);
            // xminussy(nE, x, 0.1, dx, x); // 0.01 step size
            xminussy(nE, x, 0.01, dx, x); // 0.01 step size (for network5)
            // check
            if((k = allzero(nE, R, 1e-6)) > 0) break;
            if(allzero(nE, R, 1e8) == 0){ k = 3; break; } // too large
            // debug
            // printA("J = \n", nE, J);
            // printx("R = ", nE, R);
            // printx("dx = ", nE, dx);
            // printx("x = ", nE, x);
            // getchar();
        }
        printx("R = ", nE, R);
        printx("x = ", nE, x);
        if(j == maxiter || k == 3)
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