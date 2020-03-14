#include "hydra.h"

int all_within(int n, double *x, double lowerb, double upperb){
    int i;
    for(i = 0; i < n; ++i)
        if(isnan(x[i]) || x[i] < lowerb || x[i] > upperb)
            break;
    if(i == n)              return 1; // all are within
    else if(isnan(x[i]))    return 2; // contain NaN
    else                    return 0; // some are not within
}
void rands(int n, double *x, double a, double b){
    for(int i = 0; i < n; ++i)
        x[i] = (double)rand() / RAND_MAX * a + b;
}
void zeros(int n, double *x){
    memset(x, 0, n * sizeof(double));
}
void ones(int n, double *x){
    for(int i = 0; i < n; ++i)
        x[i] = 1.0;
}
void x_plus_sy(int n, double *x, double s, double *y, double *z){
    for(int i = 0; i < n; ++i)
        z[i] = x[i] + s * y[i];
}
void x_minus_Ay(int n, double *x, double *A, double *y, double *z){
    memcpy(z, x, n * sizeof(double));
    for(int i = 0; i < n; ++i)
        for(int j = 0; j < n; ++j)
            z[i] -= A[i*n+j] * y[j];
}
void x_minus_ATy(int n, double *x, double *A, double *y, double *z){
    memcpy(z, x, n * sizeof(double));
    for(int i = 0; i < n; ++i)
        for(int j = 0; j < n; ++j)
            z[i] -= A[j*n+i] * y[j];
}
void x_minus_sAy(int n, double *x, double s, double *A, double *y, double *z){
    memcpy(z, x, n * sizeof(double));
    for(int i = 0; i < n; ++i)
        for(int j = 0; j < n; ++j)
            z[i] -= s * A[i*n+j] * y[j];
}
void x_minus_sATy(int n, double *x, double s, double *A, double *y, double *z){
    memcpy(z, x, n * sizeof(double));
    for(int i = 0; i < n; ++i)
        for(int j = 0; j < n; ++j)
            z[i] -= s * A[j*n+i] * y[j];
}
double xTy(int n, double *x, double *y){
    double ret = 0.0;
    for(int i = 0; i < n; ++i)
        ret += x[i] * y[i];
    return ret;
}
double xTAy(int n, double *x, double *A, double *y){
    double ret = 0.0;
    for(int i = 0; i < n; ++i)
        for(int j = 0; j < n; ++j)
            ret += x[i] * A[i*n+j] * y[j];
    return ret;
}
void print_x(int n, double *x, const char *fmt, const char *prefix){
    if(prefix)
        printf("%s", prefix);
    for(int i = 0; i < n; ++i)
        printf(fmt, x[i]);
    putchar('\n');
}
void print_A(int n, double *A, const char *fmt, const char *prefix){
    if(prefix)
        printf("%s", prefix);
    for(int i = 0; i < n; ++i){
        for(int j = 0; j < n; ++j)
            printf(fmt, A[i*n+j]);
        putchar('\n');
    }
}

bool bicg(int n, double *A, double *x, double *b, int maxattempts, int maxiters, double tol){
    int i, j, k;
    double *x_, *ri, *rj, *ri_, *rj_, *p, *p_, *tmp; // j = i + 1, _ = transpose
    double alpha, beta;
    bool ret = false;

    // get parameters
    if(n * 10 > maxiters) maxiters = n * 10;
    // allocate
    x_ = (double*)malloc(n * sizeof(double));
    ri = (double*)malloc(n * sizeof(double));
    rj = (double*)malloc(n * sizeof(double));
    ri_ = (double*)malloc(n * sizeof(double));
    rj_ = (double*)malloc(n * sizeof(double));
    p = (double*)malloc(n * sizeof(double));
    p_ = (double*)malloc(n * sizeof(double));
    tmp = (double*)malloc(n * sizeof(double));
    // try different initial guess several times
    for(i = 0; i < maxattempts; ++i){
        // printf("Attempt %d\n", i + 1);
        // initialize
        rands(n, x, 2.0, -1.0);
        memcpy(x_, x, n * sizeof(double));
        x_minus_Ay(n, b, A, x, ri);
        x_minus_ATy(n, b, A, x, ri_);
        memcpy(p, ri, n * sizeof(double));
        memcpy(p_, ri_, n * sizeof(double));
        // start iteration
        for(j = 0; j < maxiters; ++j){
            alpha = xTy(n, ri_, ri) / xTAy(n, p_, A, p);
            x_plus_sy(n, x, alpha, p, x);
            // check x
            x_minus_Ay(n, b, A, x, tmp); // residual
            if((k = all_within(n, tmp, -tol, tol)) > 0) break; // NaN or all within
            // update
            x_plus_sy(n, x_, alpha, p_, x_);
            x_minus_sAy(n, ri, alpha, A, p, rj);
            x_minus_sATy(n, ri_, alpha, A, p_, rj_);
            beta = xTy(n, rj_, rj) / xTy(n, ri_, ri);
            x_plus_sy(n, rj, beta, p, p);
            x_plus_sy(n, rj_, beta, p_, p_);
            memcpy(ri, rj, n * sizeof(double));
            memcpy(ri_, rj_, n * sizeof(double));
        }
        // check
        // if(j == maxiters)   puts("Cannot converge.");
        // else if(k == 2)     puts("Detect NaN.");
        // else{               printf("Converge at iteration %d.\n", j+1); ret = true; break; }
        if(j != maxiters && k != 2){ ret = true; break; }
    }
    // release
    free(x_); free(ri); free(rj); free(ri_); free(rj_); free(p); free(p_); free(tmp);
    return ret;
}
bool gaussian(int n, double *A, double *x, double *b){
    double *A_, *b_; // temp A and b to change later
    int pr, pc, pi; // pivot row, pivot col, pivot index
    int i, j, ii, jj;
    double maxabs, tmp;
    bool ret = true;

    // allocate 
    A_ = (double*)malloc(n * n * sizeof(double));
    b_ = (double*)malloc(n * sizeof(double));
    // initialize
    memcpy(A_, A, n * n * sizeof(double));
    memcpy(b_, b, n * sizeof(double));
    // To row encholon form
    for(pr = pc = 0; pr < n && pc < n; ){
        // find max abs
        pi = pr, maxabs = fabs(A_[pr*n+pc]);
        for(i = pr+1; i < n; ++i)
            if((tmp = fabs(A_[i*n+pc])) > maxabs)
                pi = i, maxabs = tmp;
        // start ellimination
        if(ISZERO(maxabs)){ // pivot is zero
            A_[pr*n+pc] = NAN; // mark as NAN, so we can recognize it later
            if(pc + 1 < n) ++pc;
            else break;
        }
        else{
            if(pi != pr){ // swap row if needed
                memcpy(x, A_+pr*n, n * sizeof(double)); // x as tmp array
                memcpy(A_+pr*n, A_+pi*n, n *sizeof(double));
                memcpy(A_+pi*n, x, n * sizeof(double));
                tmp = b_[pr], b_[pr] = b_[pi], b_[pi] = tmp;
            }
            for(i = pr + 1; i < n; ++i){ // elliminate
                tmp = -A_[i*n+pc] / A_[pr*n+pc];
                A_[i*n+pc] = 0.0;
                for(j = pc + 1; j < n; ++j)
                    A_[i*n+j] += A_[pr*n+j] * tmp;
                b_[i] += b_[pr] * tmp;
            }
            if(pc + 1 >= n || pr + 1 >= n) break;
            else{ ++pc; ++pr; }
        }
    }
    // retrieve x
    memcpy(x, b_, n * sizeof(double));
    while(pr >= 0 && pc >= 0){
        if(isnan(A_[pr*n+pc])){
            if(ISZERO(x[pc])){ // 0 x[pc] = 0
                printf("Has multiple solutions for x[%d].\n", pc);
                x[pc] = 0.0; // can be arbitrary, but set to zero to avoid updating others
            }
            else{ // 0 x[pc] = b
                puts("No solution.");
                ret = false;
                break;
            }
            --pc;
        }
        else if(ISZERO(A_[pr*n+pc])){ // just go up if zero, then become non-zero
            --pr;
        }
        else{ // update
            x[pc] /= A_[pr*n+pc];
            for(i = pr - 1; i >= 0; --i)
                x[i] -= A_[i*n+pc] * x[pc];
            --pc;
        }
    }
    // release
    free(A_); free(b_);
    return ret;
}

void load_data(const char *filename, Edge *&edgeList, Source *&srcList, int &bcType, double &n, int &nN, int &nN0, int &nE, int &nL, int &nNeq, int &nLeq){
    FILE *fp;
    int i, j, k;
    int a, b;
    double r;

    // open the file
    fp = fopen(filename, "r");
    if(fp == NULL){
        puts("Failed to load data.");
        return;
    }
    // load parameters
    fscanf(fp, "%d %d %d %d %lf", &nN, &nE, &nN0, &bcType, &n);
    if(nN <= 0 || nE <= 0 || nN <= 1 ||
        (bcType != 0 && bcType != 1) ||
        (n < 1.5 || n > 2.1)){
        puts("Wrong parameters in data");
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

void get_equations(Edge *edgeList, Source *srcList, int bcType, double n, int nN, int nN0, int nE, int nL, int nNeq, int nLeq, double *&incNode, double *&conNode, double *&incLoop, double *&conLoop){
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
        i = buffer.front(), buffer.pop(); // BFS
        // i = buffer.top(), buffer.pop(); // DFS
        for(j = 0; j < nN; ++j)
            if(adj[(ii = i*nN+j)].dir != 0 && !visited[j]){
                parent[j] = i, depth[j] = depth[i] + 1, visited[j] = true, buffer.push(j);
                visitedE[adj[ii].edge] = true;
            }
    }
    free(visited);
    // // test
    // for(i = 0; i < nN; ++i)
    //     printf("parent of node %d is node %d\n", i, parent[i]);
    // for(tmp = 0, i = 0; i < nE; ++i)
    //     if(!visitedE[i])
    //         ++tmp;
    // if(bcType == 0)
    //     printf("%d edges are not in BFS spanning tree (nLeq = %d)\n", tmp, nLeq);
    // else if(bcType == 1)
    //     printf("%d edges are not in BFS spanning tree (nLeq - %d = %d)\n", tmp, nN0-1, nLeq-nN0+1);

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
            conLoop[tmp++] = srcList[0].bc - srcList[i].bc;
            // pressure drop = -(P_target - P_start) = P_start - P_target
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

void compute_R_and_J(int nE, int nLeq, int nNeq, double *incLoop, double *conLoop, double *incNode, double *conNode, double *x, double n, double *xtmp, double *R, double *J){
    int i, j, offset;
    // Compute R
    offset = nLeq;
    for(i = 0; i < nE; ++i) // tmp for store terms related to x
        xtmp[i] = x[i] * pow(fabs(x[i]), n-1.0);
        // xtmp[i] = x[i] * fabs(x[i]); // n == 2
    for(i = 0; i < nLeq; ++i){
        R[i] = conLoop[i];
        for(j = 0; j < nE; ++j)
            R[i] += incLoop[i*nE+j] * xtmp[j];
    }
    for(i = 0; i < nNeq; ++i){
        R[offset + i] = conNode[i];
        for(j = 0; j < nE; ++j)
            R[offset + i] += incNode[i*nE+j] * x[j];
    }
    // Compute J
    offset = nLeq * nE;
    for(i = 0; i < nE; ++i)
        xtmp[i] *= n / x[i];
        // xtmp[i] = fabs(x[i]); // n == 2
    for(i = 0; i < nLeq; ++i)
        for(j = 0; j < nE; ++j)
            J[i*nE+j] = incLoop[i*nE+j] * xtmp[j];
    memcpy(J+offset, incNode, nNeq * nE * sizeof(double));
}
void solve(int nE, int nLeq, int nNeq, double *incLoop, double *conLoop, double *incNode, double *conNode, double *&x, double n, int maxiter, int maxattempts, double tol, double step){
    int i, j, k;
    double *R, *J, *xtmp, *dx;
    bool ret;

    // allocate
    R = (double*) malloc(nE * sizeof(double));
    J = (double*) malloc(nE * nE * sizeof(double));
    xtmp = (double*) malloc(nE * sizeof(double));
    dx = (double*) malloc(nE * sizeof(double));
    x = (double*) malloc(nE * sizeof(double));

    // solve by newton's method
    srand(time(NULL));
    for(i = 0; i < maxattempts; ++i){ // has multiple attempts
        printf("Attempt %d\n", i + 1);
        rands(nE, x, 2.0, -1.0);
        for(j = 0; j < maxiter; ++j){
            // get residual and jacobian
            compute_R_and_J(nE, nLeq, nNeq, incLoop, conLoop, incNode, conNode, x, n, xtmp, R, J);
            // get dx and update
            memcpy(xtmp, dx, n * sizeof(double)); // store preivous step
            ret = gaussian(nE, J, dx, R);
            // ret = bicg(nE, J, dx, R);
            if(!ret) // cannot solve dx, so just use previous step (or random step?)
                memcpy(dx, xtmp, n * sizeof(double));
            x_plus_sy(nE, x, -step, dx, x); // we change -R -> R, so step*dx -> -step*dx
            // check
            if((k = all_within(nE, R, -tol, tol)) > 0) break;
            // if(all_within(nE, R, -1e12, -1e12) != 0){ k = 3; break; } // too large
            // // debug
            // print_A(nE, J, "%.5lf ", "J = \n");
            // print_x(nE, R, "%.5lf ", "R = ");
            // print_x(nE, dx, "%.5lf ", "dx = "); 
            // print_x(nE, x, "%.5lf ", "x = ");
            // getchar();
        }
        print_x(nE, R, "%.5lf ", "R = ");
        print_x(nE, x, "%.5lf ", "x = ");
        if(j == maxiter || k == 3)  puts("Caonnot converge");
        else if(k == 2)             puts("Detect NaN");
        else{                       printf("Converge at iteration %d\n", j+1); break; }
    }
    free(R); free(J);
    free(xtmp); // why this has bug
    free(dx);
}


#ifdef USE_CUDA
// wrapper for cuda
void solve_cuda(int nE, int nLeq, int nNeq, double *incLoop, double *conLoop, double *incNode, double *conNode, double *&x, double n, int maxiter = 100000, int maxattempts = 1000, double tol = 1e-6, double step = 0.1){
    solve(nE, nLeq, nNeq, incLoop, conLoop, incNode, conNode, x, n, maxiter, maxattempts, tol, step);
}
#endif // USE_CUDA