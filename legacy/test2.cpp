// g++ test2.cpp -o test.exe -fopenmp -lpthread

#include "hydra2.h"

int main(int argc, char *argv[]){
    Edge *edgeList;
    Source *srcList;
    int bcType;
    long double n;
    int nN, nN0, nE, nL, nNeq, nLeq; // nN include src nodes
    long double *incNode, *conNode, *incLoop, *conLoop; // incidence matrix and constant vectors of node and loop equations (constant at the same side with unknowns)
    long double *x;
    int i, j;

    if(argc != 2){ fprintf(stderr, "Wrong input arguments\n"); return 1; }

    // load data
    loadData(argv[1], edgeList, srcList, bcType, n, nN, nN0, nE, nL, nNeq, nLeq);
    printf("nE = %d, nN = %d, nN0 = %d, nL = %d\n", nE, nN, nN0, nL);
    printf("nNeq = %d, nLeq = %d\n", nNeq, nLeq);
    printf("n = %Lf\n", n);
    printf("edgeList =\n");
    for(i = 0; i < nE; ++i)
        printf("%d %d %Lf\n", edgeList[i].a, edgeList[i].b, edgeList[i].r);
    printf("srcList =\n");
    for(i = 0; i < nN0; ++i)
        printf("%d %Lf\n", srcList[i].node, srcList[i].bc);

    // get equations
    getEquations(edgeList, srcList, bcType, n, nN, nN0, nE, nL, nNeq, nLeq, incNode, conNode, incLoop, conLoop);
    free(edgeList); free(srcList);
    for(i = 0; i < nNeq; ++i){ // show node equations
        printf("Node eq %d:", i);
        for(j = 0; j < nE; ++j)
            printf(" %.2Lf", incNode[i*nE+j]);
        printf(", con = %.2Lf\n", conNode[i]);
    }
    for(i = 0; i < nLeq; ++i){ // show loop equations
        printf("Loop eq %d:", i);
        for(j = 0; j < nE; ++j)
            printf(" %.2Lf", incLoop[i*nE+j]);
        printf(", con = %.2Lf\n", conLoop[i]);
    }

    // solve flowrate on each pipe
    printf("Press <Enter> to start."); getchar();
    solve(nE, nLeq, nNeq, incLoop, conLoop, incNode, conNode, x, n);
    free(incLoop); free(conLoop); free(incNode); free(conNode);

    // plot results
    printx("x = ", nE, x);

    // release
    free(x);
    return 0;
}

// int main(int argc, char *argv[]){
//     int n = argc == 2? atoi(argv[1]) : 5;
//     long double *A, *x, *b, *tmp;
//     // allocate
//     A = (long double*) malloc(n * n * sizeof(long double));
//     x = (long double*) malloc(n * sizeof(long double));
//     b = (long double*) malloc(n * sizeof(long double));
//     tmp = (long double*) malloc(n * sizeof(long double));
//     // init
//     srand(time(NULL));
//     // randA(n, A);
//     randA(n, A, 2.0, -1.0);
//     randx(n, b, 2.0, -1.0);
//     // solve Ax = b
//     // bicgstab(n, A, b, x, false, 1e-3); // sometimes cannot converge
//     // bicgstab(n, A, b, x, true, 1e-3);  // preconditioned not work?
//     gaussian(n, A, b, x);           // this work well
//     xminusAy(n, b, A, x, tmp);
//     // print
//     printA("A = \n", n, A);
//     printx("b = \n", n, b);
//     printx("x = \n", n, x);
//     printx("b-Ax = \n", n, tmp);
//     // release
//     free(A); free(x); free(b); free(tmp);
//     return 0;
// }
