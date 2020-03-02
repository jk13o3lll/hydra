#include "hydra.h"

int main(int argc, char *argv[]){
    Edge *edgeList;
    Source *srcList;
    int bcType;
    double n;
    int nN, nN0, nE, nL, nNeq, nLeq; // nN include src nodes
    double *incNode, *conNode, *incLoop, *conLoop; // incidence matrix and constant vectors of node and loop equations (constant at the same side with unknowns)
    double *x;
    int i, j;

    if(argc != 2){ fprintf(stderr, "Wrong input arguments\n"); return 1; }

    // load data
    loadData(argv[1], edgeList, srcList, bcType, n, nN, nN0, nE, nL, nNeq, nLeq);
    printf("nE = %d, nN = %d, nN0 = %d, nL = %d\n", nE, nN, nN0, nL);
    printf("nNeq = %d, nLeq = %d\n", nNeq, nLeq);
    printf("n = %lf\n", n);
    printf("edgeList =\n");
    for(i = 0; i < nE; ++i)
        printf("%d %d %lf\n", edgeList[i].a, edgeList[i].b, edgeList[i].r);
    printf("srcList =\n");
    for(i = 0; i < nN0; ++i)
        printf("%d %lf\n", srcList[i].node, srcList[i].bc);

    // get equations
    getEquations(edgeList, srcList, bcType, n, nN, nN0, nE, nL, nNeq, nLeq, incNode, conNode, incLoop, conLoop);
    free(edgeList); free(srcList);
    for(i = 0; i < nNeq; ++i){ // show node equations
        printf("Node eq %d:", i);
        for(j = 0; j < nE; ++j)
            printf(" %.2lf", incNode[i*nE+j]);
        printf(", con = %.2lf\n", conNode[i]);
    }
    for(i = 0; i < nLeq; ++i){ // show loop equations
        printf("Loop eq %d:", i);
        for(j = 0; j < nE; ++j)
            printf(" %.2lf", incLoop[i*nE+j]);
        printf(", con = %.2lf\n", conLoop[i]);
    }

    // solve flowrate on each pipe
    printf("Press <Enter> to start."); getchar();
    solve(nE, nLeq, nNeq, incLoop, conLoop, incNode, conNode, x, n);
    free(incLoop); free(conLoop); free(incNode); free(conNode);

    // plot results
    printf("x = "); print_x(nE, x); putchar('\n');

    // release
    free(x);
    return 0;
}



// int main(int argc, char *argv[]){
//     int n = argc == 2? atoi(argv[1]) : 5;
//     double *A, *x, *b, *tmp;
//     // allocate
//     A = (double*) malloc(n * n * sizeof(double));
//     x = (double*) malloc(n * sizeof(double));
//     b = (double*) malloc(n * sizeof(double));
//     tmp = (double*) malloc(n * sizeof(double));
//     // init
//     srand(time(NULL));
//     // randA(n, A);
//     rand_x(n * n, A, 2.0, -1.0);
//     rand_x(n, b, 2.0, -1.0);
//     // solve Ax = b
//     gaussian(n, A, b, x);
//     // bicg(n, A, b, x);
//     x_minus_Ay(n, b, A, x, tmp);
//     // print
//     printf("A = \n"); print_A(n, A); putchar('\n');
//     printf("b = \n"); print_x(n, b); putchar('\n');
//     printf("x = \n"); print_x(n, x); putchar('\n');
//     printf("b-Ax = \n"); print_x(n, tmp); putchar('\n');
//     // release
//     free(A); free(x); free(b); free(tmp);
//     return 0;
// }