#include "hydra.h"

int main(int argc, char *argv[]){
    Edge *edgeList;
    Source *srcList;
    int bcType;
    double n;
    int nN, nN0, nE, nL, nNeq, nLeq;
    double *incNode, *conNode, *incLoop, *conLoop;
    double *x;
    int i, j;

    // check input
    if(argc != 2){
        puts("Wrong input arguments");
        return 1;
    }
    // load data
    load_data(argv[1], edgeList, srcList, bcType, n, nN, nN0, nE, nL, nNeq, nLeq);
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
    get_equations(edgeList, srcList, bcType, n, nN, nN0, nE, nL, nNeq, nLeq, incNode, conNode, incLoop, conLoop);
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
    print_x(nE, x, "%.5lf ", "x = ");
    // release
    free(x);
    return 0;
}


// int main(int argc, char *argv[]){
//     int n = argc == 2? atoi(argv[1]) : 3;
//     double *A, *x, *b;

//     // allocate
//     A = (double*)malloc(n * n * sizeof(double));
//     x = (double*)malloc(n * sizeof(double));
//     b = (double*)malloc(n * sizeof(double));

//     // initialize
//     rands(n * n, A, 2.0, -1.0);
//     rands(n, b, 2.0, -1.0);
//     print_A(n, A, "%.8lf, ", "A = \n");
//     print_x(n, b, "%.8lf, ", "b = \n");

//     // // solve by bicg
//     bicg(n, A, x, b);
//     print_x(n, x, "%.8lf, ", "x = \n");

//     // solve by gaussian
//     gaussian(n, A, x, b);
//     print_x(n, x, "%.8lf, ", "x = \n");

//     // release
//     free(A); free(x); free(b);

//     return 0;
// }