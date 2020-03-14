#ifndef __HYDRA_H__
#define __HYDRA_H__

// We currently haven't support to have multiple sources on single node.
// You have to merge it by yourself
// (For bc is flow, ex. Q = Q1+Q2)
// (it is impossible to have multiple pressure at the same node)

// We currently haven't support parallel edges.
// That is multiple edges between two nodes.
// So you have to merge them to eqivalent edge,
// then retrieve individual results based on result of the equivalent edge.

// We currently haven't support simplify long serial edges.
// That is, all nodes between start and end are indegree = 1 and outdegree = 1.
// So you can merge them to eqivalent edge,
// then retrieve individual results based on result of the equivalent edge.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <queue>
// #include <stack>

#define ISZERO(x) (fabs((x)) < 1e-12)

extern "C" {
struct Edge{
    int a, b;  // from node a to node b
    double r;  // r value is property of the pipe (constant) (r > 0)
};
struct Source{
    int node;   // on which node (id of node): 0, 1, ..., nN-1
    double bc;  // value of boundary condition
};
struct Adjacency{
    int edge; // id of edge: 0, 1, ..., nE-1 (-1 is no edge)
    int dir;  // direction: +1 same as flow (Q), -1 opposite, 0 no edge
};
}

// n: size; s: scale;  x, y, z: vector; A: matrix (square); T: transpose;
int all_within(int n, double *x, double lowerb, double upperb);
void rands(int n, double *x, double a = 1.f, double b = 0.f); // x = [(0~1)*a+b, ...]
void zeros(int n, double *x); // x = [0, 0, ...]
void ones(int n, double *x); // x = [1, 1, ...]
void x_plus_sy(int n, double *x, double s, double *y, double *z); // z = x + sy
void x_minus_Ay(int n, double *x, double *A, double *y, double *z); // z = x - Ay
void x_minus_ATy(int n, double *x, double *A, double *y, double *z); // z = x - ATy
void x_minus_sAy(int n, double *x, double s, double *A, double *y, double *z); // z = x - sAy
void x_minus_sATy(int n, double *x, double s, double *A, double *y, double *z); // z = x - sATy
double xTy(int n, double *x, double *y); // return = xTy
double xTAy(int n, double *x, double *A, double *y); // return = xTAy
void print_x(int n, double *x, const char *fmt = "%.5lf ", const char *prefix = NULL);
void print_A(int n, double *A, const char *fmt = "%.5lf ", const char *prefix = NULL);

// biconjugate gradient (solve x, for Ax = b)
bool bicg(int n, double *A, double *x, double *b,
    int maxattempts = 100, int maxiters = 1000, double tol = 1e-6);
// gaussian ellimination (LU with partial pivoting)
bool gaussian(int n, double *A, double *x, double *b);

// Load network data
// n: constant for pressure loss equation (not size!!!)
// nN: num of nodes (include source node)
// nN0: num of source nodes
// nE: num of edge
// nL: num of loops (indepent loop, haven't include loop for source (bcType==1))
// nNeq: num of node equations
// nLeq: num of loop equations
void load_data(const char *filename, Edge *&edgeList, Source *&srcList, int &bcType, double &n, int &nN, int &nN0, int &nE, int &nL, int &nNeq, int &nLeq);
// Get incidence matrix (value is +-r value) and coefficients (value is from bc) of node equations and loop equations
// incidence matrix and coefficient vectors of node and loop equations (coefficient at the same side with unknowns)
void get_equations(Edge *edgeList, Source *srcList, int bcType, double n, int nN, int nN0, int nE, int nL, int nNeq, int nLeq, double *&incNode, double *&conNode, double *&incLoop, double *&conLoop);

// Compute residual and jacobian (required for solving dx)
void compute_R_and_J(int nE, int nLeq, int nNeq, double *incLoop, double *conLoop, double *incNode, double *conNode, double *x, double n, double *xtmp, double *R, double *J);
// Use Newton Raphson method (nolinear equations) to solve
// with higer dim, you should use smaller step size and larger tol
void solve(int nE, int nLeq, int nNeq, double *incLoop, double *conLoop, double *incNode, double *conNode, double *&x, double n, int maxiter = 100000, int maxattempts = 1000, double tol = 1e-6, double step = 0.1);


#ifdef HYDRA_USE_CUDA
#include "hydra_cuda.h"
// wrapper for cuda
void solve_cuda(int nE, int nLeq, int nNeq, double *incLoop, double *conLoop, double *incNode, double *conNode, double *&x, double n, int maxiter = 100000, int maxattempts = 1000, double tol = 1e-6, double step = 0.1);
#endif // USE_CUDA

#endif