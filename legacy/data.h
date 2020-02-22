#ifndef __DATA_H__
#define __DATA_H_

#if !defined(NETWORK1) && !defined(NETWORK2) && !defined(NETWORK3) && !defined(NETWORK4)
    #define NETWORK1 // network 1 as default
#endif

// basic data type
struct Edge {
    int a, b;
    double r;   // r is constant related to pipe diameter, pipe length, and friction factors
};
struct SourceNode {
    int id;
    double x;   // flow or pressure
};

// network 1 (boundary condition by flowrate)
#ifdef NETWORK1
const int bcType = 0; // boundary condition is flowrate to node
const double n = 2.;    // constant in pressure loss eqn, 1.5 < n < 2.0
const Edge edgeList[] = { 
    {0, 1, 5.}, // a, b, r (from a to b, with pipe radius = r)
    {0, 2, 1.}, // all node should be contain in single graph
    {1, 2, 1.}, // (should not have multple separate graphs)
    {1, 3, 1.}, // id of node should start from zero and continuous
    {2, 3, 4.}  // use original direction as assumed direction to
};              // generate adjacency matrix (adj store +1 or -1) and incidence matrix
const SourceNode srcList[] = {  // direction for node eqation is 
    {0, 10.},                   // flow into node is positive
    {3, -10.}                   // flow out from node is negative
};                              // remember to include all inlet and all outlet (should be mass conservation)
#endif // NETWORK1

// network 2 (boundary condition by flowrate)
#ifdef NETWORK2
const int bcType = 0; // boundary condition is flowrate to node
const double n = 2.;
const Edge edgeList[] = { 
    {0, 1, 1.},
    {1, 6, 5.},
    {1, 2, 1.},
    {0, 3, 5.},
    {3, 4, 3.},
    {2, 4, 1.},
    {2, 5, 1.},
    {5, 6, 3.},
    {3, 8, 3.},
    {4, 7, 1.},
    {5, 7, 1.},
    {7, 8, 2.},
    {8, 9, 3.},
    {6, 9, 3.},
};
const SourceNode srcList[] = {
    {0, 10.},
    {9, -10.}
};
#endif // NETWORK2

// network 3 (boundary condition by pressure)
#ifdef NETWORK3
const int bcType = 1; // boundary condition is pressure at node
const double n = 2.;
const Edge edgeList[] = {
    {0, 1, 1.},
    {1, 6, 5.},
    {1, 2, 1.},
    {0, 3, 5.},
    {3, 4, 3.},
    {2, 4, 1.},
    {2, 5, 1.},
    {5, 6, 3.},
    {3, 14, 1.},
    {4, 7, 1.},
    {5, 7, 1.},
    {7, 8, 2.},
    {8, 9, 3.},
    {6, 9, 3.},
    {8, 10, 1.},
    {14, 11, 2.},
    {13, 16, 5.},
    {14, 15, 5.},
    {10, 11, 3.},
    {11, 12, 2.},
    {9, 16, 1.},
    {12, 15, 2.},
    {10, 13, 2.},
    {13, 12, 1.},
    {16, 17, 1.},
    {15, 17, 2.},
    {18, 0, 1.}, // edge connect to source
    {9, 19, 1.}, // edge connect to source
    {17, 20, 1.} // edge connect to source
};
const SourceNode srcList[] = {  // direction for loop equation is
    {18, 10.},                  // same as loop direction is positive
    {19, 0.},                   // opposite to loop direction is negative
    {20, 0.}    // source node should be exclude when construct node eq
};
#endif // NETWORK3

// network 4 (boundary condition by pressure)
#ifdef NETWORK4
const int bcType = 1; // boundary condition is pressure at node
const double n = 2.;
const Edge edgeList[] = {
    {0, 1, 1.},
    {1, 6, 5.},
    {1, 2, 1.},
    {0, 3, 5.},
    {3, 4, 3.},
    {2, 4, 1.},
    {2, 5, 1.},
    {5, 6, 3.},
    {3, 14, 1.},
    {4, 7, 1.},
    {5, 7, 1.},
    {7, 8, 2.},
    {8, 9, 3.},
    {6, 9, 3.},
    {8, 10, 1.},
    {14, 11, 2.},
    {13, 16, 5.},
    {14, 15, 5.},
    {10, 11, 3.},
    {11, 12, 2.},
    {9, 16, 1.},
    {12, 15, 2.},
    {10, 13, 2.},
    {13, 12, 1.},
    {16, 17, 1.},
    {15, 17, 2.},
    {18, 0, 1.}, // edge connect to source
    {9, 19, 1.}, // edge connect to source
    {17, 20, 1.} // edge connect to source
};
const SourceNode srcList[] = {
    {18, 100.},
    {19, 0.},
    {20, 0.}
};
#endif // NETWORK4

#endif