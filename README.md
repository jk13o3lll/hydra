# hydra
Pipe network analysis of hydraulics network through both sequential and parallel algorithms.

## Good resource

* https://en.wikipedia.org/wiki/Pipe_network_analysis
* https://www.researchgate.net/publication/321294545_Independent_Loops_Selection_in_a_Hydraulic_Looped_Network

## TODOs

### I. Build sequential version

1. Modify network
    * For multiple source node, add one additional node and construct fictious loops by connecting source nodes together to the new node.
1. Build [spanning tree](https://en.wikipedia.org/wiki/Spanning_tree) by [BFS](https://en.wikipedia.org/wiki/Parallel_breadth-first_search) and obtain discarded edges.
1. Put every edge back to the spannnig tree alternatively, and then find the loop by starting from two vertices on the edge and tracing back to their parents until their [lowest common ancestor](https://en.wikipedia.org/wiki/Lowest_common_ancestor).
    * From the discarded edges, every edge that is put back will form an indepedent loop because every loop contain at least one unique edge (discarded edge).
1. Construct matrix by buondary conditions (soruce node), node equations (mass conservation) and loop equations (energy conservation).
1. Use [Newton–Raphson method](https://en.wikipedia.org/wiki/Newton%27s_method) to solve the unknown.
    * Use BiCG or [BiCGSTAB](https://en.wikipedia.org/wiki/Biconjugate_gradient_stabilized_method) to solve the matrix inverse (gradient).

### II. Build parallel versoin

After validating sequential version, write a CUDA version.
