# hydra
Pipe network analysis of hydraulics network through both sequential and parallel algorithms.

## Good resource

* https://en.wikipedia.org/wiki/Pipe_network_analysis
* https://www.researchgate.net/publication/321294545_Independent_Loops_Selection_in_a_Hydraulic_Looped_Network

## TODOs

### I. Build sequential version

1. Modify network
    * For multiple source node, add one additional node and construct fictious loops by connecting source nodes together to the new node.
    * Merge internal loop between two nodes.
1. Build [spanning tree](https://en.wikipedia.org/wiki/Spanning_tree) by [BFS](https://en.wikipedia.org/wiki/Parallel_breadth-first_search) and obtain discarded edges.
    * root of BFS start from node has highest degree (most edges connect to) (this is efficient way to decrease total depth or diameter) (or [BFS for every](https://codeforces.com/blog/entry/7372)(can parallel) or [other1](https://www.sciencedirect.com/science/article/pii/0020025595001352) or [other2](https://www.researchgate.net/publication/220617691_Minimum_Diameter_Spanning_Trees_and_Related_Problems))
1. Put every edge back to the spannnig tree alternatively, and then find the loop by starting from two vertices on the edge and tracing back to their parents until their [lowest common ancestor](https://en.wikipedia.org/wiki/Lowest_common_ancestor).
    * From the discarded edges, every edge that is put back will form an indepedent loop because every loop contain at least one unique edge (discarded edge).
1. Construct matrix by buondary conditions (soruce node), node equations (mass conservation) and loop equations (energy conservation).
1. Use [Newton–Raphson method](https://en.wikipedia.org/wiki/Newton%27s_method) to solve the unknown.
    * Use BiCG or [BiCGSTAB](https://en.wikipedia.org/wiki/Biconjugate_gradient_stabilized_method) to solve the matrix inverse (gradient).

### II. Build parallel versoin

After validating sequential version, write a CUDA version.

### Notes

1. If using flowrate as boundary condition, both input and output flow should be provide (one output will be depedent, so one node equation contain one output should be removed) (pressure and radius are decided implicitly) In the end you will get N-1 node eq & E-(N-1) loop eq (1 is sum in = sum out).
2. If using pressure as boundary condition, every additional inputs and outputs should add one edge, which means we need (inputs+outputs) indepedent loop equations (just select one input and go through the path contain other inputs and outpus.). In the end you will get N node eq (sum in and sum out is include in source) + L loop eq + (S-1) source loop eq (node and sourve are different here). (For one source connect to two nodes, it can be separate into two clone nodes).
