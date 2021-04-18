

<p style="text-align:left; font-weight:bold;">Course Name: Algorithms, Part I</p>

<a style="text-align:left; font-weight:bold;" href = "https://www.coursera.org/learn/algorithms-part1"> Course Website: Coursera</a>

<p style="text-align:left; font-weight:bold;">Instructor: Professor Kevin Wayne, Professor Robert Sedgewick</p>

<p style="text-align:left; font-weight:bold;">Books: Algorithms(4th edition)</p>



# Week 1

1. Introduction

    * Course Structure:
        ![](/images/Algorithm/CourseStructure.png)

2. Analysis of Algorithms

    * Introduction
        * Running time is the key concern
        * Reasons to analyze algorithms:
            * Predict performance
            * Compare algorithms
            * Provide guarantees
            * Understand theoretical basis
        * Scientific method applied to analysis of algorithms
            * observe => hypothesize => predict => verify => validate
            * Experiments must be reproducible, hypotheses must be falsifiable
    * Observations
        * Use plots to analyze algorithms:
            * Standard plot: `(problem size N, running time T(n))`
            * Log-log plot: `(lgN, lgT(N))`
        * Factors affecting algorithm performance
            * System independent effects: algorithm, input
            * System dependent effects: hardware, software, system
    * Mathematical Models
        * Total running time: sum of cost * frequency for all operations
        * Simplification
            * Cost model: use some basic operation as a proxy for running time
            * Tilde notation: estimate running time as a function of input size N, and ignore lower order terms
                $$1 + 2 + ... + N  = \sum_{i = 1}^{N} i  \sim \int_{x = 1}^{N} xdx \sim  \frac{1}{2}N^{2}$$
                $$1 + \frac{1}{2} + ... + \frac{1}{N}  = \sum_{i = 1}^{N} \frac{1}{i}  \sim \int_{x = 1}^{N} \frac{1}{x}dx \sim  lnN$$
                $$\sum_{i = 1}^{N}\sum_{j = 1}^{N}\sum_{k = 1}^{N} 1  \sim \int_{x = 1}^{N} \int_{y = 1}^{N}\int_{z = 1}^{N}dxdydz \sim \frac{1}{6}N^{3}$$
        * In general, $T_{N} = c_{1}A +c_{1}A + c_{2}B + c_{3}C + c_{4}D + c_{5}E $, where A = array access, B = integer add, C = integer compare, D = increment, E = variable assignment
    * Order of Growth Classification
        * Order-of-growth of typical algorithms: $1(constant)$, $logN(logarithmic)$, $N(linear)$, $NlogN(linearithmic)$, $n^{2}(quadratic)$, $N^{3}(cubic)$, $2^{N}(exponential)$
    * Theory of Algorithms
        * Types of analyses
            * Best case: lower bound on cost, determined by easiest input, provides a goal for all inputs
            * Worst case: upper bound on cost, determined by most difficult input, provides a guarantee for all inputs
            * Average case: expected cost for random input, need a model for random input, provides a way to predict performance
        * Theory of algorithms
            * Goals: establish difficulty of a problem, develop optimal algorithms
            * Approach: suppress details in analysis, eliminate variability in input by focusing on the worst case
            * Optimal algorithm: performance guarantee for any input, no algorithm can provide a better performance guarantee
            * Commonly-used notations:
                * Big Theta: asymptotic order of growth, $\Theta(f(N))$
                * Big Oh: $\Theta(f(N))$ and smaller, $O(f(N))$
                * Big Omega: $\Theta(f(N))$ and larger, $\Omega(f(N))$
    * Memory
        * Basics
            * Bit: 0,1; byte: 8 bits; MB: 1million or $2^{30}$ bytes; GB: 1 billion or $2^{30}$ bytes
            * Old machine: a 32-bit machine with 4 byte pointers
            * Modern machine: a 64-bit machine with 8 byte pointers
        * Object memory usage:
            * boolean: 1 byte; byte: 1 byte; char: 2 bytes; int: 4 bytes; float: 4 bytes; long: 8 bytes; double: 8 bytes
            * One-dimensional arrays: char[]: 2N + 24 bytes; int[]: 4N + 24 bytes; double[]: 8N + 24 bytes
            * Two-dimensional arrays: char[][]: ~ 2MN bytes; int[][]: ~ 4MN bytes; double[][]: ~ 8MN bytes
            * Object overhead: 16 bytes
            * Reference: 8 bytes
            * Padding: each object uses a multiple of 8 bytes

3. Union Find

    * Dynamic Connectivity

        * Problem description
            * Given a set of N objects, we can union two objects by connecting them. Two objects are connected if there is a path through which you can travel from one point to another. What we want to know is that given two objects in the set, whether they are connected or not.
        * Application: pixels in a digital photo, computers in a network, friends in a social network, variable names in a program, etc.
        * Equivalence relation of the connectivity:
            * Reflexive: p is connected to p
            * Symmetric: if p is connected to q, then q is connected to p
            * Transitive: if p is connected to q, and q is connected to r, then p is connected to r
        * Connected components: maximal set of objects that are mutually connected

    * Quick Find

        * It's an eager algorithm to solve the dynamic connectivity problem
        * Data structure description
            * The data structure is an integer array `id[]` of size N， each entry is the id that represents the connected component to which this point belongs
            * Connected command: if two points have the same id, then they are connected
            * Union command: when union point `p` to point `q`, change all elements whose id equals to `id[p]` to `id[q]`
        * Algorithm

         ```java
         public class QuickFindUF {
           private int[] id;
        
           public QuickFindUF(int N) {
             id = new int[N];
             for (int i = 0; i < N; i++) {
               id[i] = i;
             }
           }
        
           public boolean connected(int p, int q) {
             return id[p] == id[q];
           }
        
           public void union(int p, int q) {
             int idp = id[p];
             int idq = id[q];
             for (int i = 0; i < N; i++) {
               if (id[i] == idp) {
                 id[i] = idq;
               }
             }
           }
         }
         ```

         * Performance analysis: N union objects on N objects takes quadratic time, which is quite unbearable fo large datasets

     * Quick Union

        * It's a lazy approach  to solve the dynamic connectivity problem
        * Data structure description
            * The data structure is an integer array `id[]` of size N, each entry is the id that represents the parent point of the current point. All points construct a forest containing connected trees and unconnected leaves
            * Root of the tree containing point `i` is `id[id[id[...id[i]...]]]`
            * Connected command: if two points have the same root, then they are connected
            * Union command: when union point `p` to point `q`, set the parent of p's root be q's root
        * Algorithm

         ```java
         public class QuickUnionUF {
           private int[] id;
        
           public QuickUnionUF(int N) {
             for (int i = 0; i < N; i++) {
               id[i] = i;
             }
           }
        
           private int root(int node) {
             while(node != id[node]) {
               node = id[node];
             }
             return node;
           }
          
           public boolean connected(int p, int q) {
             return root(p) == root(q);
           }
        
           public void union(int p, int q) {
             int rootp = root(p);
             int rootq = root(q);
             id[rootp] = rootq;
           }
         }
         ```

         * Performance analysis: trees can be too tall, then find command can be slow since it's a linear operation

      * Quick Union Improvement

        * Improvement 1: Weighting

            * Keep track of size of each tree, and balance by linking root of smaller tree to root of larger tree in union command. When there are equal size trees, put the former tree under the latter tree
            * Data structure description: same as quick union, but maintain extra array `size[i]` to count number of objects in the tree rooted at `i`
            * Find command: identical to quick union
            * Union link root of smaller tree to root of larger tree, update the size[] array

            ```java
             public void union(int p, int q) {
             int rootp = root(p);
             int rootq = root(q);
             if (size(rootp) <= size(rootq)) {
               id[rootp] = rootq;
               size[rootq] += size[rootp];
             } else {
               id[rootq] = rootp;
               size[rootp] += size[rootq];
             }
             }
            ```

            * Theorem: depth of any node x is at most lgN
            * Performance analysis
                ![](/images/Algorithm/UFPerformance.png)

        * Improvement 2: Path Compression

            * Idea: after computing the root of p, set the id of each examined node to point to that root, so that the tree is flattened
                ![](/images/Algorithm/PathCompression1.png)
                 Convert the above tree into the following one
                 ![](/images/Algorithm/PathCompression2.png)

            * Java implementation 

                * Two-pass implementation: add a second loop to `root()` to set the `id[]` of each examined node to the root
                * Simpler one-pass variant: make every other node in path point to its grandparent

                 ```java
                 private int root(int node) {
                   while (node != id[node]) {
                     id[node] = id[id[node]];
                     node = id[node];
                   }
                   return node;
                 }
                 ```

             * Theorem: starting from an empty data structure, any sequence of M union-find operations on N objects makes $\leq c(N + M lgN)$ array accesses. Actually there is no linear algorithm for the union find problem

         * Summary
            ![](/images/Algorithm/UFSummary.png)

       * Union Find Application

         * There are many applications for the union find problem

         * Percolation

             * It's a model for many physical systems

             * There is N-by-N grid of sites, each site is open with probability $p$ or blocked with probability $1 - p$. The system percolates iff the top and bottom are connected by open sites
                 ![](/images/Algorithm/Percolation.png)

             * Likelihood of percolation

                 * When N is large, theory guarantees a sharp threshold $p^{*}$

                     If $p \ge p^{*}$: almost certainly percolates

                     If $p \le p^{*}$: almost certainly does not percolate
                      ![](/images/Algorithm/ThresholdPercolation.png)

                 * No analytical solution of $p^{*}$ could be found, but it could be estimated using Monte Carlo simulation

                     1. Initialize N-by-N whole grid to be blocked
                     2. Declare random sites open until top connected to bottom
                     3. Vacancy percentage estimates $p^{*}$

                 * We can use dynamic connectivity to estimate the percolation threshold. For $N = 100$, the estimated $p^{*} \approx 0.592746$