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





# Week 2

1. Stacks and Queues

    * Fundamental data types:

        * Value: collection of objects
        * Operations: insert, remove, iterate, test if empty
        * Stack: LIFO, insert is called push, remove is called pop
        * Queue: FIFO, insert is called enqueue, remove is called dequeue
        * Module programming: completely separate interface and implementation

    * Stack 

        * Stack API

            ```java
            public class StackOfStrings {
              StackOfStrings() 
              void push(String item)
              String pop()
              boolean isEmpty()
              int size()
            }
            ```

         * LinkedList Implementation

            ```java
            public class LinkedListStack {
              private Node first = null;
            
              private class Node {
                String item;
                Node next;
              }
            
              public String pop() {
                String item = first.item;
                first = first.next;
                return item;
              }
            
              public void push(String item) {
                Node oldFirst = first;
                first = new Node();
                first.item = item;
                first.next = oldFirst;
              }
            
              public boolean isEmpty() {
                return first == null;
              }
            }
            ```

            * Every operation takes constant time in the worst case

         * Array Implementation

            ```java
            public class ArrayStack {
              private String[] items;
              private int N = 0;
            
              public ArrayStack(int capacity) {
                items = new String[capacity];
              }
            
              public boolean isEmpty() {
                return N ==0;
              }
            
              public void push(String item) {
                items[N++] = item;
              }
              /* Loitering
              public String pop() {
                return items[--N];
              }
              */
              public String pop() {
                String item = items[--N];
                items[N] = null;
                return item; 
              }
            
            }
            ```

              * Don't address underflow(pop with an empty stack) and overflow(push with a full stack)
              * Loitering: holding a reference to an object when it is no longer needed

    * Resizing Arrays

        * First try: change size by 1 for push and pop

            * It's too expensive, you need to copy all items to the new array. In order to insert $N$ items into an empty stack, it costs $N^2$ time(To insert into a full stack with $k$ size, you need to create a $k + 1$ array and copy all items, so that's $1 + ... + N \sim N^2$)

        * Repeated doubling: double array size when full, halve array size when one-quarterfull

            ```java
            public class ArrayStack {
              private String[] items = new String[1];
              private N = 0;
            
              public void push(String item) {
                if (N == items.length) {
                  resize(2 * items.length);
                }
                items[N++] = item;
              }
            
              public String pop() {
                String item = items[--N];
                items[N] = null;
                if  (N > 0 && N == items.length /4) {
                  resize(items.length / 2)
                }
                return item;
              }
            
              private void resize(int capacity) {
                String[] copy = new String[capacity];
                for (int i = 0; i < N; i++) {
                  copy[i] = items[i];
                }
                items = copy;
              }
            }
            ```

              * To insert N items into an empty array, it costs $N + 2 + 4 + ... + N \sim 3N$
              * Amortized analysis: average running time per operation over a worst-case sequence of operations
              * Starting from an empty stack, any sequence of $M$ push and pop operations takes time proportional to $M$
            * LinkedList v.s. resizing array for stack implementation
                * LinkedList: every operation takes constant time in the worst case, use extra time and space to deal with the links
                * Resizing-array: every operation takes constant amortized time, less wasted space

    * Queues

        * Queue API

            ```java
            public class QueueOfStrings {
              QueueOfStrings()
              void enqueue(String item)
              String dequeue()
              boolean isEmpty()
              int size()
            }
            ```

         * LinkedList Implementation

            ```java
            public class LinkedListQueue {
              private Node first;
             private Node last;
            
              private class Node {
                String item;
                Node next;
             }
            
              public String dequeue() {
                String item = first.item;
                first = first.next;
                if (isEmpty()) {
                  last = null;
                }
                return item;
             }
            
              public void enqueue(String item) {
                Node oldLast = last;
                last = new Node();
                last.item = item;
                last.next = null;
                if (isEmpty()) {
                  first = last;
                } else {
                  oldLast.next = last;
                }
            
              public boolean isEmpty() {
                return first == null;
              }
            }
            ```

         * Resizing Array Implementation

            ```java
            public class ArrayQueue {
              private String[] items = new String[1];
              private int front = 0;
              private int rear = 0;
            
              public void enqueue(String item) {
                if (rear - front == items.length) {
                  resize(2 * items.length);
                }
                items[rear] = item;
                rear += 1;
             }
            
              public String dequeue() {
                String item = items[front];
                items[front] = null;
                front += 1;
                if ( rear - front > 0 && rear - front == items.length / 4) {
                  resize(items.length / 2);
                }
             }
            
              public void resize(int capacity) {
                String[] copy = new String[capacity];
                for (int i = 0; i < rear - front; i++) {
                  copy[i] = items[front + i];
                }
                items = copy;
              }
            }
            ```

     * Generics

        * Collections that can contain various types data

            * Implement a separate collection class for each type: tedious
            * Implement a collection for object type and use type casting like: Object variable = new Object(); T target = (T) variable
            * Implement a collection with generics: avoid casting when using collections, discover type mismatch errors at compile-time rather than run-time

        * LinkedList Stack with Generics

            ```java
            public class LinkedListStackGeneric<T> {
              private Node first = null;
            
              private class Node {
                T item;
                Node next;
              }
              public void push(T item) {
                Node oldFirst = first;
                first = new Node();
                first.item = item;
                first.next = oldFirst;
              }
            
              public T pop() {
                T item = first.item;
                first = first.next;
                return item;
              }
            
              public boolean isEmpty() {
                return first == null;
              }
            }
            ```

        * ArrayStack with Generics

            ```java
            public class ArrayStackGeneric<T> {
              private T[] items;
              private int N = 0;
            
              public ArrayStackGeneric(int capacity) {
                // It's not allowed in java to create generic arrays, you need type casting here
                s = (T[]) new Object[capacity];
              }
            
              public boolean isEmpty() {
                return N ==0;
              }
            
              public void push(T item) {
                items[N++] = item;
              }
            
              public T pop() {
                T item = items[--N];
                items[N] = null;
                return item;
              }
            }
            ```

            * Autoboxing
                * For each primitive type there is a wrapper object type, automatic cast ca be done between primitive type and wrapper type
                * You need to use wrapper types rather than primitive types in collections
                * Convert from String to other type: `int aInt = Integer.parseInt("aString");`

    * Iterators

        * Java supports iteration over stack items by client, without revealing the internal representation of the collections

        * Iterable: a interface that has a method that returns an Iterator

            ```java
            public interface Iterable<Item> {
              Iterator<Item> iterator();
            }
            ```

        * Iterator: an interface that has methods `hasNext()` and `next()`

            ```java
            public interface Iterator<Item> {
              boolean hasNext();
              Item next();
            }
            ```

        * Java supports elegant enhanced-for loop with iterable data structures

            ```java
            for(String element: stack) {
              System.out.println(element);
            }
            ```

        * LinkedList Stack with Iterator

            ```java
             
            public class LinkedListStack<T> implements Iterable<T> {
              ...
              public Iterator<T> iterator() {
                return new ListIterator();
              }
             
              private class ListIterator implements java.util.Iterator<T> {
                private Node current = first;
                    
                public boolean hasNext() {
                  return current != null;
                }
             
                public T next() {
                  T item = current.item;
                  current = current.next;
                  return item;
                }
              }
            }
            ```

         * Array List Stack with Iterator

            ```java
            public class ArrayStack<T> implements Iterable<T> {
              
              public Iterator<T> iterator() {
                return new ArrayIterator();
              }
            
              private class ArrayIterator implements Iterator<T> {
                private int i = N;
                
                public boolean hasNext() {
                  return i > 0;
                }
                public T next() {
                  return items[--i];
                }
              }
            }
            ```

         * Bag 

            * Application: adding items to a collection and iterating, order doesn't matter

            * Bag API

                ```java
                public class Bag<T> implements Iterable<T> {
                  Bag();
                  void add(T item);
                  int size();
                  Iterable<T> iterator();
                }
                ```

    * Stack and Queue Applications

        * Java collections library
            * List
                * List interface: `java.util.List: public interface List<E> extends Collection<E>`
                * Implementations: `java.util.ArrayList: public class LinkedList<E> extends AbstractSequentialList<E> implements List<E>, Deque<E>, Cloneable, java.io.Serializable`, `java.util.LinkedList: public class LinkedList<E> extends AbstractSequentialList<E> implements List<E>, Deque<E>, Cloneable, java.io.Serializable`
            * Stack: `java.util.Stack: public class Stack<E> extends Vector<E>`
            * Queue: `java.util.Queue: public interface Queue<E> extends Collection<E>`
        * Stack applications:
            * Parsing in a complier
                * Dijkstra's two-stack algorithm
                    * Operand: push into the operand stack
                    * Operator: push into the operator stack
                    * Left parenthesis: ignore
                    * Right parenthesis: pop operator and two operands, calculate, push the result into the operand stack
            * Java virtual machine
            * Undo in a word processor
            * Back button in a web browser
            * ...

    * Problems

        * Create a queue using two stack

            ```java
            // You should not use java.util.Stack, you should use java.util.Deque instead:
            // A more complete and consistent set of LIFO stack operations is provided by the  
            // Deque interface and its implementations, which should be used in preference to this class. 
            //Space complexity: O(N), time complexity: O(N)
            class MyQueue {
              Deque<Integer> out, in;
            
              public MyQueue() {
                  in = new ArrayDeque<>();
                  out = new ArrayDeque<>();
              }
            
              public void push(int x) {
                  while (!out.isEmpty()) {
                    in.addLast(out.pollLast());
                  }
                  in.addLast(x);
              }
              
              public int pop() {
                  while (!in.isEmpty()) {
                    out.addLast(in.pollLast());
                  }
                  return out.pollLast();
              }
              
              public int peek() {
                  while (!in.isEmpty()) {
                    out.addLast(in.pollLast());
                  }
                  return out.peekLast();
              }
              
              public boolean empty() {
                  return out.isEmpty() && in.isEmpty();
              }
            }
            //Time complexity: amortized O(1), space complexity: O(N)
            class MyQueue {
              Deque<Integer> out, in;
              public MyQueue() {
                  in = new ArrayDeque<>();
                  out = new ArrayDeque<>();
              }
              
              public void push(int x) {
                  in.addLast(x);
              }
              
              public int pop() {
                  if (out.isEmpty()) {
                      while (!in.isEmpty()) {
                        out.addLast(in.pollLast());
                      }
                  }
                  return out.pollLast();
              }
              
              public int peek() {
                  if (out.isEmpty()) {
                      while (!in.isEmpty()) {
                        out.addLast(in.pollLast());
                      }
                  }
                  return out.peekLast();
              }
              
              public boolean empty() {
                  return out.isEmpty() && in.isEmpty();
              }
            } 
            ```

         * Max/Min Stack: use another stack to trace the max/min value



2. Elementary Sorts

    * Sorting Introduction

        * Callbacks

            * Goal: sort any type of data
            * Client passes array of objects to `sort()` function, the sort function calls back object's `compareTo()` method as needed

        * Implementing callbacks

            * Java: interface

                ```java
                public interface Comparable<Item> {
                  public int compareTo(Item that);
                }
                // v.compareTo(w) returns negative integer, zero, or positive if v < w, v == w, or v > w
                // if v and w are not comparable, throws an Exception
                ```

            * C: function pointers

            * C++: class-type pointers

            * Python, Perl, Javascript: first-class functions

        * Total order:

            * Antisymmetry: if $v \leq w$ and $w \leq v$, then $v = w$
            * Transitivity: if $v \leq w$ and $w \leq x$, then $v \leq x$
            * Totality: either $v \leq w$ or $w \leq v$ or both

    * Selection Sort

        * Idea to get an ascending array

            * In iteration i, find index `min` of smallest remaining entry
            * Swap `array[i]` and `array[min]`
            * The part on the left side of i is an invariant, that is, the entries are sorted

        * Implementation

            ```java
            public class SelectionSort {
              public static void sort(Comparable[] items) {
                int N = items.length;
                for (int i = 0; i < N; i++) {
                  int min = i;
                  for (int j = i + 1; j < N; j++) {
                    if (items[j].compareTo(items[min]) < 0) {
                      min = j;
                    }
                  }
                  Comparable item = items[min];
                  items[min] = items[i];
                  items[i] = item;
                }
              }
            }
            ```

        * For an array with size $N$, selection sort takes about $\frac{N^2}{2}$ comparisons and about $N$ exchanges. No matter whether the input is sorted or totally unsorted, the complexity of selection sort is always quadratic

    * Insertion Sort

        * Idea to get an ascending array

            * In iteration i, swap `array[i]` with each larger entry to its left
            * The part on the left side of i is an invariant, that is, the entries are sorted 

        * Implementation

            ```java
            public class InsertionSort {
              public static void sort(Comparable[] items) {
                int N = items.length;
                for (int i = 0; i < N; i++) {
                  for (int j = i; j > 0 && items[j].compareTo(items[j - 1]) < 0; j--) {  
                      Comparable item = items[j];
                      items[j] = items[i];
                      items[i] = item;
                  }
                }
              }
            }
            ```

        * For a randomly-ordered array with $N$ distinct keys, insertion sort uses about $\frac{N^2}{4}$ comparisons and about $\frac{N^2}{4}$ exchanges on average. 

        * For an ordered array, insertion sort makes $N - 1$ comparisons and $0$ exchanges

        * For an reversely-ordered array, insertion sort makes about $\frac{N^2}{2}$ comparisons and about $\frac{N^2}{4}$ exchanges

    * Shell Sort

        * Idea to get an ascending array

            * Move entries more than one position at a time by h-sorting the array
            * Use insertion sort in h-sorting: for large h, the subarray is small, so sorting can be fast, for small h, the array is nearly sorted, to insertion sort is fast
            * A m-sorted array is still m-sorted after being n-sorted
            * The general chosen h is $h = 3 * h + 1$

        * Implementation

            ```java
            public class ShellSort {
              public static void sort(Comparable[] items) {
                int N = items.length;
                int h = 1;
                while (h < N /3) {
                  h = 3 * h + 1;
                }
                while (h >= 1) {
                  for (int i = h; i < N; i++) {
                    for (int j = i; j >= h && items[j].compareTo(items[j - h]) < 0; j -= h) {
                      Comparable item = items[j];
                      items[j] = items[j - h];
                      items[j - h] = item;
                    }
                    h = h / 3;
                  }
                }
              }
            }
            ```

        * The worst-case number of comparisons used by shell sort with $3 * h + 1$ increments is $O(N^{\frac{3}{2}})$

        * Shell sort is faster than selection sort and insertion sort, the larger the array is, the better the performance is

    * Shuffling Sort

        * Initial idea to get a randomized shuffle of an array

            * Generate a random real number for each array entry
            * Sort the array

        * Shuffle sort produces a uniformly random permutation of the input array, provided no duplicate values

        * Knuth shuffle in linear time

            * In iteration i, pick integer r between 0 and i uniformly at random

            * Swap `array[r]` and `array[i]`

            * The part on the left side of i is an invariant, that is, the entries are uniformly randomized

            * Implementation

                ```java
                public class ShuffleSort {
                  public static void sort(Comparable[] items) {
                    int N = items.length;
                    for (int i = 0; i < N; i++) {
                      int r = StdRandom.uniform(i + 1); // Uniform random number in [0, i]
                      Comparable item = items[r];
                      items[r] = items[i];
                      items[i] = item;
                    }
                  }
                }
                ```

    * Convex Hull

        * Some equivalent definitions of a convex hull
            * The smallest perimeter fence enclosing a set of points
            * The smallest convex set containing all the points
            * The smallest area convex polygon enclosing the points
            * The convex polygon enclosing the points, whose vertices are points in set 
        * Output: sequence of vertices in counterclockwise order
        * Fact
            * Can traverse the convex hull by making only counterclockwise turns
            * The vertices of convex hull appear in increasing order of polar angle with respect to point p with lowest y-coordinate
        * Graham Scan
            * Choose point p with smallest y-coordinate
            * Sort points by polar angle with p
            * Consider points in order, discard unless it creates a counterclockwise turn



# Week 3

1. Merge Sort

    * Introduction 

        * Two classic sorting algorithm

            * Merge sort: java sort for objects, Perl, C stable sort, Python stable sort, Firefox javascript, ...
            * Quick sort: java sort for primitive types, C qsort, Unix, C++, Python, Matlab, Chrome javascript, ...

        * Idea is divide and conquer

            * Divide the array into two halves
            * Recursively sort each half
            * Merge the two halves using an auxiliary array

        * Implementation

            ```java
            public class MergeSort {
              private static Comparable[] aux;
              public static void sort(Comparable[] items) {
                aux = new Comparable[a.length];     // Do not create aux in the private recursive sorts
                sort(items, 0, items.length - 1);
              }
            
              private static void sort(Comparable[] items, int lo, int hi) {
                if (hi <= lo) {
                  return;
                }
            
                int mid = lo + (hi - lo) / 2;
                sort(items, lo, mid);
                sort(items, mid + 1, hi);
                merge(items, lo, mid, hi);
              }
            
              private static void merge(Comparable[] items, int lo, int mid, int hi) {
                // assert isSorted(items, lo, mid);
                // assert isSorted(items, mid + 1, hi);
            
                for (int i = 0; i <= hi; i++) {
                  aux[i] = items[i];
                }
            
                int i = lo; j = mid + 1;
                for(int k = lo; k <= hi; k++) {
                  if (i > mid) {
                    items[k] = aux[j++];              // items[lo...mid] has already been copied
                  } else if (j > hi) {
                    items[k] = items[i++];              // items[mid + 1...hi] has already been copied
                  } else if (aux[j].compareTo(aux[i]) < 0) {              
                    items[k] = aux[j++];                  // For stability reason, we cannot revert this else if statement and the else statement
                  } else {
                    items[k] = aux[i++];
                  }
                }
                
                // assert isSorted(items, lo, hi)
              }
            } 
            ```

        * Facts

            * Merge sort uses at most $NlogN$ comparisons and $6NlogN$ array access to sort any array of size $N$. If $D(N)$ satisfies $D(N) = 2D(\frac{N}{2}) + N$ for $N > 1$, with $D(1) = 0$, then $D(N) = N logN$
            * Merge sort uses extra space proportional to $N$
            * A sorting algorithm is in-place if it uses $\leq clogN$ extra memory. Ex. insertion sort, selection sort, shell sort. Merge sort can be in-place, [Kronrod, 1969]

        * Improvements

            * Merge sort is recursive, so it will be slow for small subarrays, we can use insertion sort for small subarrays

                ```java
                private int FUTOFF = 7;       // you can set your own cutoff to decide how small is an small array
                private static void sort(Comparable[] items, int lo, int hi) {
                  if (hi <= lo + CUTOFF - 1) {
                    InsertionSort.sort(items, lo, hi);
                    return;
                  }
                
                  int mid = lo + (hi - lo) / 2;
                  sort(items, aux, lo, mid);
                  sort(items, aux, mid + 1, hi);
                  merge(items, aux, lo, mid, hi);
                }
                ```

            * If `items[mid] <= items[mid + 1]`, then don't need to merge, the whole array is already sorted

                ```java
                private static void sort(Comparable[] items, int lo, int hi) {
                  if (hi <= lo + CUTOFF - 1) {
                    InsertionSort.sort(items, lo, hi);
                    return;
                  }
                
                  int mid = lo + (hi - lo) / 2;
                  sort(items, aux, lo, mid);
                  sort(items, aux, mid + 1, hi);
                  if (items[mid].compareTo(items[mid + 1]) < 0) {
                    return;
                  }
                  merge(items, aut, lo, mid, hi);
                }
                ```

            * Eliminate the copy to the auxiliary array

                ```java
                private static void merge(Comparable[] items, Comparable[] aux, int lo, int mid, int hi) {
                
                  int i = lo; 
                  int j = mid + 1;
                  for(int k = lo; k <= hi; k++) {
                    if (i > mid) {+
                      aux[k] = items[j++]              
                    } else if (j > hi) {
                      aux[k] = items[i++]              
                    } else if (items[j].compareTo(items[i]) < 0) {
                      aux[k] = items[j++]
                    } else {
                      aux[k] = items[i++];
                    }
                  }
                } 
                
                private static void sort(Comparable[] items, Comparable[] aux, int lo, int hi) {
                  if (hi <= lo + CUTOFF - 1) {
                    InsertionSort.sort(items, lo, hi);
                    return;
                  }
                
                  int mid = lo + (high - lo) / 2;
                  sort(aux, items, lo, mid);
                  sort(aux, items, mid + 1, hi);
                  if (aux[mid].compareTo(aux[mid + 1]) < 0) {
                    return;
                  } 
                  merge(aux, items, lo, mid, hi);
                }
                ```

                

    * Bottom-up Merge Sort

        * No recursion needed plan

            * Pass through array, merging subarrays of size 1
            * Repeat for subarrays of size 2, 4, 6, 8, 16, ..., array.length

        * Implementation

            ```java
            public class MergeSort {
              private static Comparable[] aux;
            
            
              private static void merge(Comparable[] items, int lo, int mid, int hi) {
                int i = lo; 
                int j = mid + 1;
                for (int k = lo; k <= hi; k++) {
                  if (i > mid) {
                    items[k] = aux[j++];
                  } else if (j > hi) {
                    items[k] = aux[i++];
                  } else if (aux[j].compareTo(aux[i]) < 0) {
                    items[k] = aux[j++];
                  } else {
                    items[k] = aux[i++];
                  }
                }
              }
            
              public static void sort(Comparable[] a) {
                int N = items.length;
                aux = new Comparable[N];
                for(int sz = 1; sz < N; sz *= 2) {
                  for(int lo = 0; lo < N - sz; lo += sz * 2) {
                    merge(items, lo, lo + sz - 1, Math.min(lo + sz * 2 - 1, N -1));
                  }
                }
              }
            }
            ```

    * Sorting Complexity

        * Complexity of sorting
            * Computational complexity: framework to study efficiency of algorithms for solving a particular problem X
            * Model of computation: allowable operations
            * Cost model: operation counts
            * Upper bound: cost guarantee provided by some algorithm for X
            * Lower bound: proven limit on cost guarantee of all algorithms for X
            * Optimal algorithm: algorithm with best possible cost guarantee for X
            * Sorting Ex.
                * Model of computation: decision tree
                    * The height of the tree is the worst-case number of comparisons
                    * There is at least one leaf for each possible ordering
                * Cost model: number of comparisons
                * Upper bound: about $NlogN$ from merge sort
                * Lower bound: Any compare-based sorting algorithm must at least $lg(N!) \sim NlgN (Stirling's formula)$ comparisons in the worst-case
                * Optimal algorithm: merge sort is an optimal algorithm
        * Merge sort is optimal with respect to number of comparisons but it's not optimal with respect to space usage
        * Complexity results in context: if more information is provided in terms of the initial order of input, the distribution of key values, or the representation of the keys, lower bound may not hold
            * Partially-ordered array: insertion sort requires only $N - 1$ comparisons if input is sorted
            * Duplicate keys: three-way quick sort can be used
            * Digital properties of keys: use digit/character comparisons instead of key comparisons for number and strings, radix sort is used here

    * Comparators

        * Comparable interface: sort using a type's natural order

        * Comparator interface: sort using an alternate order, the order must be a total order

            ```java
             public interface Comparator <Key> {
               int compare(Key v, Key w);
               // If v is less than w, then compare(v, w) < 0ssss
               // If v equals to w, then compare(v, w) == 0
               // If v is larger than w, then compare(v, w) > 0
            }
            ```

        * Usage in `Arrays.sort()`

            * Create Comparator instance, the instance is an object of a class that implements the Comparator interface
            * Pass as second argument to `Arrays.sort()`

    * Stability

        * A stable sort preserves the relative order of items with equal keys. If in the original array we have `array[i] == array[j], and i < j`, then in the sorted array, `array[i]` should still appear ahead of `array[j]`
        * Insertion sort and merge sort are stable, selection sort and shell sort are not stable. You need to check how the `compareTo()` method or the `compare()` method is implemented  
        * Insertion sort: we check `j > 0 && items[j].compareTo(items[j - 1]) < 0`, so no equal items never move past each other
        * Selection sort: there may be a long-distance exchange that might move an item past some equal item, here we `swap(array[i], array[min]), where min is in {i + 1, ..., array.length - 1}`
        * Merge sort: it's stable as long as merge operation is stable, and in the merge function, we have 

         ```java
         else if (aux[j].compareTo(aux[i]) < 0) {              
           items[k] = aux[j++];                  
         } else {
           items[k] = aux[i++];
         }
         ```

         hence if there is an equal situations(`aux[j] == aux[i]`), we use the left part `aux[i]` so that we will not change the relative order of `aux[i]` and `aux[j]`, therefore it's stable

        * Shell sort: there is a long-distance exchange in shell sort so it's not stable, we `swap(array[j], array[j - h])`, there may be values that are equal to `array[j - h]` in the interval of `[j - h + 1, , j - 1]`, so it's not stable
        * Question:
            ![](/images/Algorithm/Stability.png)
            * We only need to care about situation where $x_0 = x_1, and \; y_0 = y_1$, if these points can be put near each other, then it's fine
            * Choice A: A is fine
            * Choice B: quick sort x, so x is in order, merge sort y will not change the sorted x since it's stable
            * Choice C: merge sort x, so x is in order, quick sort y may change the sorted x since it's unstable
            * Choice D: merge sort is stable 

    * Problems:

        * Suppose that the subarray `a[0]` to `a[n - 1]` is sorted, and the subarray `a[n]` to `a[2 * n - 1]` is sorted. How can you merge the two subarrays so that array a is sorted using an auxiliary array of length n instead of 2n?

            ```java
              public void mergeWithShortAux(Comparable[] a, Comparable[] aux) {
                int N = a.length / 2;
                for (int i = 0; i < N; i++) {
                  aux[i] = a[i];
                }
                int i = 0; 
                int j = N;
                for (int k = 0; i < a.length; k++) {
                  if (i > N) {
                    a[k++] = a[j++];
                  } else if (j > a.length) {
                    a[k++] = a[i++];
                  } else if(a[j].CompareTo(aux[i]) < 0) {
                    a[k++] = a[j++];
                  } else {
                    a[k++] = aux[i++];
                  }
                }
              }
            ```

        * An inversion in an array `a[]` is a pair of entries `a[i]` and `a[j]` such that `i < j` but `a[i] > a[j]`. Given an array, design a linearithmic algorithm to count the number of inversions.

            ```java
            //Obviously there is an O(n^2) solution, for i in range(len(a)), if a[j] > a[i] for j in range(i + 1, len(a)), then count ++
            // The following one is from https://www.youtube.com/watch?v=Vj5IOD7A6f8
            //Idea inversion in a = recursive count (left half) + recursive count (right half) + merge count(left + right)
            private int merge(Comparable[] a, Comparable[] aux, int lo, int mid, int hi) {
              for (int k = lo; k <= hi; k++) {
                aux[k] = a[k];
              }
            
              int i = lo, j = mid + 1, count = 0;
              for (int k = lo; k <= hi; k++) {
                  if (i > mid)
                      a[k] = aux[j++];
                  else if (j > hi)
                      a[k] = aux[i++];
                  else if (aux[i].compareTo(aux[j]) > 0) { 
                      a[k] = aux[j++];
                      count += mid - i + 1;
                  } else { 
                      a[k] = aux[i++];
                  }
              }
              return count;
            }
            
            private int sort(Comparable[] a, Comparable[] aux, int lo, int hi) {
                
                if (lo >= hi) return 0;
                int mid = lo + (hi - lo)/2;
                int count1 = sort(a, aux, lo, mid);
                int count2 = sort(a, aux, mid + 1, hi);
                int count3 = merge(a, aux, lo, mid, hi);
                return count1 + count2 + count3;
            }
            
            public int inverseionCount(Comparable[] a) {
                Comparable[] aux = new Comparable[a.length];
                return sort(a, aux, 0, a.length - 1);
            }
            ```

1. Quick Sort

    * Introduction

        * Ieda

            * Randomly shuffle the array
            * Partition the array, so that for some j, `array[i]` is in place, no larger entry to the left of j, no smaller entry to the right of j
            * Sort each piece recursively

        * Partition

            * Phase I: repeat until i and j pointers cross
                * Scan i from left to right so long as `array[i] < array[lo]`
                * Scan j from right to left so long as `array[j] > array[lo]`
                * Exchange `array[i]` with `array[j]`
            * Phase II: when pointers cross.
                * Exchange `array[j]` with `array[lo]`

        * Java implementation

            ```java 
            public class QuickSort{
                private static int partition(Comparable[] items, int lo, int hi) {
                  int i = lo, j = hi + 1;
                  while (true) {
                    while (items[++i].compareTo(items[lo]) < 0) {
                      if (i == hi) {
                        break;
                      }
                    }
            
                    while (items[--j].compareTo(items[lo]) > 0) {
                      if (j == lo) {
                        break;
                      }
                    }
                    
                    if (i >= j) {
                      break;
                    }
                  }
                  Comparable item = items[j];
                    items[j] = items[lo];
                    items[lo] = item;
                    return j;
                }
            
                private static void sort(Comparable[] items, int lo, int hi) {
                  if (hi <= lo) {
                    return;
                  }
                  int j = partition(items, lo, hi);
                  sort(items, lo, j - 1);
                  sort(items, j + 1, hi);
                }
            
                public static void sort(Comparable[] items) {
                  StdRandom.shuffle(items);
                  sort(items, 0, items.length - 1);
                }
            }
            ```

        * Implementation Comments

            * The partition in quick sort is in-place, compared with merge sort, quick sort does not take any extra space
            * Testing whether the pointers cross is a bit trickier than it might seem, especially when there are duplicate keys
            * Testing `j ==lo` is redundant, but testing `i == hi` is not
            * Random shuffling is needed for performance guarantee
            * When duplicates are present, it is better to stop on keys equal to the partitioning item's key

        * Performance analysis

            * On average, quick sort makes $\sim 2NlogN$ comparisons and $\sim \frac{1}{3}NlogN$ exchanges
            * In best cases, quick sort makes $\sim NlogN$ comparisons
            * In worst cases, quick sort makes $\sim \frac{1}{2} N^{2}$ comparisons, randomly shuffling can avoid these situations

        * Properties

            * Quick sort is an in-place sorting algorithm
            * Quick sort is not stable

    * Selection

        * Goal: given an array of $N$ items, find the $k^{th}$ largest

        * Applications: order statistics, find the top k

        * Quick-select

            * Partition array so that:
                * Entry `array[j]` is in place
                * No larger entry to the left of j
                * No larger entry to the right of j
            * Repeat in one subarray, depending on j; finished when `j == k`

        * Java implementation

            ```java
            public static Comparable select(Comparable[] items, int k) {
              StdRandom.shuffle(items);
              int lo = 0, hi = items.length - 1;
              while (hi > lo) {
                int j = partition(items, lo, hi);
                if (j < k) {
                  lo = j + 1;
                } else if (j > k) {
                  hi = j - 1;
                } else {
                  return items[k]
                }
              }
              return items[k];
            }
            
            ```

        * Quick-select takes linear time on average

    * Duplicate Keys

        * Merge sort with duplicate keys: always between $\frac{1}{2}NlogN$ and $NlogN$ comparisons

        * Quick sort with duplicate keys: algorithm goes quadratic unless partitioning stops on equal keys

        * Dijkstra three-way partitioning quick sort

            * Let v be partitioning item `items[lo]`

            * Scan i from left to right:

                * If `items[i] < v`: exchange `items[lt]` with `items[i]`, increase both lt and i
                * If `items[i] > v`: exchange `items[gt]` with `items[i]`, decrease gt
                * If `items[i] == v`: increase i

            * Java implementation

                ```java
                public class Quick3Way {
                  private static void sort(Comparable[] items, int lo, int hi) {
                    if (hi <= lo) {
                      return;
                    }
                    int lt = lo, i = lo + 1, gt = hi;
                    Comparable v = items[lo];
                    while (i <= gt) {
                      int cmpt = items[i].compareTo(v);
                      if (cmp < 0) {
                        Comparable item = items[i];
                        items[i] = items[lt];
                        items[lt] = item;
                        lt++;
                        i++;
                      } else if (cmp > 0) {
                        Comparable item = items[i];
                        items[i] = items[gt];
                        items[gt] = item;
                        gt--;
                      } else {
                        i++;
                      }
                    }
                    sort(items, lo, lt - 1);
                    sort(items, gt + 1, hi);
                  }
                }
                ```

            * Quick sort with 3-way partitioning is entropy-optimal

            * Randomized quick sort with 3-way partitioning reduces running time from linearithmic to linear in board class of applications

    * System Sorts

        * `Arrays.sort()`
            * Uses tuned quick sort for primitive types; uses tuned merge sort for objects. For objects, merge sort is stable and guarantees $NlogN$ performance, for primitive types, quick sort uses less memory and is faster in generals. If you deliberately create primitive inputs without proper shuffling, the `Arrays.sort()` with quick sort can fail. 
            * Has different method for each primitive type
            * Has a method for data types that implements Comparable
            * Has a method that uses a Comparator
        * Which sort to choose
            * Internal sort:
                * Insertion sort, selection sort, bubble sort, shaker sort
                * Quick sort, merge sort, heap sort, sample sort, shell sort
                * Solitaire sort, red-black sort, splay sort, Yaroslavskiy sort, psort
            * External sort: Poly-phase merge sort, cascade-merge, oscillating sort
            * String/radix sorts: distribution, MSD, LSD, 3-way string quick sort
            * Parallel sort:
                * Bitonic sort, Batcher even-odd sort
                * Smooth sort, cube sort, column sort
                * GPU sort
        * Application attributes to consider: stable, parallel, deterministic, keys all distinct, multiple key types, linked list or arrays, large or small items, large or small items, is your array randomly ordered, need guarantee performance
        * Sort summary
            ![](/images/Algorithm/SortSummary.png)


    * Problems 
      * Selection in two sorted arrays. Given two sorted arrays `a[]` and `b[]`, of length $n_1$ and $n_2$ and an integer $0 \leq k \leq n_1 + n_2$, design an algorithm to find a key of rank k.
          * Solution 1: merge `a[]` and `b[]` into `c[]` where c[] is ordered with length $n_1 + n_2$, pick c[k - 1]
          * Solution 2: compare `a[i]` and `b[j]`, use `count` to track number of comparisons, if `count == k`, return `a[count]` or `b[count]`





# Week 4

1. Priority Queues
   * API and Elementary Implementations
     * API
        ```java
        public class MaxPQ<Key extends Comparable<Key>> {
          MaxPQ();                   // create an empty priority queue
          MaxPQ(Key[] items);        // create a priority queue with given keys
          void insert(Key v);        // insert a key into the priority queue
          Key delMax();              // return and remove the largest key
          boolean isEmpty();         // is the priority queue empty
          Key max();                 // return the largest key
          int size();                // number of entries in the priority queue
        }
        ```
     * Applications
       * Event-driven simulation: customers in a line
       * Data compression: huffman codes
       * Graph searching: Dijkstra's algorithm, Prim's algorithm
       * AI: A* search
       * Operating systems: load balancing, interrupt handling
       * Spam filtering: Bayesian spam filter
     * Unordered Array Implementation
        ```java
        public class UnorderedMaxPQ<Key extends Comparable<Key>> {
          private Key[] pq;
          private int N;

          public UnorderedMaxPQ(int capacity) {
            pq = (Key[])new Comparable[capacity];
          }

          public boolean isEmpty() {
            return N == 0;
          }

          public void insert(Key key) {
            pq[N++] = key;
          }

          public Key delMax() {
            int max = 0;
            for (int i = 1; i < N; i++) {
              if (pq[i].compareTo(pq[max]) > 0) {
                max = i;
              }
            }
            Key item = pq[max];
            pq[max] = pq[N - 1];
            pq[N - 1] = item;
            return pq[--N];
          }

          public Key max() {
            int max = 0;
            for (int i = 1; i < N; i++) {
              if (pq[i].compareTo(pq[max]) > 0) {
                max = i;
              }
            }
            return pq[max];
          }
        }
        ```
     * Challenge: implement all operations efficiently

        | Implementation  | insert | del max |   max  |
        |       :---:     |  :---: |  :---:  |  :---: |
        | unordered array |   1    |   $N$   |   $N$  |
        |  ordered array  |  $N$   |    1    |    1   |
        |       goal      |   1    |    1    |    1   |
   
   * Binary Heaps
     * Binary tree: empty or node with links to left and right binary trees
     * Complete tree: perfectly balanced, except for bottom level
       
       * Height of complete tree with $N$ nodes is $\left[ logN \right]$
     * Binary heap: array representation of a heap-ordered complete binary tree
       * Heap-ordered binary tree
         * A binary tree is heap-ordered if the key in each node is larger than (or equal to) the keys in that nodes two children (if any)
         * The largest key in a heap-ordered binary tree is found at the root
       * Array representation
         * Indices start at 1
         * Take nodes in level order
         * No explicit links needed
       * Can use array indices to move through tree
         * Parent of node at `k` is at `k/2`
         * Children of node `k` are at `2 * k` and `2 * k + 1`
       * Promotion in a heap: child's key becomes than larger than its parent's key
         * Exchange key in a child with key in parent
         * Repeat until heap order restored
            ```java
            private void swim(int k) {
              while (k > 1 && pq[k /2].compareTo(pq[k]) < 0) {
                Key item = pq[k / 2];
                pq[k / 2] = pq[k];
                pq[k] = item;
                k = k / 2;
              }
            }
            ```
       * Insertion in a heap
         * Add node at end, then make promotion if needed
         * At most $1 + logN$ comparisons
            ```java
            public void insert(Key value) {
              pq[++N] = value;
              swim(N);
            }
            ```
       * Demotion in a heap: parent's key becomes smaller than one(or both) of its children's
         * Exchange key in parent with key in larger child
         * Repeat until heap order restored
            ```java
            private void sink(int k) {
              while (2 * k <= N) {
                int j = 2 * k;
                if (j < N && pq[j].compareTo(pq[j + 1]) < 0) {
                  j++;
                }
                if (pq[k].compareTo(pq[j]) > 0) {
                  break;
                }
                Key item = pq[k];
                pq[k] = pq[j];
                pq[j] = item;
                k = j;
              }
            }
            ```
       * Delete the maximum in a heap
         * Exchange root with node at the end(the smallest one), and then sink it down
         * At most $2logN$ comparisons, for each level, 2 comparisons are needed, one to find the bigger child, one to determine whether a demotion is needed
            ```java
            public Key delMax() {
              Key max = pq[1];
              pq[1] = pq[N];
              pq[N--] = item;
              sink(1);
              pq[N + 1] = null;
              return max;
            }
            ```
     * Binary Heal Java Implementation
        ```java
        public class MaxPQ<Key extends Comparable<Key>> {
          private Keyp[] pq;
     private int N;
     
          public MaxPQ(int capacity) {
            pq = (Key[]) new Comparable[capacity + 1];
     }
     
          public boolean isEmpty() {
            return N == 0;
     }
     
          private void swim(int k) {
              while (k > 1 && pq[k /2].compareTo(pq[k]) < 0) {
                Key item = pq[k / 2];
                pq[k / 2] = pq[k];
                pq[k] = item;
                k = k / 2;
              }
       }
     
          public void insert(Key key) {
            pq[N++] = key;
            swim(N);
     }
     
          private void sink(int k) {
              while (2 * k <= N) {
                int j = 2 * k;
                if (j < N && pq[j].compareTo(pq[j + 1]) < 0) {
                  j++;
                }
                if (pq[k].compareTo(pq[j]) > 0) {
                  break;
                }
                Key item = pq[k];
                pq[k] = pq[j];
                pq[j] = item;
                k = j;
              }
       }
     
          public Key delMax() {
            Key max = pq[1];
            pq[1] = pq[N];
            pq[N--] = max;
            sink(1);
            pq[N + 1] = null;
            return max;
          }
        }
        ```
     * Implementation considerations
       * Immutability of keys
         * Immutable data type: cannot change the data type once created
         * Immutable types: String, Integer, Double, Color, Vector, ...
         * Mutable: StringBuilder, Stack, Counter, Java Array, ...
       * Underflow and overflow: remove when empty, add when full
       * Minimum-oriented priority queue
       * Other operations: remove an arbitrary item, change the priority of an item
     * Priority queue implementation cost summary
       
        | Implementation  |   insert  |     del max     |  max  |
        |      :---:      |   :---:   |      :---:      | :---: |
        | unordered array |     1     |       $N$       |  $N$  |
        |  ordered array  |   $N$     |        1        |   1   |
        |   binary heap   |  $logN$   |      $logN$     |   1   |
        |   d-ray heap    | $log_d N$ |  $d \; log_d N$ |   1   |
        |    Fibonacci    |     1     |      $logN$     |   1   |
   |   impossible    |     1     |        1        |   1   |
     
   * Heap Sort
     * Idea
       * Create max-heap with all N keys using bottom-up method, entries are indexed 1 to N
       * Repeatedly remove the maximum key
     * First pass: build heap using bottom-up method
        ```java
        private static void sink(Comparable[] items, int k, int n) {
          while (2 * k <= n) {
            int j = 2 * k;
            if (j < n && pq[j].compareTo(pq[j + 1]) < 0) {
              j++;
            }
            Comparable item = pq[k];
            pq[k] = pq[j];
            pq[j] = item;
            k = j;
          }
        }
        // heapify phase, we omit leaves, only consider root of sub-heap
        for (int k = n /2; k >= 1; k--) {
          sink(pq, k, n);
        }
        ```
     * Second pass: remove the maximum, once at a time, leave it in array instead of nulling out
       ```java
        while (N > 1) {
          Comparable item = pq[1];
          pq[1] = pq[N];
          pq[N] = item;
          sink(pq, 1, --N);
        }
       ```
     * Heap Sort Java Implementation
        ```java
        public class HeapSort {
          public static void sort(Comparable[] pq) {
            int N = pq.length;
            for (int k = N / 2; k >= 1; k--) {
              sink(pq, k, N);
            }
            while (N > 1) {
              Comparable item = pq[1];
              pq[1] = pq[N];
              pq[N] = item;
              sink(pq, 1, --N);
            }
       }
   
          private static void sink(Comparable[] items, int k, int n) {
            while (2 * k <= n) {
              int j = 2 * k;
              if (j < n && pq[j].compareTo(pq[j + 1]) < 0) {
                j++;
              }
              Comparable item = pq[k];
              pq[k] = pq[j];
              pq[j] = item;
              k = j;
            }
          }
        }
        ```
     * Proposition
       * Heap construction uses $\leq 2N$ comparisons and $\leq N$ exchanges
       * Heap sort uses $\leq 2 N log N + 2 N$ comparisons and $\leq N log N + N$ exchanges
       * Heap sort is the a in-place sorting algorithm with $N log N$ worst case(given a heap already), merge sort uses linear extra space, quick sort uses quadratic time in worst case
       * Heap sort is optimal for both time and space, but inner loop is longer than quick sort, it makes poor use of cache memory and it's not stable
     * Sort Summary
    ![](/images/Algorithm/SortSummary2.png)
   
   * Event-Driven Simulation
     * Goal: simulate the motion of $N$ moving particles that behave according to the laws of elastic collision
     * Hard sphere model(billiard ball model)
       * N particles in motion, confined in the unit box.
       * Particle i has known position $(rx_i, ry_i)$, velocity $(vx_i, vy_i)$, mass $m_i$, and radius $\sigma_i$.
       * Particles interact via elastic collisions with each other and with the reflecting boundary.
       * No other forces are exerted. Thus, particles travel in straight lines at constant speed between collisions.
     * Simulation
       * Time-driven simulation: Discretize time into quanta of size $\Delta t$. Update the position of each particle after every dt units of time and check for overlaps. If there is an overlap, roll back the clock to the time of the collision, update the velocities of the colliding particles, and continue the simulation. For $N$ particles during $\Delta t$ time, the time needed for simulation is proportional to $\frac{N^{2}}{\Delta t}$
       * Event-driven simulation: We focus only on those times at which interesting events occur, that is to determine the ordered sequence of particle collisions. We maintain a priority queue of future events, ordered by time. At any given time, the priority queue contains all future collisions that would occur, assuming each particle moves in a straight line trajectory forever. As particles collide and change direction, we use lay strategy to remove the invalid events. The simulation process is:
         * Delete the impending event, i.e., the one with the minimum priority t.
         * If the event corresponds to an invalidated collision, discard it. The event is invalid if one of the particles has participated in a collision since the time the event was inserted onto the priority queue.
         * If the event corresponds to a physical collision between particles i and j:
           * Advance all particles to time t along a straight line trajectory.
           * Update the velocities of the two colliding particles i and j according to the laws of elastic collision.
           * Determine all future collisions that would occur involving either i or j, assuming all particles move in straight line trajectories from time t onwards. Insert these events onto the priority queue.
         * If the event corresponds to a physical collision between particles i and a wall, do the analogous thing for particle i.
     * Collision prediction
       * Between particle and wall. If a particle with velocity $(vx, vy)$ collides with a wall perpendicular to x-axis, then the new velocity is $(-vx, vy)$; if it collides with a wall perpendicular to the y-axis, then the new velocity is $(vx, -vy)$
    * Between two particles.
   
         At $t$, we have $(rx_i, ry_i)$, $(vx_i, vy_i)$, $\sigma_i$, $m_i$, $(rx_j, ry_j)$, $(vx_j, vy_j)$, $\sigma_j$, $m_j$. At $t + \Delta t$, we have $(rx_{i}^{'}, ry_{i}^{'})$, $(vx_{i}^{'}, vy_{i}^{'})$, $\sigma_{i}$, $m_i$, $(rx_{j}^{'}, ry_{j}^{'})$, $(vx_{j}^{'}, vy_{j}^{'})$, $\sigma_j$, $m_j$
      
   
           $$rx_{i}^{'} = rx_{i} + \Delta t vx_{i}, \; ry_{i}^{'} = ry_{i} + \Delta t vy_{i} $$
           $$rx_{j}^{'} = rx_{j} + \Delta t vx_{j}, \; ry_{j}^{'} = ry_{j} + \Delta t vy_{j} $$
           $$vx_{i}^{'} = vx_{i} + \frac{J_{x}}{m_{i}}, \; vy_{i}^{'} = vy_{i} + \frac{J_{y}}{m_{i}} $$
        $$vx_{j}^{'} = vx_{j} + \frac{J_{x}}{m_{j}}, \; vy_{j}^{'} = vy_{j} + \frac{J_{y}}{m_{j}} $$
   
           where $J = \frac{2m_{i}m_{j}(\Delta v \cdot \Delta r)}{\sigma (m_{i} + m_{j})}$, $J_{x} = \frac{J\Delta rx}{\sigma}$, $J_{y} = \frac{J \Delta ry}{\sigma}$, $\Delta v = (\Delta vx, \Delta vy)$, $\Delta r = (\Delta rx, \Delta ry)$, $\sigma = \sigma_{i} + \sigma_{y}$
     * Java Implementation: I don't put the code here because they are too long, you can find them in the book or on the webpage for the book
       * Particle Class
       * Event class
       * CollisionSystem
     * Proposition: for $N$ particles, it uses at most $N^{2}$ priority queue operations at initialization, it uses at most $N$ priority queue operations during collision


2. Elementary Symbol Tables
   * Symbol Table API
     * Key-value pair abstraction
       * Insert a value with specified key
       * Given a key, search for the corresponding value
       * Examples: dictionary, book index, file share, web search, complier, routing table, DNS, reverse DNS, file system
     * Associative array abstraction
        ```java
        public class ST<Key, Value> {
          ST()                                // create a symbol table
          void put(Key key, Value val)        // put key-value pair into the table,remove key from table if value is null
          Value get(Key key)                  // get  value paired with key, null if key is absent
          void delete(Key key)                // remove key and its value from table
          boolean contains(Key key)           // is there is a value paired with key
          boolean isEmpty()                   // is the table empty
          int size()                          // number of key-value pairs in the table
          Iterable<Key> keys()                // all the keys in the table
        }
        ```
     * Conventions
       * Values are not null
       * Method `get()` returns null if key is not present
       * Method `put()` overwrites old value with new value
          ```java
          public boolean contains(Key key) {
            return get(key) != null;
          }

          public void delete(Key key) {
            put(key, null);
          }
          ```
     * Keys and values
       * Value type: any generic type
       * Key type:
         * Assume keys are `Comparable`, use `compareTo()`
         * Assume keys are any generic type, use `equals()` to test equality
         * Assume keys are any generic type, use `equals()` to test equality, use `hashCode()` to scramble key, this is built-in in Java
       * Best practices: use immutable types for symbol table keys
     * Equality test
       * All Java classes inherit a method `equals()` from `Object` class
       * Java requirements, for any references x, y and z
         * Reflexive: `x.equals(x)` is true
         * Symmetric: `x.equals(y)` iff `y.equals(x)`
         * Transitive: if `x.equals(y)` and `y.equals(z)`, then `x.equals(z)`
         * Non-null: `x.equals(null)` is false
         * An equivalence relation: reflexive, symmetric, transitive
       * Implementation
         * Default implementation: `x == y` in `Object` class
         * Customized implementation: `Integer`, `Double`, `String`, `File`, `URL`
         * User-defined implementation: 
            ```java
            public boolean equals(Object that) {
              
              if (this == that) {
                return true;
              }

              if (that == null) {
                return false;
              }

              if (this.getClass() != that.getClass()) {
                return false;
              }

              Type typeThat = (Type) that;
              
              // compare instance fields here
              // for primitive fields, use ==
              // for object fields, use equals
              // for array fields, check equivalence for each entry in the array
            }
            ```
   * Elementary Implementations
     * Linked list with sequential search
       * Maintain a linked list of key-value pairs
       * Search: scan through all keys until find a match
       * Insert: scan through all keys until find a match, if not match add to the front
       * Performance: on average, search $\sim O(N)$, insert $\sim O(N)$ 
     * Ordered array with binary search
       * Maintain an ordered array of key-value pairs
       * Search: use binary search to find a key
       * Insert: need to shift all greater keys over
       * Performance: on average, search $\sim O(logN)$, insert $\sim O(N)$
   * Ordered Operations
     * Ordered symbol table API
        ```java
        public class ST<Key extends Comparable<Key>, Value> {
          ST()                                         // create an ordered symbol table
          void put(Key key, Value val)                 // put key-value pair into the table,remove key from table is value is null
          Value get(Key key)                           // value paired with key, null if key absent
          void delete(Key key)                         // remove key and its value from table
         boolean contains(Key key)                     // is there a value paired with key
          boolean isEmpty()                            // is the table empty
          int size()                                   // number of key-value pairs
          Key min()                                    // smallest key
          Key max()                                    // largest key
          Key floor(Key key)                           // largest key less than or equal to key
          Key ceiling(Key key)                         // smallest key greater than or equal to key
          int rank(Key key)                            // number or keys less than key
          Key select(int k)                            // key of rank k
          void deleteMin()                             // delete smallest key
          void deleteMax()                             // delete largest key
          int size(Key lo, Key hi)                     // number of keys in [lo..hi]
          Iterable<Key> keys(Key lo, Key hi)           // keys in [lo..hi] in sorted order
          Iterable<Key> keys()                         // all keys in the table in sorted order
        }
        ```
     * Performance 

        |                   | sequential search | binary search   |
        |       :---:       |       :---:       |      :---:      |
        |       search      |        $N$        |     $logN$      |
        |    insert/delete  |        $N$        |       $N$       |
        |      min/max      |        $N$        |       $1$       |
        |   floor/ceiling   |        $N$        |     $logN$      |
        |       rank        |        $N$        |     $logN$      |
        |      select       |        $N$        |       $1$       |
        | ordered iteration |       $NlogN$     |       $N$       |
   * Binary Search Trees
     * BST
       * Definition: A BST is a binary tree in symmetric order
       * It's either empty or two disjoint binary trees(left tree and right tree)
       * Symmetric order: each node has a key, and every node's key is larger than all keys in its left subtree and smaller than all keys in its right subtree
     * BST in Java
       * Definition: A BSt is a reference to a root Node
       * A Node is comprised of four fields: a Key, a Value, a reference to left and right subtree
          ```java
          public class BST<Key extends Comparable<Key>, Value> {
            private Node root;

            public class Node {
              private Key key;                    // sorted by key
              private Value val;                  // associated data
              private Node left, right;           // left and right subtree
              private int N;                      // number of nodes in the tree rooted with current node
              
              public Node(Key key, Value val, int size) {
                this.key = key;
                this.value = val;
                this.size = size;
              }
            }  
            public boolean isEmpty() {
              return size() == 0;
            }

            private int size(Node x) {
              if (x == null) {
                return 0;
              } else {
                return x.size;
              }
            }

            public int size() {
              return size(root);
            }
          
            // helper function, if key not in BST, put key-value pair into BST, else update the value associated with key
            private void put(Node x, key key, Value val) {
              if (x == null) {
                return new Node(key, val);
              }
              int cmp = keu.compareTo(x.key);
              if (cmp < 0) {
                x.left = put(x.left, key, val);
              } else if (cmp > 0) {
                x.right = put(x.right, key, val);
              } else if (cmp == 0) {
                x.val = val;
              }
              // the recursive calls in the above if-else statements update size with the corresponding node
              x.N = size(x.left) + size(x.right) + 1;
              return x;
            }

            public void put(Key key, Value val) {
              root = put(root, key, val);
            }

            // helper function, is key is not int BST, return null, else find the key recursively and return the value paired with key
            private Value get(Node x, Key key) {
              if (x == null) {
                return null;
              }
              int cmp = key.compareTo(x.key);
              if (cmp < 0) {
                return get(x.left, key);
              } else if(cmp > 0) {
                return get(x.right, key);
              } else {
                return x.val;
              }
            }

            public Value get(Key key) {
              return get(root, key);
            }

            private Node min(Node x) {
              if (x.left == null) {
                return x;
              } 
              return min(x.left);
            }

            public Key min() {
              return min(root).key;
            }

            private Node max(Node x) {
              if (x.right == null) {
                return x;
              }
              return max(x.right);
            }

            public Key max() {
              return max(root).key;
            }
            
            private Node floor(Node x, Key key) {
              if (x == null) {
                return null;
              }
              int cmp = key.compareTo(x.key);
              if (cmp == 0) {
                return x;
              } else if (cmp < 0) {
                return floor(x.left, key);
              }
              Node t = floor(x.right, key);
              if (t != null) {
                return t;
              } else {
                return x;
              }
            }

            public Key floor(Key key) {
              Node x = floor(root, key);
              if (x == null) {
                return null;
              }
              return x.key;
            }

            private Node ceiling(Node x, Key key) {
              if (x == null) {
                return null;
              }
              int cmp = key.compareTo(x.key);
              if (cmp == 0) {
                return x;
              } else if(cmp > 0) {
                return ceiling(x.right, key);
              }
              Node t = ceiling(x.left, key);
              if (t!= null) {
                return t;
              } else {
                return x;
              }
            }

            public Key floor(Key key) {
              Node x = ceiling(root, key);
              if (x == null) {
                return null;
              }
              return x.key;
            }

            private int rank(Node x, Key key) {
              if (x == null) {
                return 0;
              }
              int cmp = key.compareTo(x.key);
              if (cmp < 0) {
                return rank(x.left, key);
              } else if (cmp == 0) {
                return size(x.left);
              } else if (cmp > 0) {
                return rank(x.right, key) + size(x.left) + 1;
              }
            }

            public int rank(Key key) {
              return rank(root, key);
            }

            private Node deleteMin(Node x) {
              if (x.left == null) {
                return x.right;
              }
              x.left = deleteMin(x.left);
              x.N = size(x.left) + size(s.right) + 1;
              return x;
            }

            public void deleteMin() {
              root = deleteMin(root);
            }

            private Node deleteMax(Node x) {
              if (x.right == null) {
                return x.left;
              }
              x.right = deleteMax(x.right);
              x.N = size(x.left) + size(x.right) + 1;
              return x;
            }

            public void deleteMax() {
              root = deleteMax(root);
            }

            private Node(Node x, Key key) {
              if (x == null) {
                return null;
              }
              int cmp = key.compareTo(x.key);
              if (cmp < 0) {
                x.left = delete(x.left, key);
              } else if (cmp > 0) {
                x.right = delete(x.right, key);
              } else {
                if (x.left == null) {
                  return x.right;
                }
                if (x.right == null) {
                  return x.left;
                }
                Node temp = x;
                x = min(temp.right);
                x.right = deleteMin(temp.right);
                x.left = temp.left;
              }
              x.N = size(x.left) + size(x.right) + 1;
              return x;
            }

            public void delete(Key key) {
              root = delete(root, key);
            }

            private void inorder(Node x, Queue<Key> queue) {
              if (x == null) {
                return;
              }
              inorder(x.left, queue);
              queue.enqueue(x);
              inorder(x.right, queue);
            }

            public Iterable<Key> iterator() {
              Queue<Key> queue = new Queue<<Key>();
              inorder(root, queue);
              return queue;
            }
          }
          ```
       * Search: if less, go left; if greater, go right; if equal, search hit. If null but not hit, search miss
       * Insert: if less, go left; if greater, go right; if null, insert, if old value exists, update it to new value
       * Performance: on average, search $\sim O(logN) (\sim 1.39 logN)$, insert $\sim O(logN) (\sim 1.39 logN)$
   * Ordered Operations in BSTs
     * Operations
       * Minimum: smallest key in the table, the leftmost one
       * Maximum: largest key in the table, the rightmost one
       * Floor: largest key $\leq$ to a given key, if `key == root.key`, result is `root.key`, if `key < root.key`, result is in the left subtree, if `key > root.key`, result is in the right subtree
       * Ceiling: smallest key $\geq$ to a given key, if `key == root.key`, result is `root.key`, if `key > root.key`, result is in the left subtree, if `key < root.key`, result is in the right subtree
       * Rank: number of keys in the symbol table strictly less than
       * In-order traversal: maintain a queue, traverse left subtree, enqueue key, traverse right subtree
   * Deletion in BSTs
     * Lazy approach deletion
       * To remove a node with a give key:
         * Set its value to `null` 
         * Leave key in tree to guide searches
       * Cost: $\sim 2logN$ per insertion, search, and deletion
     * Delete the minimum:
       * To delete the minimum key:
         * Go left until finding a node with a null left link
         * Replace that node by its right link
         * Update subtree counts
     * Hibbard deletion: delete a node with key k: search for node t containing key k
         * Node t has 0 children: delete t by setting parent link to null
         * Node t has 1 child: delete t by replacing parent link to t with parent link to t's child
         * Node t has 2 children:
           * Find the minNode in the right subtree of t
           * Replace t with minNode
           * New t's right subtree is the right subtree of old t without minNode
           * New t's left subtree is the left subtree of old t
           * Update size of t
         * Performance:  on average, if deletion is allowed, then search, insertion and deletion all become $\sim \sqrt{N}$. But for Red-Black BSTs guarantee logarithm cost for all the three operations.
   * Problems
     * Auto-boxing 
        ```java
        double a = 0.0;
        double b = -0.0;
        System.out.println(a == b); // true
        System.out.println(Double.valueOf(a) == Double.valueOf(b)); // false

        double c = Double.NaN;
        double d = Double.NaN;

        System.out.println(c == d); // false
        System.out.println(Double.valueOf(c) == Double.valueOf(d)); // true
        ```
     * isBST
        ```java
        public boolean isBST() {
          return isBSTHelper(root, Integer.MIN_VALUE, Integer.MAX_VALUE);
        }

        private isBSTHelper(Node node, int minVal, int maxVal) {
          if (node == null) {
            return true;
          }

          if (node.key < minVal || node.key > maxVal) {
            return false;
          }

          return isBSTHelper(node.left, minVal - 1, maxVal) && isBSTHelper(node.right, minVal, maxVal - 1);
        }
        ```
     * Inorder traversal
       * Recursion approach
          ```java
          public void inorder(Node x) {
            if (x == null) {
              return;
            }
            inorder(x.left);
            System.out.println(x);
            inorder(x.right);
          }
          ```
       * Iteration with stack approach
          ```java
          public void inorder(Node root) {
            Deque<Node> stack = new ArrayDeque<Node>();
            Node current = root;
            while (!stack.isEmpty() || current != null) {
              if (node != null) {
                stack.addLast(node);
                node = node.left;
              } else {
                Node temp = stack.removeLast();
                System.out.println(node);
                node = node.right;
              }
            }
          }
          ```
       * Morris traversal: no recursion, no stack
         * Algorithm
           * Initialize the root as the current node `curr`
           * While `curr` is not null, check if `curr` has a left child
           * If `curr` does not have a left child, print `curr` and update it to `curr.right`
           * Else, make `curr` the right child of th rightmost node in `curr`'s left subtree(the rightmost node in `curr`'s left subtree is the inorder predecessor, we give it a link to `curr` if there is no link, and we remove the link if there has already been one)
           * Update `curr` to his left node
         * Video: <a href = "https://www.youtube.com/watch?v=wGXB9OWhPTg">https://www.youtube.com/watch?v=wGXB9OWhPTg</a>
            ```java
            public void inorder(Node root) { 
              Node curr, prev; 

              if (root == null) 
                return; 

              curr = root; 

              while (curr != null) { 
                if (curr.left == null) { 
                  System.out.print(curr.key); 
                  curr = curr.right; 
                } 
                else { 
                  // find the prev of curr 
                  prev = curr.left; 
                  // for the first time, prev.right !=null is termination, for the second time, prev.right != current is termination
                  while (prev.right != null && prev.right != curr) 
                    prev = prev.right; 

                  // first time, link prev to current
                  if (prev.right == null) { 
                    prev.right = curr; 
                    curr = curr.left; 
                  } 
                  // second time, remove the link
                  else { 
                    prev.right = null; 
                    System.out.print(curr.key); 
                    curr = curr.right; 
                  } 
                } 
              }
            } 
            ```