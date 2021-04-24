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
              for (int j = i; j > 0; j--) {
                if (items[j].compareTo(items[i]) > 0) {
                  Comparable item = items[j];
                  items[j] = items[i];
                  items[i] = item;
                }
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

