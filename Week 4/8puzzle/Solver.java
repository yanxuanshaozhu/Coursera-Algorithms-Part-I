/* *****************************************************************************
 *  Name:  yanxuanshaozhu
 *  Date:  05/05/2021
 *  Description: coursera algorithm week 4 project
 **************************************************************************** */

import edu.princeton.cs.algs4.In;
import edu.princeton.cs.algs4.MinPQ;
import edu.princeton.cs.algs4.Stack;
import edu.princeton.cs.algs4.StdOut;

public class Solver {
    private SearchNode lastNode;
    private boolean isSolvable;

    private class SearchNode implements Comparable<SearchNode> {
        Board board;
        int moves;
        SearchNode prev;
        int priority;

        public SearchNode(Board board, int moves, SearchNode prev) {
            this.board = board;
            this.moves = moves;
            this.prev = prev;
            this.priority = this.board.manhattan() + this.moves;
        }

        public int compareTo(SearchNode that) {
            return this.priority - that.priority;
        }
    }

    // find a solution to the initial board (using the A* algorithm)
    public Solver(Board initial) {
        if (initial == null) {
            throw new IllegalArgumentException("null constructor argument");
        }

        SearchNode initialNode = new SearchNode(initial, 0, null);
        SearchNode initialTwin = new SearchNode(initial.twin(), 0, null);


        MinPQ<SearchNode> checkGoal = new MinPQ<>();
        MinPQ<SearchNode> twinCheckGoal = new MinPQ<>();

        checkGoal.insert(initialNode);
        twinCheckGoal.insert(initialTwin);

        while (true) {
            SearchNode current = checkGoal.delMin();
            SearchNode currentTwin = twinCheckGoal.delMin();

            if (current.board.isGoal()) {
                this.lastNode = current;
                this.isSolvable = true;
                break;
            } else if (currentTwin.board.isGoal()) {
                this.isSolvable = false;
                break;
            } else {
                for (Board neighbor : current.board.neighbors()) {
                    if (current.prev == null) {
                        checkGoal.insert(new SearchNode(neighbor, current.moves + 1, current));
                    } else {
                        if (!current.prev.board.equals(neighbor)) {
                            checkGoal.insert(new SearchNode(neighbor, current.moves + 1, current));
                        }
                    }
                }

                for (Board neighbor : currentTwin.board.neighbors()) {
                    if (currentTwin.prev == null) {
                        twinCheckGoal.insert(new SearchNode(neighbor, currentTwin.moves + 1, currentTwin));
                    } else {
                        if (!currentTwin.prev.board.equals(neighbor)) {
                            twinCheckGoal.insert(new SearchNode(neighbor, currentTwin.moves + 1, currentTwin));
                        }
                    }
                }
            }
        }
    }

    // is the initial board solvable? (see below)
    public boolean isSolvable() {
        return this.isSolvable;
    }

    // min number of moves to solve initial board; -1 if unsolvable
    public int moves() {
        if (!this.isSolvable) {
            return -1;
        }
        return this.lastNode.moves;
    }

    // sequence of boards in a shortest solution; null if unsolvable
    // According to the autograder, the first in the solution should be the initial board,
    // so we use a stack here
    public Iterable<Board> solution() {
        if (!isSolvable()) {
            return null;
        }
        Stack<Board> boards = new Stack<>();
        SearchNode node = this.lastNode;
        while (node != null) {
            boards.push(node.board);
            node = node.prev;
        }
        return boards;
    }

    // test client (see below)
    public static void main(String[] args) {
        // create initial board from file
        In in = new In(args[0]);
        int n = in.readInt();
        int[][] tiles = new int[n][n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                tiles[i][j] = in.readInt();
        Board initial = new Board(tiles);
        // solve the puzzle
        Solver solver = new Solver(initial);

        // print solution to standard output
        if (!solver.isSolvable())
            StdOut.println("No solution possible");
        else {
            StdOut.println("Minimum number of moves = " + solver.moves());
            for (Board board : solver.solution())
                StdOut.println(board);
        }
    }
}
