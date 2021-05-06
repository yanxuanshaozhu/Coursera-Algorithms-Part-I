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
    private final MinPQ<SearchNode> checkGoal;
    private final MinPQ<SearchNode> twinCheckGoal;

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


        checkGoal = new MinPQ<>();
        twinCheckGoal = new MinPQ<>();

        checkGoal.insert(initialNode);
        twinCheckGoal.insert(initialTwin);

        while (!checkGoal.min().board.isGoal() && !twinCheckGoal.min().board.isGoal()) {
            SearchNode node = checkGoal.delMin();
            SearchNode twin = twinCheckGoal.delMin();
            for (Board neighbor : node.board.neighbors()) {
                if (!neighbor.equals(node.board)) {
                    SearchNode newNode = new SearchNode(neighbor, node.moves + 1, node);
                    checkGoal.insert(newNode);
                }
            }
            for (Board twinNeighbor : twin.board.neighbors()) {
                if (!twinNeighbor.equals(twin.board)) {
                    SearchNode twiNewNode = new SearchNode(twinNeighbor, twin.moves + 1, twin);
                    twinCheckGoal.insert(twiNewNode);
                }
            }
        }
    }

    // is the initial board solvable? (see below)
    public boolean isSolvable() {
        return checkGoal.min().board.isGoal();
    }

    // min number of moves to solve initial board; -1 if unsolvable
    public int moves() {
        if (!isSolvable()) {
            return -1;
        }
        return checkGoal.min().moves;
    }

    // sequence of boards in a shortest solution; null if unsolvable
    // According to the autograder, the first in the solution should be the initial board,
    // so we use a stack here
    public Iterable<Board> solution() {
        if (!isSolvable()) {
            return null;
        }
        Stack<Board> boards = new Stack<>();
        SearchNode node = checkGoal.min();
        while (node.prev != null) {
            boards.push(node.board);
            node = node.prev;
        }
        boards.push(node.board);
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
