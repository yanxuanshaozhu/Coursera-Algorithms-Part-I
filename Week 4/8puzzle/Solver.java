/* *****************************************************************************
 *  Name:  yanxuanshaozhu
 *  Date:  05/24/2021
 *  Description: coursera algorithm week 4 project
 **************************************************************************** */

import edu.princeton.cs.algs4.In;
import edu.princeton.cs.algs4.MinPQ;
import edu.princeton.cs.algs4.StdOut;

import java.util.ArrayList;

public class Solver {
    private MinPQ<SearchNode> checkGoal;
    private MinPQ<SearchNode> twinCheckGoal;

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
            Iterable<Board> neighbors = node.board.neighbors();
            for (Board neighbor : neighbors) {
                SearchNode newNode = new SearchNode(neighbor, node.moves + 1, node);
                checkGoal.insert(newNode);
            }
            Iterable<Board> twinNeighbors = twin.board.neighbors();
            for (Board twinNeighbor : twinNeighbors) {
                SearchNode twiNewNode = new SearchNode(twinNeighbor, twin.moves + 1, twin);
                twinCheckGoal.insert(twiNewNode);
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
    public Iterable<Board> solution() {
        if (!isSolvable()) {
            return null;
        }
        ArrayList<Board> boards = new ArrayList<>();
        SearchNode node = checkGoal.min();
        while (node.prev != null) {
            boards.add(node.board);
            node = node.prev;
        }
        boards.add(node.board);
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
