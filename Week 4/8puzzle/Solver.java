/* *****************************************************************************
 *  Name:  yanxuanshaozhu
 *  Date:  05/24/2021
 *  Description: coursera algorithm week 4 project
 **************************************************************************** */

import edu.princeton.cs.algs4.MinPQ;

public class Solver {

    private class SearchNode implements  Comparable<SearchNode> {
        Board board;
        int moves;
        SearchNode prev;


        public int compareTo(SearchNode o) {
            int distThis = this.board.hamming() + this.moves;
            int distThan = o.board.dimension() + o.moves;
        }
    }

    // find a solution to the initial board (using the A* algorithm)
    public Solver(Board initial) {
        if (initial == null) {
            throw new IllegalArgumentException("null constructor argument");
        }

        SearchNode initialNode = new SearchNode();
        initialNode.board = initial;
        initialNode.moves = 0;
        initialNode.prev = null;

        MinPQ<SearchNode> pq = new MinPQ<>();
        pq.insert(initialNode);
        SearchNode node = pq.delMin();
        while (!node.board.isGoal()) {
            Iterable<Board> neighbors = node.board.neighbors();
            for (Board neighbor : neighbors) {
                SearchNode newNode = new SearchNode();
                newNode.board = neighbor;
                newNode.moves = node.moves + 1;
                newNode.prev = node;
                pq.insert(newNode);
            }
            node = pq.delMin();
        }


    }

    // is the initial board solvable? (see below)
    public boolean isSolvable() {
        return false;
    }

    // min number of moves to solve initial board; -1 if unsolvable
    public int moves() {
        return 0;
    }

    // sequence of boards in a shortest solution; null if unsolvable
    public Iterable<Board> solution() {
        return null;
    }

    // test client (see below)
    public static void main(String[] args) {

    }

}
