/* *****************************************************************************
 *  Name:  yanxuanshaozhu
 *  Date:  05/24/2021
 *  Description: coursera algorithm week 4 project
 **************************************************************************** */

import edu.princeton.cs.algs4.StdRandom;

import java.util.ArrayList;
import java.util.Arrays;

public class Board {
    private int[][] board;
    private final int size;

    // create a board from an n-by-n array of tiles,
    // where tiles[row][col] = tile at (row, col)
    public Board(int[][] tiles) {
        board = new int[tiles.length][tiles.length];
        for (int i = 0; i < tiles.length; i++) {
            for (int j = 0; j < tiles.length; j++) {
                this.board[i][j] = tiles[i][j];
            }
        }
        this.size = tiles.length;
    }

    // string representation of this board
    public String toString() {
        StringBuilder s = new StringBuilder();
        s.append(this.size + "\n");
        for (int i = 0; i < this.size; i++) {
            for (int j = 0; j < this.size; j++) {
                s.append(String.format("%2d ", this.board[i][j]));
            }
            s.append("\n");
        }
        return s.toString();
    }

    // board dimension n
    public int dimension() {
        return this.size;
    }

    // number of tiles out of place
    public int hamming() {
        int count = 0;
        for (int i = 0; i < this.size; i++) {
            for (int j = 0; j < this.size; j++) {
                if (i == this.size - 1 && j == this.size - 1) {
                    continue;
                } else if (this.board[i][j] != i * this.size + j + 1) {
                    count += 1;
                }

            }
        }
        return count;
    }

    // sum of Manhattan distances between tiles and goal
    public int manhattan() {
        int count = 0;
        for (int i = 0; i < this.size; i++) {
            for (int j = 0; j < this.size; j++) {
                if (this.board[i][j] != 0) {
                    int trueRow = (this.board[i][j] - 1) / this.size;
                    int trueCol = this.board[i][j] - 1 - trueRow * this.size;
                    int dist = Math.abs(j - trueCol) + Math.abs(i - trueRow);
                    count += dist;
                }
            }
        }
        return count;
    }

    // is this board the goal board?
    public boolean isGoal() {
        if (manhattan() == 0 && hamming() == 0) {
            return true;
        } else {
            return false;
        }
    }

    // does this board equal y?
    public boolean equals(Object y) {
        if (!(y instanceof Board)) {
            return false;
        } else {
            if (this.size != ((Board) y).size) {
                return false;
            } else {
                for (int i = 0; i < this.size; i++) {
                    for (int j = 0; j < this.size; j++) {
                        if (this.board[i][j] != ((Board) y).board[i][j]) {
                            return false;
                        }
                    }
                }
            }
        }
        return true;
    }

    // all neighboring boards
    public Iterable<Board> neighbors() {
        int[][] temp = new int[this.board.length][];
        for (int i = 0; i < this.board.length; i++) {
            temp[i] = Arrays.copyOf(this.board[i], this.board[i].length);
        }
        ArrayList<Board> boards = null;
        // two neighbours: 0 is at one of the 4 corners
        if (this.board[0][0] == 0) {
            boards = new ArrayList<>(2);
            temp[0][0] = temp[0][1];
            temp[0][1] = 0;
            Board board1 = new Board(temp);
            boards.add(board1);
            temp[0][1] = temp[0][0];
            temp[0][0] = temp[1][0];
            temp[1][0] = 0;
            Board board2 = new Board(temp);
            boards.add(board2);
            return boards;

        }
        if (this.board[0][this.size - 1] == 0) {
            boards = new ArrayList<>(2);
            temp[0][this.size - 1] = temp[0][this.size - 2];
            temp[0][this.size - 2] = 0;
            Board board1 = new Board(temp);
            boards.add(board1);
            temp[0][this.size - 2] = temp[0][this.size - 1];
            temp[0][this.size - 1] = temp[1][this.size - 1];
            temp[1][this.size - 1] = 0;
            Board board2 = new Board(temp);
            boards.add(board2);
            return boards;
        }
        if (this.board[this.size - 1][0] == 0) {
            boards = new ArrayList<>(2);
            temp[this.size - 1][0] = temp[this.size - 2][0];
            temp[this.size - 2][0] = 0;
            Board board1 = new Board(temp);
            boards.add(board1);
            temp[this.size - 2][0] = temp[this.size - 1][0];
            temp[this.size - 1][0] = temp[this.size - 1][1];
            temp[this.size - 1][1] = 0;
            Board board2 = new Board(temp);
            boards.add(board2);
            return boards;
        }
        if (this.board[this.size - 1][this.size - 1] == 0) {
            boards = new ArrayList<>(2);
            temp[this.size - 1][this.size - 1] = temp[this.size - 1][this.size - 2];
            temp[this.size - 1][this.size - 1] = 0;
            Board board1 = new Board(temp);
            boards.add(board1);
            temp[this.size - 1][this.size - 2] = temp[this.size - 1][this.size - 1];
            temp[this.size - 1][this.size - 1] = temp[this.size - 2][this.size - 1];
            temp[this.size - 2][this.size - 1] = 0;
            Board board2 = new Board(temp);
            boards.add(board2);
            return boards;
        }

        // three neighbours: 0 is on one of the 4 sides
        for (int i = 1; i < this.size - 1; i++) {
            boards = new ArrayList<>(3);
            if (this.board[i][0] == 0) {
                temp[i][0] = temp[i - 1][0];
                temp[i - 1][0] = 0;
                Board board1 = new Board(temp);
                boards.add(board1);
                temp[i - 1][0] = temp[i][0];
                temp[i][0] = temp[i + 1][0];
                temp[i + 1][0] = 0;
                Board board2 = new Board(temp);
                boards.add(board2);
                temp[i + 1][0] = temp[i][0];
                temp[i][0] = temp[i][1];
                temp[i][1] = 0;
                Board board3 = new Board(temp);
                boards.add(board3);
                return boards;
            }
            if (this.board[i][this.size - 1] == 0) {
                temp[i][this.size - 1] = temp[i - 1][this.size - 1];
                temp[i - 1][this.size - 1] = 0;
                Board board1 = new Board(temp);
                boards.add(board1);
                temp[i - 1][this.size - 1] = temp[i][this.size - 1];
                temp[i][this.size - 1] = temp[i + 1][this.size - 1];
                temp[i + 1][this.size - 1] = 0;
                Board board2 = new Board(temp);
                boards.add(board2);
                temp[i + 1][this.size - 1] = temp[i][this.size - 1];
                temp[i][this.size - 1] = temp[i][this.size - 2];
                temp[i][this.size - 2] = 0;
                Board board3 = new Board(temp);
                boards.add(board3);
                return boards;
            }
            if (this.board[0][i] == 0) {
                temp[0][i] = temp[0][i - 1];
                temp[0][i - 1] = 0;
                Board board1 = new Board(temp);
                boards.add(board1);
                temp[0][i - 1] = temp[0][i];
                temp[0][i] = temp[0][i + 1];
                temp[0][i + 1] = 0;
                Board board2 = new Board(temp);
                boards.add(board2);
                temp[0][i + 1] = temp[0][i];
                temp[0][i] = temp[1][i];
                temp[1][i] = 0;
                Board board3 = new Board(temp);
                boards.add(board3);
                return boards;
            }
            if (this.board[this.size - 1][i] == 0) {
                temp[this.size - 1][i] = temp[this.size - 1][i - 1];
                temp[this.size - 1][i - 1] = 0;
                Board board1 = new Board(temp);
                boards.add(board1);
                temp[this.size - 1][i - 1] = temp[this.size - 1][i];
                temp[this.size - 1][i] = temp[this.size - 1][i + 1];
                temp[this.size - 1][i + 1] = 0;
                Board board2 = new Board(temp);
                boards.add(board2);
                temp[this.size - 1][i + 1] = temp[this.size - 1][i];
                temp[this.size - 1][i] = temp[this.size - 2][i];
                temp[this.size - 2][i] = 0;
                Board board3 = new Board(temp);
                boards.add(board3);
                return boards;
            }

        }
        // 4 neighbours: 0 is at the interior of the board
        for (int i = 1; i < this.size - 1; i++) {
            for (int j = 1; j < this.size - 1; j++) {
                boards = new ArrayList<>(4);
                if (temp[i][j] == 0) {
                    temp[i][j] = temp[i - 1][j];
                    temp[i - 1][j] = 0;
                    Board board1 = new Board(temp);
                    boards.add(board1);
                    temp[i - 1][j] = temp[i][j];
                    temp[i][j] = temp[i + 1][j];
                    temp[i + 1][j] = 0;
                    Board board2 = new Board(temp);
                    boards.add(board2);
                    temp[i + 1][j] = temp[i][j];
                    temp[i][j] = temp[i][j + 1];
                    temp[i][j + 1] = 0;
                    Board board3 = new Board(temp);
                    boards.add(board3);
                    temp[i][j + 1] = temp[i][j];
                    temp[i][j] = temp[i][j - 1];
                    temp[i][j - 1] = 0;
                    Board board4 = new Board(temp);
                    boards.add(board4);
                    return boards;
                }
            }
        }
        return null;
    }

    // a board that is obtained by exchanging any pair of tiles
    public Board twin() {
        int n = StdRandom.uniform(this.size);
        int m = 0;
        if (n == 0) {
            m = StdRandom.uniform(this.size - 1);
            m += 1;
        } else {
            m = StdRandom.uniform(n);
        }
        int[][] temp = new int[this.board.length][];
        for (int i = 0; i < this.board.length; i++) {
            temp[i] = Arrays.copyOf(this.board[i], this.board[i].length);
        }
        int val = temp[n][n];
        temp[n][n] = temp[m][m];
        temp[m][m] = val;
        return new Board(temp);
    }

    // unit testing (not graded)
    // public static void main(String[] args) {
    //     int[][] nums = { { 1, 4, 2 }, { 3, 0, 5 }, { 6, 7, 8 } };
    //     Board board = new Board(nums);
    //     System.out.println("size: " + board.size);
    //     System.out.println("print: \n " + board);
    //     System.out.println("Hamming:" + board.hamming());
    //     System.out.println("Manhattan: " + board.manhattan());
    //     System.out.println("self: " + board);
    //     Iterable<Board> neighbors = board.neighbors();
    //     for (Board neighbor : neighbors) {
    //         System.out.println("neighbor:\n" + neighbor);
    //     }
    //     System.out.println("self: \n" + board);
    //     Board twin = board.twin();
    //     System.out.println("twin: " + twin);
    // }
}
