/* *****************************************************************************
 *  Name:              yanxuanshaozhu
 *  Coursera User ID:  yanxuanshaozhu
 *  Last modified:     04/18/2021
 **************************************************************************** */

import edu.princeton.cs.algs4.WeightedQuickUnionUF;

public class Percolation {
    private final WeightedQuickUnionUF uf;
    private boolean[] grid;
    private int openNum;
    private final int size;

    /**
     * Constructor of the percolation system
     *
     * @param n side length of the grid
     */
    public Percolation(int n) {
        if (n <= 0) {
            throw new IllegalArgumentException("Invalid inputs!");
        }
        this.size = n;
        this.openNum = 0;
        // n^2 elements from the grid + virtualFront + virtualRear
        this.uf = new WeightedQuickUnionUF(n * n + 2);
        grid = new boolean[n * n + 2];
        for (int i = 0; i < n * n; i++) {
            grid[i] = false;
        }
    }

    /**
     * Convert index of the grid into that used in the UF
     *
     * @param row grid row index
     * @param col grid col index
     * @return index of internal parent array in the UF
     */
    private int index(int row, int col) {
        return (row - 1) * this.size + col;
    }

    /**
     * Validate that the index pair is appropriate
     *
     * @param row row index of a site
     * @param col col index of a site
     * @return if the indices are valid, return true, else false
     */
    private boolean validate(int row, int col) {
        if (row < 1 || row > this.size || col < 1 || col > this.size) {
            return false;
        }
        return true;
    }

    /**
     * Convert a site from closed to open, then add it to a connected component
     *
     * @param row row index of the site
     * @param col col index of the site
     */
    public void open(int row, int col) {
        if (!validate(row, col)) {
            throw new IllegalArgumentException("Invalid inputs!");
        }
        grid[index(row, col)] = true;
        this.openNum += 1;

        // if the site is in the first row, connect it to the virtualFront
        if (row == 1) {
            uf.union(index(row, col), 0);
        }
        // if the site is in the last row, connect it to the virtualRear
        if (row == this.size) {
            uf.union(index(row, col), this.size * this.size + 1);
        }
        // Union openSites in the top, bottom, left, right:
        // topSite
        if (validate(row - 1, col) && isOpen(row - 1, col)) {
            uf.union(index(row, col), index(row - 1, col));
        }
        // rightSite
        if (validate(row, col + 1) && isOpen(row, col + 1)) {
            uf.union(index(row, col), index(row, col + 1));
        }
        // bottomSite
        if (validate(row + 1, col) && isOpen(row + 1, col)) {
            uf.union(index(row, col), index(row + 1, col));
        }
        // leftSite
        if (validate(row, col - 1) && isOpen(row, col - 1)) {
            uf.union(index(row, col), index(row, col - 1));
        }
    }

    /**
     * Check whether a site is open
     *
     * @param row row index of the site
     * @param col col index of the site
     * @return if the site is open, return true, else false
     */
    public boolean isOpen(int row, int col) {
        if (!validate(row, col)) {
            throw new IllegalArgumentException("Invalid inputs!");
        }
        return grid[index(row, col)] == true;
    }

    /**
     * Check whether a site is full, i.e., whether a site is connected to an open site in the first
     * row
     *
     * @param row row index of the site
     * @param col col index of the site
     * @return if the site is full, return true, else false
     */
    public boolean isFull(int row, int col) {
        if (!validate(row, col)) {
            throw new IllegalArgumentException("Invalid inputs!");
        }
        // the connected method in WeightedQuickUnionUF is deprecated, but uf.connected(p, q) <=> uf.find(p) == uf.find(q)
        return uf.find(index(row, col)) == uf.find(0) && isOpen(row, col);
    }

    /**
     * Find the number of open sites in the grid
     *
     * @return the number of open sites
     */
    public int numberOfOpenSites() {
        return this.openNum;
    }

    /**
     * Check if the percolation system percolates, i.e., whether there is a path connecting an open
     * site in the first row to an open site in the last row
     *
     * @return if the system percolates, return true, else false
     */
    public boolean percolates() {
        return uf.find(0) == uf.find(this.size * this.size + 1);
    }

}
