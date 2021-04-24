/* *****************************************************************************
 *  Name:              yanxuanshaozhu
 *  Coursera User ID:  yanxuanshaozhu
 *  Last modified:     04/18/2021
 **************************************************************************** */

import edu.princeton.cs.algs4.StdRandom;
import edu.princeton.cs.algs4.StdStats;

public class PercolationStats {
    private final double[] results;
    private final int trials;

    /**
     * Calculate the p threshold for the percolation system with size n
     *
     * @param n size of the percolation system
     * @return p threshold for the percolation system
     */
    private double simulation(int n) {
        int threshold = 0;
        Percolation percolation = new Percolation(n);
        while (!percolation.percolates()) {
            // generate a random number in [1, n]
            // int row = 1 + (int) (Math.random() * n);
            // int col = 1 + (int) (Math.random() * n);
            int row = StdRandom.uniform(1, n + 1);
            int col = StdRandom.uniform(1, n + 1);
            if (!percolation.isOpen(row, col)) {
                percolation.open(row, col);
                threshold += 1;
            }
        }
        return 1.0 * threshold / (n * n);
    }
    /**
     * Constructor for the simulation process
     *
     * @param n      size for the percolation system in each simulation
     * @param trials amount of simulation conducted
     */
    public PercolationStats(int n, int trials) {
        if (n <= 0 || trials <= 0) {
            throw new IllegalArgumentException("invalid inputsss");
        }
        results = new double[trials];
        this.trials = trials;
        for (int i = 0; i < trials; i++) {
            results[i] = simulation(n);
        }

    }
    // sample mean of percolation threshold
    public double mean() {
        return StdStats.mean(results);
    }

    // sample standard deviation of percolation threshold
    public double stddev() {
        return StdStats.stddev(results);
    }

    // low endpoint of 95% confidence interval
    public double confidenceLo() {
        return StdStats.mean(results) - 1.96 * stddev() / Math.sqrt(this.trials);
    }

    // high endpoint of 95% confidence interval
    public double confidenceHi() {
        return StdStats.mean(results) + 1.96 * stddev() / Math.sqrt(this.trials);
    }

    // test client (see below)
    public static void main(String[] args) {
        int n = Integer.parseInt(args[0]);
        int trials = Integer.parseInt(args[1]);
        PercolationStats percolationStats = new PercolationStats(n, trials);
        System.out.println("mean                    = " + percolationStats.mean());
        System.out.println("stddev                  = " + percolationStats.stddev());
        System.out.println("95% confidence interval = " + "[" + percolationStats.confidenceLo() + ", " + percolationStats.confidenceHi() + "]");
    }
}
