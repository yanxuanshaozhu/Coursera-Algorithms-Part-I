/* *****************************************************************************
 *  Name:              yanxuanshaozhu
 *  Coursera User ID:  yanxuanshaozhu
 *  Last modified:     April 16 2021
 **************************************************************************** */

import edu.princeton.cs.algs4.StdIn;
import edu.princeton.cs.algs4.StdOut;
import edu.princeton.cs.algs4.StdRandom;

public class RandomWord {
    public static void main(String[] args) {
        double prob = 0;
        String champion = null;
        int count = 0;
        while (!StdIn.isEmpty()) {
            String word = StdIn.readString();
            count += 1;
            prob = 1.0 / count;
            if (StdRandom.bernoulli(prob)) {
                champion = word;
            }
        }
        StdOut.println(champion);
    }
}
