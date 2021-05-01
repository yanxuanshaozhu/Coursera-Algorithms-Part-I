/* *****************************************************************************
 *  Name: yanxuanshaozhu
 *  Date: 04/27/2021
 *  Description: coursera algorithm week 3 project
 **************************************************************************** */

import java.util.Arrays;

public class BruteCollinearPoints {

    private final Point[] points;
    private LineSegment lineSegments[];
    private int numberOfSegments;


    public BruteCollinearPoints(Point[] points) {

        int N = points.length;

        if (points == null) {
            throw new IllegalArgumentException(
                    "argument to BruteCollinearPoints constructor is null");
        }

        for (int i = 0; i < points.length; i++) {
            if (points[i] == null) {
                throw new IllegalArgumentException("input array contains null points");
            }
        }

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                if (i != j && points[i].compareTo(points[j]) == 0) {
                    throw new IllegalArgumentException("input array contains identical points");
                }
            }
        }

        this.points = points.clone();
        this.lineSegments = new LineSegment[2];
        this.numberOfSegments = 0;


        Arrays.sort(this.points);

        for (int i = 0; i < this.points.length - 3; i++) {
            for (int j = i + 1; j < this.points.length - 2; j++) {
                for (int k = j + 1; k < this.points.length - 1; k++) {
                    for (int m = k + 1; m < this.points.length; m++) {
                        if (this.points[i].slopeTo(this.points[j]) == this.points[i]
                                .slopeTo(this.points[k]) &&
                                this.points[i].slopeTo(this.points[k]) == this.points[i]
                                        .slopeTo(this.points[m])) {


                            LineSegment lineSegment = new LineSegment(this.points[i],
                                                                      this.points[m]);
                            if (this.numberOfSegments == this.lineSegments.length) {
                                resize(this.lineSegments.length * 2);
                            }
                            this.lineSegments[this.numberOfSegments++] = lineSegment;

                        }
                    }
                }
            }
        }
    }

    private void resize(int capacity) {
        LineSegment[] temp = new LineSegment[capacity];
        System.arraycopy(this.lineSegments, 0, temp, 0, this.numberOfSegments);
        this.lineSegments = temp;
    }

    public int numberOfSegments() {
        return this.numberOfSegments;
    }


    public LineSegment[] segments() {
        return Arrays.copyOf(this.lineSegments, this.numberOfSegments);
    }


}