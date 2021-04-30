/* *****************************************************************************
 *  Name: yanxuanshaozhu
 *  Date: 04/27/2021
 *  Description:
 **************************************************************************** */

import java.util.Arrays;

public class BruteCollinearPoints {
    private LineSegment[] lineSegments;
    private int numberOfSegments;
    private Point[] points;


    private void resize(int capacity) {
        LineSegment[] copy = new LineSegment[capacity];
        System.arraycopy(this.lineSegments, 0, copy, 0, this.lineSegments.length);
        this.lineSegments = copy;
    }


    // finds all line segments containing 4 points
    public BruteCollinearPoints(Point[] points) {
        this.points = points;
        this.numberOfSegments = 0;
        this.lineSegments = new LineSegment[1];

        if (this.points == null) {
            throw new IllegalArgumentException(
                    "argument to BruteCollinearPoints constructor is null");
        }

        for (Point point : this.points) {
            if (point == null) {
                throw new IllegalArgumentException("input array contains null points");
            }
        }

        for (int i = 1; i < this.points.length; i++) {
            if (points[i - 1].compareTo(points[i]) == 0) {
                throw new IllegalArgumentException("input array contains identical points");
            }
        }

        Arrays.sort(this.points);

        for (int i = 0; i < this.points.length - 3; i++) {
            for (int j = i + 1; j < this.points.length - 2; j++) {
                for (int k = j + 1; k < this.points.length - 2; k++) {
                    for (int l = k + 1; l < this.points.length; l++) {
                        if (this.points[i].slopeTo(this.points[j]) == this.points[i]
                                .slopeTo(this.points[k])
                                && this.points[i].slopeTo(this.points[k]) == this.points[i]
                                .slopeTo(this.points[l])) {
                            LineSegment lineSegment = new LineSegment(points[i], points[l]);

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

    // the number of line segments
    public int numberOfSegments() {
        return this.numberOfSegments;
    }

    // the line segments
    public LineSegment[] segments() {
        return lineSegments;
    }


}
