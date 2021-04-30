/* *****************************************************************************
 *  Name: yanxuanshaozhu
 *  Date: 04/30/2021
 *  Description:
 **************************************************************************** */

import java.util.Arrays;

public class FastCollinearPoints {
    private LineSegment[] lineSegments;
    private int numberOfSegments;
    private Point[] points;


    private void resize(int capacity) {
        LineSegment[] copy = new LineSegment[capacity];
        System.arraycopy(this.lineSegments, 0, copy, 0, this.lineSegments.length);
        this.lineSegments = copy;
    }


    // finds all line segments containing 4 or more points
    public FastCollinearPoints(Point[] points) {
        this.lineSegments = new LineSegment[1];
        this.numberOfSegments = 0;
        this.points = points;

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
        for (int i = 0; i < this.points.length; i++) {
            Arrays.sort(this.points, points[i].slopeOrder());
            int segmentLength = 1;



        }
    }

    // the number of line segments
    public int numberOfSegments() {
        return this.numberOfSegments;
    }


    // the line segments
    public LineSegment[] segments() {
        return this.lineSegments;
    }
}
