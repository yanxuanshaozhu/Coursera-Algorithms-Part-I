/* *****************************************************************************
 *  Name: yanxuanshaozhu
 *  Date: 04/30/2021
 *  Description: coursera algorithm week 3 project
 **************************************************************************** */


import java.util.Arrays;
import java.util.LinkedList;

public class FastCollinearPoints {
    private int numberOfSegments = 0;
    private LinkedList<LineSegment> lineSegments = new LinkedList<LineSegment>();

    public FastCollinearPoints(Point[] points) {
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

        for (int i = 0; i < N; i++) {
            Point[] temp = new Point[N - i];
            for (int k = 0; k < N - i; k++) {
                temp[k] = points[k + i];
            }
            Arrays.sort(temp, points[i].slopeOrder());
            int j = 0;
            while (j < N - i - 2) {
                if (points[i].slopeTo(temp[j]) == points[i].slopeTo(temp[j + 1])
                        && points[i].slopeTo(temp[j]) == points[i].slopeTo(temp[j + 2])) {
                    LinkedList<Point> linkedsegments = new LinkedList<Point>();
                    linkedsegments.add(points[i]);
                    int k = j + 2;
                    int g = 3;
                    while (k < (points.length - i - 1) && points[i].slopeTo(temp[k]) == points[i]
                            .slopeTo(temp[k + 1])) {
                        k++;
                        g++;
                    }
                    while (k >= j) {
                        linkedsegments.add(temp[k]);
                        k--;
                    }
                    Point[] arraysegments = linkedsegments.toArray(new Point[g + 1]);
                    Arrays.sort(arraysegments);
                    LineSegment segment = new LineSegment(arraysegments[0], arraysegments[g]);
                    lineSegments.add(segment);
                    numberOfSegments++;
                    j = g + j;
                } else
                    j++;
            }
        }
    }

    // finds all line segments containing 4 or more points
    public int numberOfSegments() {
        return this.numberOfSegments;
    }

    // the number of line segments
    public LineSegment[] segments() {
        LineSegment[] segments = lineSegments.toArray(new LineSegment[this.numberOfSegments]);
        return segments;
    }

}