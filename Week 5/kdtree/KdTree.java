/* *****************************************************************************
 *  Name:  yanxuanshaozhu
 *  Date:  05/18/2021
 *  Description: coursera algorithm week 5 project
 **************************************************************************** */

import edu.princeton.cs.algs4.Point2D;
import edu.princeton.cs.algs4.RectHV;
import edu.princeton.cs.algs4.SET;
import edu.princeton.cs.algs4.StdDraw;

public class KdTree {
    private Node root;
    private int size;

    // construct an empty set of points
    public KdTree() {
        root = null;
        size = 0;
    }

    // is the set empty?
    public boolean isEmpty() {
        return size == 0;
    }

    // number of points in the set
    public int size() {
        return size;
    }

    private Node insert(Node node, Point2D point, boolean horizontal) {
        if (node == null) {
            size += 1;
            return new Node(point, horizontal);
        }
        if (node.point.equals(point)) {
            return node;
        }
        if (horizontal) {
            if (point.y() < node.point.y()) {
                node.left = insert(node.left, point, !horizontal);
            }
            else {
                node.right = insert(node.right, point, !horizontal);
            }
        }
        if (!horizontal) {
            if (point.x() < node.point.x()) {
                node.left = insert(node.left, point, horizontal);
            }
            else {
                node.right = insert(node.right, point, horizontal);
            }
        }
        return node;
    }

    // add the point to the set (if it is not already in the set)
    public void insert(Point2D p) {
        if (p == null) {
            throw new IllegalArgumentException("null argument");
        }
        root = insert(root, p, false);
    }

    // does the set contain point p?
    public boolean contains(Point2D p) {
        if (p == null) {
            throw new IllegalArgumentException("null argument!");
        }
        Node curr = root;
        while (curr != null) {
            if (curr.point.equals(p)) {
                return true;
            }
            if (curr.horizontal) {
                if (p.y() < curr.point.y()) {
                    curr = curr.left;
                }
                else {
                    curr = curr.right;
                }
            }
            if (!curr.horizontal) {
                if (p.x() < curr.point.x()) {
                    curr = curr.left;
                }
                else {
                    curr = curr.right;
                }
            }
        }
        return false;
    }

    private void draw(Node node, Node parent) {
        if (node == null) {
            return;
        }
        if (node.horizontal) {
            StdDraw.setPenColor(StdDraw.BLUE);
            if (parent == null) {
                node.point.drawTo(new Point2D(0, node.point.y()));
                node.point.drawTo(new Point2D(1, node.point.y()));
            }
            else if (node.point.x() < parent.point.x()) {
                node.point.drawTo(new Point2D(0, node.point.y()));
                node.point.drawTo(new Point2D(parent.point.x(), node.point.y()));
            }
            else {
                node.point.drawTo(new Point2D(1, node.point.y()));
                node.point.drawTo(new Point2D(parent.point.x(), node.point.y()));
            }
        }
        if (!node.horizontal) {
            StdDraw.setPenColor(StdDraw.BLUE);
            if (parent == null) {
                node.point.drawTo(new Point2D(node.point.x(), 0));
                node.point.drawTo(new Point2D(node.point.x(), 1));
            }
            else if (node.point.y() < parent.point.y()) {
                node.point.drawTo(new Point2D(node.point.x(), 0));
                node.point.drawTo(new Point2D(node.point.x(), parent.point.y()));
            }
            else {
                node.point.drawTo(new Point2D(node.point.x(), 1));
                node.point.drawTo(new Point2D(node.point.x(), parent.point.y()));
            }
        }
        StdDraw.setPenColor(StdDraw.BLACK);
        node.point.draw();
        draw(node.left, node);
        draw(node.right, node);
    }

    // draw all points to standard draw
    public void draw() {
        draw(root, null);
    }

    private void helper(RectHV rectHV, Node node, SET<Point2D> set) {
        if (node == null) {
            return;
        }
        if (rectHV.contains(node.point)) {
            set.add(node.point);
        }
        helper(rectHV, node.left, set);
        helper(rectHV, node.right, set);
    }

    // all points that are inside the rectangle (or on the boundary)
    public Iterable<Point2D> range(RectHV rect) {
        if (rect == null) {
            throw new IllegalArgumentException("null argument!");
        }
        SET<Point2D> points = new SET<>();
        helper(rect, root, points);
        return points;
    }

    private Point2D nearest(Point2D point, Point2D nearest, Node node) {
        if (node == null) {
            return nearest;
        }

        if (point.distanceTo(node.point) < point.distanceTo(nearest)) {
            nearest = node.point;
        }
        if (node.horizontal) {
            if (point.y() < node.point.y()) {
                nearest = nearest(point, nearest, node.left);
                nearest = nearest(point, nearest, node.right);
            }
            else {
                nearest = nearest(point, nearest, node.right);
                nearest = nearest(point, nearest, node.left);
            }
        }
        if (!node.horizontal) {
            if (point.x() < node.point.x()) {
                nearest = nearest(point, nearest, node.left);
                nearest = nearest(point, nearest, node.right);
            }
            else {
                nearest = nearest(point, nearest, node.right);
                nearest = nearest(point, nearest, node.left);
            }
        }
        return nearest;
    }

    // a nearest neighbor in the set to point p; null if the set is empty
    public Point2D nearest(Point2D p) {
        return nearest(p, root.point, root);
    }

    private class Node {
        private final Point2D point;
        private Node left;
        private Node right;
        private final boolean horizontal;

        public Node(Point2D point, boolean horizontal) {
            this.point = point;
            this.horizontal = horizontal;
            this.left = null;
            this.right = null;
        }

    }
}
