/* *****************************************************************************
 *  Name: yanxuanshaozhu
 *  Date: 04/22/2021
 *  Description: coursera algorithm week 2 project
 **************************************************************************** */

import java.util.Iterator;
import java.util.NoSuchElementException;

public class Deque<Item> implements Iterable<Item> {
    private Node first;
    private Node last;
    private int size;

    private class Node {
        Item item;
        Node prev;
        Node next;
    }

    /**
     * Construct an empty deque
     */
    public Deque() {
        first = null;
        last = null;
        size = 0;
    }

    /**
     * Check whether the deque is empty
     *
     * @return return true is the deque is empty, false otherwise
     */
    public boolean isEmpty() {
        return size == 0;
    }


    /**
     * Find the amount of nodes in the deque
     *
     * @return return the number of nodes in the deque
     */
    public int size() {
        return size;
    }

    /**
     * Add a node at the front of the deque
     *
     * @param item item of the node
     */
    public void addFirst(Item item) {
        if (item == null) {
            throw new IllegalArgumentException("Illegal argument!");
        }
        if (size == 0) {
            first = new Node();
            first.item = item;
            last = first;
        }
        else {
            Node node = new Node();
            node.item = item;
            node.next = first;
            first.prev = node;
            first = node;
        }
        size += 1;
    }

    /**
     * Add a node at the rear of the deque
     *
     * @param item item of the node
     */
    public void addLast(Item item) {
        if (item == null) {
            throw new IllegalArgumentException("Illegal argument!");
        }
        if (size == 0) {
            last = new Node();
            last.item = item;
            first = last;
        }
        else {
            Node node = new Node();
            node.item = item;
            node.prev = last;
            last.next = node;
            last = node;
        }
        size += 1;
    }

    /**
     * Remove a node at the front of the deque
     *
     * @return return the item of the removed node
     */
    public Item removeFirst() {
        if (size == 0) {
            throw new NoSuchElementException("No such element!");
        }
        else {
            Item item = first.item;
            first = first.next;
            if (first == null) {
                last = null;
            }
            else {
                first.prev = null;
            }
            size -= 1;
            return item;
        }
    }

    /**
     * Remove a node at the rear of the deque
     *
     * @return return the item of the removed node
     */
    public Item removeLast() {
        if (size == 0) {
            throw new NoSuchElementException("No such element!");
        }
        else {
            Item item = last.item;
            last = last.prev;
            if (last == null) {
                first = null;
            }
            else {
                last.next = null;
            }
            size -= 1;
            return item;
        }
    }


    private class DequeIterator implements Iterator<Item> {
        private Node current = first;

        public boolean hasNext() {
            return current != null;
        }

        public Item next() {
            if (!hasNext()) {
                throw new NoSuchElementException("No such element!");
            }
            else {
                Item item = current.item;
                current = current.next;
                return item;
            }
        }

        public void remove() {
            throw new UnsupportedOperationException("Unsupported operation!");
        }
    }

    /**
     * Construct an iterator for the deque
     *
     * @return return the generated iterator
     */
    public Iterator<Item> iterator() {
        return new DequeIterator();
    }

    public static void main(String[] args) {
        Deque<String> deque = new Deque<String>();
        System.out.println(deque);
        deque.addFirst("item0");
        deque.addLast("item1");
        System.out.println(deque.size());
        Iterator<String> iterator = deque.iterator();
        while (iterator.hasNext()) {
            System.out.println(iterator.next());
        }
        System.out.println(deque.removeFirst());
        System.out.println(deque.isEmpty());
        System.out.println(deque.size());
        System.out.println(deque.removeLast());
        System.out.println(deque.isEmpty());
        System.out.println(deque.size());
    }
}