/* *****************************************************************************
 *  Name: yanxuanshaozhu
 *  Date: 04/22/2021
 *  Description: coursera algorithm week 2 project
 **************************************************************************** */

import edu.princeton.cs.algs4.StdRandom;

import java.util.Iterator;
import java.util.NoSuchElementException;

public class RandomizedQueue<Item> implements Iterable<Item> {
    private Item[] items;
    private int index;

    // construct an empty randomized queue
    public RandomizedQueue() {
        items = (Item[]) new Object[2];
        index = 0;
    }

    // is the randomized queue empty?
    public boolean isEmpty() {
        return index == 0;
    }

    // return the number of items on the randomized queue
    public int size() {
        return index;
    }

    private void resize(int capacity) {
        Item[] copy = (Item[]) new Object[capacity];
        System.arraycopy(items, 0, copy, 0, index);
        items = copy;
    }

    // add the item
    public void enqueue(Item item) {
        if (item == null) {
            throw new IllegalArgumentException("Illegal argument!");
        }
        if (index == items.length) {
            resize(items.length * 2);
        }
        items[index] = item;
        index += 1;
    }

    // remove and return a random item
    public Item dequeue() {
        if (index == 0) {
            throw new NoSuchElementException("No such element!");
        }
        int randomNum = StdRandom.uniform(index);
        Item item = items[randomNum];
        items[randomNum] = items[index - 1];
        items[index - 1] = null;
        index -= 1;
        if (index == items.length / 4) {
            resize(items.length / 2);
        }
        return item;
    }

    // return a random item (but do not remove it)
    public Item sample() {
        if (index == 0) {
            throw new NoSuchElementException("No such element!");
        }
        int randomNum = StdRandom.uniform(index);
        return items[randomNum];
    }

    private class RandomizedQueueIterator implements Iterator<Item> {
        private int[] indices;
        private int remaining;

        public RandomizedQueueIterator() {
            indices = new int[index];
            for (int i = 0; i < index; i++) {
                indices[i] = i;
            }
            remaining = index;
        }

        public boolean hasNext() {
            if (remaining > 0) {
                return true;
            }
            return false;
        }

        public Item next() {
            if (!hasNext()) {
                throw new NoSuchElementException("No such element!");
            }
            int randNum = StdRandom.uniform(index);
            while (indices[randNum] < 0) {
                randNum = StdRandom.uniform(index);
            }
            Item item = items[randNum];
            indices[randNum] = remaining - 1;
            indices[remaining - 1] = -1;
            remaining -= 1;
            return item;
        }

        public void remove() {
            throw new UnsupportedOperationException("Unsupported operation!");
        }
    }

    // return an independent iterator over items in random order
    public Iterator<Item> iterator() {
        return new RandomizedQueueIterator();
    }

    // unit testing (required)
    public static void main(String[] args) {
        RandomizedQueue<String> randomizedQueue = new RandomizedQueue<String>();
        System.out.println(randomizedQueue.size());
        randomizedQueue.enqueue("item0");
        System.out.println(randomizedQueue.size());
        randomizedQueue.enqueue("item1");
        System.out.println(randomizedQueue.size());
        randomizedQueue.enqueue("item2");
        System.out.println(randomizedQueue.size());
        System.out.println("------------------------------");
        Iterator<String> iterator = randomizedQueue.iterator();
        while (iterator.hasNext()) {
            System.out.println(randomizedQueue.dequeue());
            System.out.println(randomizedQueue.size());
            if (randomizedQueue.size() == 0) {
                System.out.println("Empty que!");
                break;
            }
        }
    }

}