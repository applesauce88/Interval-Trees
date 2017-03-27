package structures;

import java.util.ArrayList;

/**
 * Encapsulates an interval tree.
 * 
 * @author runb-cs112
 */
public class IntervalTree {

	/**
	 * The root of the interval tree
	 */
	IntervalTreeNode root;

	/**
	 * Constructs entire interval tree from set of input intervals. Constructing
	 * the tree means building the interval tree structure and mapping the
	 * intervals to the nodes.
	 * 
	 * @param intervals
	 *            Array list of intervals for which the tree is constructed
	 */
	public IntervalTree(ArrayList<Interval> intervals) {

		// make a copy of intervals to use for right sorting
		ArrayList<Interval> intervalsRight = new ArrayList<Interval>(intervals.size());
		for (Interval iv : intervals) {
			intervalsRight.add(iv);
		}

		// rename input intervals for left sorting
		ArrayList<Interval> intervalsLeft = intervals;

		// sort intervals on left and right end points
		sortIntervals(intervalsLeft, 'l');
		sortIntervals(intervalsRight, 'r');

		// get sorted list of end points without duplicates
		ArrayList<Integer> sortedEndPoints = getSortedEndPoints(intervalsLeft, intervalsRight);
		System.out.println("sortedEndPoints.size" + sortedEndPoints.size());
		// build the tree nodes
		root = buildTreeNodes(sortedEndPoints);

		// map intervals to the tree nodes
		mapIntervalsToTree(intervalsLeft, intervalsRight);
	}

	/**
	 * Returns the root of this interval tree.
	 * 
	 * @return Root of interval tree.
	 */
	public IntervalTreeNode getRoot() {
		return root;
	}

	/**
	 * Sorts a set of intervals in place, according to left or right endpoints.
	 * At the end of the method, the parameter array list is a sorted list.
	 * 
	 * @param intervals
	 *            Array list of intervals to be sorted.
	 * @param lr
	 *            If 'l', then sort is on left endpoints; if 'r', sort is on
	 *            right endpoints
	 */
	public static void sortIntervals(ArrayList<Interval> intervals, char lr) { // COMPLETE
																				// THIS
																				// METHOD
		ArrayList<Interval> temp = new ArrayList<Interval>(intervals.size());

		if (lr == 'l') {
			while (intervals.size() > 0) {
				int t = intervals.get(0).leftEndPoint;
				int r = intervals.get(0).rightEndPoint;
				String desc = intervals.get(0).desc;
				int index = 0;
				for (int i = 1; i < intervals.size(); i++) {
					if (intervals.get(i).leftEndPoint < t) {
						t = intervals.get(i).leftEndPoint;
						r = intervals.get(i).rightEndPoint;
						desc = intervals.get(i).desc;
						index = i;
					}
				}
				temp.add(new Interval(t, r, desc));
				intervals.remove(index);
			}

			for (int i = 0; i < temp.size(); i++) {
				System.out.println("-" + temp.get(i).toString());
				intervals.add(temp.get(i));
			}
		} else if (lr == 'r') {
			while (intervals.size() > 0) {
				int t = intervals.get(0).rightEndPoint;
				int r = intervals.get(0).leftEndPoint;
				String desc = intervals.get(0).desc;
				int index = 0;
				for (int i = 1; i < intervals.size(); i++) {
					if (intervals.get(i).rightEndPoint < t) {
						t = intervals.get(i).rightEndPoint;
						r = intervals.get(i).leftEndPoint;
						desc = intervals.get(i).desc;
						index = i;
					}
				}
				temp.add(new Interval(r, t, desc));
				intervals.remove(index);
			}

			for (int i = 0; i < temp.size(); i++) {
				System.out.println("-" + temp.get(i).toString());
				intervals.add(temp.get(i));
			}

		}

	}

	/**
	 * Given a set of intervals (left sorted and right sorted), extracts the
	 * left and right end points, and returns a sorted list of the combined end
	 * points without duplicates.
	 * 
	 * @param leftSortedIntervals
	 *            Array list of intervals sorted according to left endpoints
	 * @param rightSortedIntervals
	 *            Array list of intervals sorted according to right endpoints
	 * @return Sorted array list of all endpoints without duplicates
	 */
	public static ArrayList<Integer> getSortedEndPoints(ArrayList<Interval> leftSortedIntervals,
			ArrayList<Interval> rightSortedIntervals) { // COMPLETE THIS METHOD
		// THE FOLLOWING LINE HAS BEEN ADDED TO MAKE THE PROGRAM COMPILE
		ArrayList<Integer> temp = new ArrayList<Integer>();
		int i = 0, j = 0;
		System.out.println("size left, size right" + leftSortedIntervals.size() + ", " + rightSortedIntervals.size());
		while (i < leftSortedIntervals.size() && j < rightSortedIntervals.size()) {
			if (leftSortedIntervals.get(i).leftEndPoint < rightSortedIntervals.get(j).rightEndPoint) {
				temp.add(leftSortedIntervals.get(i).leftEndPoint);
				i++;
				// System.out.println("A");
			} else if (leftSortedIntervals.get(i).leftEndPoint > rightSortedIntervals.get(j).rightEndPoint) {
				temp.add(rightSortedIntervals.get(j).rightEndPoint);
				j++;
				// System.out.println("B");
			} else if (leftSortedIntervals.get(i).leftEndPoint == rightSortedIntervals.get(j).rightEndPoint) {
				temp.add(rightSortedIntervals.get(j).rightEndPoint);
				i++;
				j++;
				// System.out.println("C");
			}

		}
		if (i == leftSortedIntervals.size()) {
			while (j < rightSortedIntervals.size()) {
				temp.add(rightSortedIntervals.get(j).rightEndPoint);
				j++;
			}
		} else if (j == rightSortedIntervals.size()) {
			while (i < leftSortedIntervals.size()) {
				temp.add(leftSortedIntervals.get(i).leftEndPoint);
				i++;
			}
		}
		for (int k = 1; k < temp.size(); k++) {
			// System.out.println("compare: "+temp.get(k)+" :
			// "+(temp.get(k-1)));
			if (temp.get(k).equals(temp.get(k - 1))) {
				temp.remove(k);
				k--;
				// System.out.println(temp.size());
			}
		}
		System.out.print("temp result: ");
		for (int k = 0; k < temp.size(); k++) {
			System.out.print(temp.get(k) + ", ");
		}
		System.out.println();
		return temp;
	}

	/**
	 * Builds the interval tree structure given a sorted array list of end
	 * points without duplicates.
	 * 
	 * @param endPoints
	 *            Sorted array list of end points
	 * @return Root of the tree structure
	 */
	public static IntervalTreeNode buildTreeNodes(ArrayList<Integer> endPoints) { // takes
																					// in
																					// sortedEndPoints
		// COMPLETE THIS METHOD
		// THE FOLLOWING LINE HAS BEEN ADDED TO MAKE THE PROGRAM COMPILE
		return buildTNodes(endPoints);
	}

	// use next method for recursion

	private static IntervalTreeNode buildTNodes(ArrayList<Integer> endPoints) {
		Queue<IntervalTreeNode> q = new Queue<IntervalTreeNode>();
		IntervalTreeNode T;
		System.out.println("build the first row");
		for (int i = 0; i < endPoints.size(); i++) {
			T = new IntervalTreeNode(endPoints.get(i).floatValue(), endPoints.get(i).floatValue(),
					endPoints.get(i).floatValue());
			System.out.print(endPoints.get(i).floatValue() + ", ");
			T.leftIntervals = new ArrayList<Interval>();
			T.rightIntervals = new ArrayList<Interval>();
			q.enqueue(T);
		}

		System.out.println();

		int temps = q.size;
		while (temps != 1) {
			IntervalTreeNode T1;
			IntervalTreeNode T2;
			float v1, v2;
			System.out.println("start next layer with this many things: " + q.size);
			System.out.println("temps: " + temps);
			while (temps > 1) {
				T1 = q.dequeue();
				T2 = q.dequeue();
				v1 = T1.maxSplitValue;
				System.out.print(" max: " + v1);
				v2 = T2.minSplitValue;
				System.out.print(" min: " + v2 + ", ");
				float x = (v1 + v2) / 2;
				IntervalTreeNode N = new IntervalTreeNode(x, T1.minSplitValue, T2.maxSplitValue);
				N.leftIntervals = new ArrayList<Interval>();
				N.rightIntervals = new ArrayList<Interval>();
				N.leftChild = T1;
				N.rightChild = T2;
				System.out.print(N.splitValue + ": ");
				q.enqueue(N);
				// N as root, T1 as left child of N, and T2 as right child of N
				temps = temps - 2;
			}
			if (temps == 1) {
				IntervalTreeNode a = q.dequeue();
				System.out.print(a.splitValue + ", " + a.maxSplitValue + ", " + a.minSplitValue + ": ");
				q.enqueue(a);
			}
			System.out.println();
			temps = q.size;
		}
		System.out.println("finished");
		T = q.dequeue();
		return T;

	}

	/**
	 * Maps a set of intervals to the nodes of this interval tree.
	 * 
	 * @param leftSortedIntervals
	 *            Array list of intervals sorted according to left endpoints
	 * @param rightSortedIntervals
	 *            Array list of intervals sorted according to right endpoints
	 */
	public void mapIntervalsToTree(ArrayList<Interval> leftSortedIntervals, ArrayList<Interval> rightSortedIntervals) {
		// COMPLETE THIS METHOD
		IntervalTreeNode T = root;

		for (int i = 0; i < leftSortedIntervals.size(); i++) {
			float a = leftSortedIntervals.get(i).leftEndPoint;
			float b = leftSortedIntervals.get(i).rightEndPoint;
			System.out.println(i + " " + leftSortedIntervals.get(i));
			while (T.splitValue < a || T.splitValue > b) {
				System.out.println("A " + T.splitValue);
				if (T.splitValue > a) {
					T = T.leftChild;
					System.out.println("B " + T.splitValue);
				} else if (T.splitValue < b) {
					T = T.rightChild;
					System.out.println("C " + T.splitValue);
				}
			}
			System.out.println("GONE THROUGH with i=:" + i);
			System.out.println(i + " " + leftSortedIntervals.get(i));
			Interval e = leftSortedIntervals.get(i);
			T.leftIntervals.add(e);
			T = root;
		}
		for (int i = 0; i < rightSortedIntervals.size(); i++) {
			float a = rightSortedIntervals.get(i).leftEndPoint;
			float b = rightSortedIntervals.get(i).rightEndPoint;
			System.out.println(i + " " + rightSortedIntervals.get(i));
			while (T.splitValue < a || T.splitValue > b) {
				System.out.println("A " + T.splitValue);
				if (T.splitValue > a) {
					T = T.leftChild;
					System.out.println("B " + T.splitValue);
				} else if (T.splitValue < b) {
					T = T.rightChild;
					System.out.println("C " + T.splitValue);
				}
			}
			System.out.println("GONE THROUGH with i=:" + i);
			System.out.println(i + " " + rightSortedIntervals.get(i));
			Interval e = rightSortedIntervals.get(i);
			T.rightIntervals.add(e);
			T = root;
		}
		// print out everything
		//printTree(root);
	}

	private void printTree(IntervalTreeNode q) {
		System.out.print(q.leftIntervals+ "            ");
		System.out.println(q.rightIntervals);
		if (q.leftChild == null) {
			return;
		}

		printTree(q.leftChild);
		printTree(q.rightChild);
	}

	/**
	 * Gets all intervals in this interval tree that intersect with a given
	 * interval.
	 * 
	 * @param q
	 *            The query interval for which intersections are to be found
	 * @return Array list of all intersecting intervals; size is 0 if there are
	 *         no intersections
	 */
	public ArrayList<Interval> findIntersectingIntervals(Interval q) {
		IntervalTreeNode R = root;
		ArrayList<Interval> resultList = new ArrayList<Interval>();
		findII(R, q, resultList);
		return resultList;
		// COMPLETE THIS METHOD
	}

	private void findII(IntervalTreeNode R, Interval q, ArrayList<Interval> resultList) {

		ArrayList<Interval> Llist;
		ArrayList<Interval> Rlist;

		if ((R.leftChild == null || R.rightChild == null)) {
			return;
		}

		Llist = R.leftIntervals;
		Rlist = R.rightIntervals;
		if (R.splitValue > q.leftEndPoint && R.splitValue < q.rightEndPoint) {
			System.out.println("A: " + R.splitValue + ": " + q.leftEndPoint + ", " + q.rightEndPoint);
			add(resultList, Llist);
			System.out.println("add: all here");
			findII(R.leftChild, q, resultList);
			findII(R.rightChild, q, resultList);
		} else if (R.splitValue < q.leftEndPoint) {
			System.out.println("B: " + R.splitValue + ": " + q.leftEndPoint + ", " + q.rightEndPoint);
			int i = Rlist.size() - 1;
			while (i >= 0 && intersect(Rlist.get(i), q)) {
				resultList.add(Rlist.get(i));
				System.out.println("add: " + Rlist.get(i));
				i--;
			}
			findII(R.rightChild, q, resultList);
		} else if (R.splitValue > q.rightEndPoint) {
			System.out.println("C: " + R.splitValue + ": " + q.leftEndPoint + ", " + q.rightEndPoint + "i<"+Llist.size());
			int i = 0;
			while (i < Llist.size() && intersect(Llist.get(i), q)) {
				resultList.add(Llist.get(i));
				System.out.println("add: " + Llist.get(i));
				i++;
			}
			findII(R.leftChild, q, resultList);
		}

		return;
	}

	private void add(ArrayList<Interval> resultList, ArrayList<Interval> otherList) {
		for (int i = 0; i < otherList.size(); i++) {
			resultList.add(otherList.get(i));
		}
	}

	private boolean intersect(Interval a, Interval b) {
		float al = a.leftEndPoint;
		float ar = a.rightEndPoint;
		float bl = b.leftEndPoint;
		float br = b.rightEndPoint;
		if (al <= bl && ar >= br) { //encapsulation
			return true;
		} else if (al >= bl && ar <= br) { //encapsulation
			return true;
		} else if (al <= bl && bl <= ar && ar <= br) {
			return true;
		} else if (bl <= al && al <= br && br <= ar) {
			return true;
		} else {
			return false;
		}
	}

}
