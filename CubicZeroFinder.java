import java.util.Scanner;

public class CubicZeroFinder {
	public static void main(String[] args) {
		Scanner sc = new Scanner(System.in);

		// Read in a, b, c, and d for our
		// equation ax^3 + bx^2 + cx + d = 0
		int a = Integer.parseInt(sc.nextLine());
		int b = Integer.parseInt(sc.nextLine());
		int c = Integer.parseInt(sc.nextLine());
		int d = Integer.parseInt(sc.nextLine());

		// Find the zeroes
		String[] zeroes = findZeroes(a, b, c, d);

		// Print the zeroes, one per line
		for (int i = 0; i < zeroes.length; i++) {
			System.out.println(zeroes[i]);
		}
	}

	/**
	 * Given a, b, c, and d for a cubic equation ax^3 + bx^2 + cx + d = 0, returns a
	 * sorted list of zeroes of this cubic equation.
	 */
	private static String[] findZeroes(int a, int b, int c, int d) {
		/*
		 * Given the coefficients of a cubic Returns the zeroes of the cubic accurate to
		 * 5 decimal places
		 */
		double doubleA = (double) a, doubleB = (double) b, doubleC = (double) c;
		double[] criticalPoints = findCriticalPoints(doubleA, doubleB, doubleC);
		String[] finalRoots = new String[3]; // solution
		String criticalPointType = "";
		double localMax = 0, localMin = 0;
		double rightCritPoint = Integer.MIN_VALUE, leftCritPoint = Integer.MAX_VALUE;
		int minCount = 0, maxCount = 0;

		for (double criticalPoint : criticalPoints) { // Parses each critical point
			// Finds type of critical point it is
			criticalPointType = findCriticalPointType(doubleA, doubleB, criticalPoint);

			// determines rightmost/leftmost critical points
			// essentially min/max finder
			if (criticalPoint > rightCritPoint) {
				rightCritPoint = criticalPoint;
			}
			if (criticalPoint < leftCritPoint) {
				leftCritPoint = criticalPoint;
			}

			if (criticalPointType.equals("sad")) {
				if (testIfRoot(a, b, c, d, criticalPoint) == 0) {
					// if saddle point on x axis, all roots are same
					String roots = String.format("%.5f", criticalPoint);
					return new String[] { roots, roots, roots };
				} else {
					// Terminates loop and goes to default argument:
					// looks for the 1 real root, then uses that root to find the other 2
					break;
				}
			} else if (criticalPointType.equals("max")) {
				maxCount++;
				localMax = criticalPoint;
				if (testIfRoot(a, b, c, d, criticalPoint) == 0) {
					// if the max is on the x axis, 2 are roots, tell program to find last root
					finalRoots[0] = String.format("%.5f", criticalPoint);
					finalRoots[1] = String.format("%.5f", criticalPoint);
					Object finalRoot = findBoundedRoot(rightCritPoint, leftCritPoint, a, b, c, d,
							Double.parseDouble(finalRoots[0]));
					finalRoots[2] = String.format("%.5f", finalRoot);
					finalRoots = sortFinalSolutions(finalRoots);
					return finalRoots;
				}

			} else if (criticalPointType.equals("min")) {
				minCount++;
				localMin = criticalPoint;
				if (testIfRoot(a, b, c, d, criticalPoint) == 0) {
					// if the min is on the x axis, 2 are roots, tell program to find last root
					finalRoots[0] = String.format("%.5f", criticalPoint);
					finalRoots[1] = String.format("%.5f", criticalPoint);
					Object finalRoot = findBoundedRoot(rightCritPoint, leftCritPoint, a, b, c, d,
							Double.parseDouble(finalRoots[0]));
					finalRoots[2] = String.format("%.5f", finalRoot);
					finalRoots = sortFinalSolutions(finalRoots);
					return finalRoots;
				}

			}
		}

		if (minCount > 0 && maxCount > 0) {
			// If there are 2 local extrema continue to find the roots
			if ((yValue(a, b, c, d, localMax) > 0 && yValue(a, b, c, d, localMin) > 0)) {
				// If both extrema are above the x axis
				// Find the 1 real root at the extrema,

				if (localMax > localMin) {
					finalRoots[0] = String.format("%.5f", binarySearchRoots(localMax, 100, a, b, c, d));
				} else {
					finalRoots[0] = String.format("%.5f", binarySearchRoots(-100, localMin, a, b, c, d));
				}

			} else if ((yValue(a, b, c, d, localMax) < 0 && yValue(a, b, c, d, localMin) < 0)) {
				// If both extrema are below the x axis
				// Find the 1 real extrema
				if (localMax < localMin) {
					finalRoots[0] = String.format("%.5f", binarySearchRoots(localMin, 100, a, b, c, d));
				} else {
					finalRoots[0] = String.format("%.5f", binarySearchRoots(-100, localMax, a, b, c, d));
				}

			} else {
				// Otherwise, do the finding 3 real roots protocol
				double max = rightCritPoint, min = leftCritPoint;
				if (yValue(a, b, c, d, leftCritPoint) > yValue(a, b, c, d, rightCritPoint)) {
					max = leftCritPoint;
					min = rightCritPoint;
				}
				finalRoots[0] = String.format("%.5f", binarySearchRoots(-100, leftCritPoint, a, b, c, d));
				finalRoots[1] = String.format("%.5f", binarySearchRoots(min, max, a, b, c, d));
				finalRoots[2] = String.format("%.5f", binarySearchRoots(rightCritPoint, 100, a, b, c, d));
				finalRoots = sortFinalSolutions(finalRoots);
				return finalRoots;
			}
		}
		// Default protocol:
		// Targeted at cubics with no extrema
		if (finalRoots[0] == null) {// if no other operation found the 1 real root, find it
			finalRoots[0] = String.format("%.5f", binarySearchRoots(-100, 100, a, b, c, d));
		}
		// Synthetic divide the cubic to a quadratic with the real root found
		double[] quadraticReduction = syntheticDivideValue(a, b, c, d, Double.parseDouble(finalRoots[0]));
		// Use quadratic formula to find imaginary roots, return them
		Object[][] quadraticRoots = quadraticFormula(quadraticReduction[0], quadraticReduction[1],
				quadraticReduction[2]);
		finalRoots[1] = String.format("%.5f", quadraticRoots[0][0]) + quadraticRoots[0][1];
		finalRoots[2] = String.format("%.5f", quadraticRoots[1][0]) + "+" + quadraticRoots[1][1];
		return finalRoots;

	}

	private static double testIfRoot(int a, int b, int c, int d, double root) {
		/*
		 * Returns true if the possible root is a root, false if not
		 */
		if (root < 1 && root > -1) {
			// Higher precision for likely tiny values
			return Double.parseDouble(String.format("%.8f", yValue(a, b, c, d, root)));
		} else {
			// Lower decimal precision for likely higher values
			return Double.parseDouble(String.format("%.5f", yValue(a, b, c, d, root)));
		}

	}

	private static double yValue(int a, int b, int c, int d, double x) {
		// returns f(x) where f is a cubic function with coefficients a, b, c, d
		return a * power(x, 3) + b * power(x, 2) + c * x + d;
	}

	private static String findCriticalPointType(double a, double b, double criticalPoint) {
		/*
		 * Given coefficients for 2nd derivative, returns min for local minimum max for
		 * local maximum, sad for saddlepoint
		 */
		double secondDerivativeOfX = a * 6 * criticalPoint + b * 2;

		if (secondDerivativeOfX > 0) {
			return "min";
		} else if (secondDerivativeOfX < 0) {
			return "max";
		} else {
			return "sad";
		}

	}

	private static double[] findCriticalPoints(double a, double b, double c) {
		/*
		 * Given coefficients of a quadratic (1st derivative), returns its roots
		 * (critical points)
		 */
		Object[][] rootsOfQuadratic = quadraticFormula(a * 3, b * 2, c);

		if ((rootsOfQuadratic[0][1].toString()).contains("i")) { // no real roots
			return new double[] {};
		} else if (rootsOfQuadratic[0][0] == rootsOfQuadratic[1][0]) {// multiplicity 2 root
			return new double[] { (double) rootsOfQuadratic[0][0] };
		} else { // unique roots
			return new double[] { (double) rootsOfQuadratic[0][0], (double) rootsOfQuadratic[1][0] };
		}

	}

	private static Object findBoundedRoot(double rightMost, double leftMost, int a, int b, int c, int d,
			double initRoot) {
		/*
		 * Given a rightmost and leftmost root and the initial root Checks if another
		 * unique root is between -100 and the leftmost root or the rightmost root and
		 * 100 Returns the the unique root
		 */

		Object possibleRoot = binarySearchRoots(rightMost, 100, a, b, c, d); // checks for root between rightmost and
																				// 100
		if (possibleRoot != null) { // if the root is there, ensure it isn't the multiplicity 2 root
			double doubleCopyPossibleRoot = (double) possibleRoot;
			if (Double.parseDouble(String.format("%.2f", doubleCopyPossibleRoot)) == initRoot) {
				// accounts for slight inaccuracy in binary search of discrete numbers
				// sometimes if the lower bound is the root,
				return binarySearchRoots(-100, leftMost, a, b, c, d);
			} else {
				// if this is a unique root, return it
				return possibleRoot;
			}
		}
		// otherwise return the root between -100 and the left most root
		return binarySearchRoots(-100, leftMost, a, b, c, d);
	}

	private static Object binarySearchRoots(double min, double max, int a, int b, int c, int d) {
		/*
		 * Binary searches for a root between the minimum and maximum value
		 */
		int numberIterations = 0;

		double midPoint = ((max - min) / 2.0 + min);
		while (numberIterations < 100) { // cuts off after 10^30 possible values accounted for
			numberIterations++;
			midPoint = ((max - min) / 2.0 + min);
			double approximateXValue = testIfRoot(a, b, c, d, midPoint);
			if (aboutEquals(approximateXValue, 0)) {
				return Double.parseDouble(String.format("%.5f", midPoint));
			} else if (approximateXValue > 0) {
				max = midPoint;
			} else if (approximateXValue < 0) {
				min = midPoint;
			}
		}
		return null;

	}

	private static boolean aboutEquals(double guess, double whatItShouldEqual) {
		/*
		 * Given a value and what it should equal, returns true if it about equals it
		 */
		final double epsilon = 0.00000001;
		if (guess == whatItShouldEqual) {
			return true;
		} else if (guess + epsilon == whatItShouldEqual) {
			return true;
		} else if (guess - epsilon == whatItShouldEqual) {
			return true;
		}
		return false;
	}

	private static double[] syntheticDivideValue(int a, int b, int c, int d, double divisor) {
		/*
		 * Given the coefficients of a cubic and a value to be divided out Returns a
		 * double[] of size three of a quadratic with value synthetically divided out a
		 * is double[0], b is double[1] c is double[2]
		 */
		double coeff1, coeff2, coeff3;

		coeff1 = a;
		coeff2 = b + divisor * coeff1;
		coeff3 = c + divisor * coeff2;
		// System.out.println(coeff1 + " " + coeff2 + " " + coeff3);
		return new double[] { round(coeff1), round(coeff2), round(coeff3) };
		// remainder must be 0 b/c value is a predetermined root of the cubic

	}

	private static double round(double a) {
		/*
		 * Given a double, returns a double with the accuracy of 5 (rounded) decimal
		 * places
		 */
		return Double.parseDouble(String.format("%.5f", a));
	}

	private static String[] sortFinalSolutions(String[] roots) {
		/*
		 * Returns the roots in sorted order
		 */
		double[] doubleRoots = new double[3];
		for (int i = 0; i < roots.length; i++) { // Converts each root to a sortable double
			try {
				doubleRoots[i] = Double.parseDouble(roots[i]);
			} catch (Exception e) {// should never happen, used for debugging
				System.out.println(i);
				continue;
			}
		}
		doubleRoots = ascendingInsertionSort(doubleRoots); // Sorts the double version of roots in ascending order
		for (int i = 0; i < roots.length; i++) { // Properly formats doubles as strings
			roots[i] = String.format("%.5f", doubleRoots[i]);
		}

		return roots;
	}

	private static double[] ascendingInsertionSort(double[] roots) {
		/*
		 * Given an array of values Sorts them in ascending order Returns the sorted
		 * array
		 */
		for (int i = 1; i < roots.length; i++) {
			for (int n = i; n > 0; n--) {
				if (roots[n] < roots[n - 1]) {
					double previousValue = roots[n - 1];
					roots[n - 1] = roots[n];
					roots[n] = previousValue;
				}
			}
		}
		return roots;
	}

	private static Object[][] quadraticFormula(double a, double b, double c) {
		/*
		 * Performs the quadratic formula on coefficients a, b, c Returns a
		 * dictionary-esk 2-d array with String.length Representing the number of roots,
		 * and String[].length representing A dictionary with the key being the root and
		 * value The type "r" if real, or "(value)i" where value is the imaginary
		 * component if imaginary
		 */
		int insideRadical = (int) (power(b, 2) - 4 * a * c);
		String[] radicalType = new String[2];
		double[] roots = new double[2];
		Object[][] answer = new Object[2][2];

		if (insideRadical < 0) { // checks for problematic imaginary numbers, marks if imaginary
			radicalType[0] = radicalType[1] = "i";
			insideRadical *= -1;
		} else {
			radicalType[0] = radicalType[1] = "r";
		}

		double radical = Double.parseDouble(findNthRoot(insideRadical, 2, 5)); // finds the sqrt of b^2 - 4ac
		if (radicalType[0].equals("r")) { // all is a real number
			roots[0] = ((-1 * b + radical) / (2 * a));
			roots[1] = ((-1 * b - radical) / (2 * a));
		} else { // part imaginary, part real
			roots[0] = roots[1] = -1 * b / (2 * a);
			radicalType[0] = String.format("%.5f", (-1 * radical / (2.0 * a))) + "i";
			radicalType[1] = String.format("%.5f", (radical / (2.0 * a))) + "i";
		}
		answer[0][0] = roots[0]; // could use a loop but this works
		answer[0][1] = radicalType[0];
		answer[1][0] = roots[1];
		answer[1][1] = radicalType[1];

		return answer;

	}

	private static String findNthRoot(int number, int n, int precision) {
		/*
		 * Given a int number, the number being rooted, int n, the power by which it is
		 * being rooted and int precision the number of decimals the root should be
		 * precise to Returns a string of the nth root of number to precision decimal
		 * places
		 */
		double upperBound = number;
		double lowerBound = 0.0;
		double midPoint;

		while (true) {
			midPoint = (upperBound - lowerBound) / 2 + lowerBound;
			if (Double.parseDouble(
					String.format("%." + Integer.toString(precision) + "f", power(midPoint, n))) == number) {
				return String.format("%." + Integer.toString(precision) + "f", midPoint);
			} else if (power(midPoint, n) > number) {
				upperBound = midPoint;
			} else {
				lowerBound = midPoint;
			}
		}
	}

	private static double power(double n, int rootNumber) {
		/*
		 * Given int n, the number to be raised to a power And int rootNumber, the power
		 * it is to be raised to Returns number to that power
		 */
		double result = n;
		for (int i = 1; i < rootNumber; i++) {
			result *= n;
		}
		return result;
	}
}