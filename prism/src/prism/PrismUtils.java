//==============================================================================
//	
//	Copyright (c) 2002-
//	Authors:
//	* Dave Parker <david.parker@comlab.ox.ac.uk> (University of Oxford, formerly University of Birmingham)
//	
//------------------------------------------------------------------------------
//	
//	This file is part of PRISM.
//	
//	PRISM is free software; you can redistribute it and/or modify
//	it under the terms of the GNU General Public License as published by
//	the Free Software Foundation; either version 2 of the License, or
//	(at your option) any later version.
//	
//	PRISM is distributed in the hope that it will be useful,
//	but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//	GNU General Public License for more details.
//	
//	You should have received a copy of the GNU General Public License
//	along with PRISM; if not, write to the Free Software Foundation,
//	Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//	
//==============================================================================

package prism;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.PrimitiveIterator;
import java.util.Map.Entry;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Various general-purpose utility methods in Java
 */
public class PrismUtils
{
	// Threshold for comparison of doubles
	public static double epsilonDouble = 1e-12;

	/**
	 * Compute logarithm of x to base b.
	 */
	public static double log(double x, double b)
	{
		// If base is <=0 or ==1 (or +Inf/NaN), then result is NaN
		if (b <= 0 || b == 1 || (Double.isInfinite(b)) || Double.isNaN(b))
			return Double.NaN;

		// Otherwise, log_b (x) is log(x) / log(b)
		return Math.log(x) / Math.log(b);
	}

	/**
	 * Compute logarithm of x to base 2.
	 */
	public static double log2(double x)
	{
		return Math.log(x) / Math.log(2);
	}

	/**
	 * See if two doubles are within epsilon of each other (absolute error).
	 */
	public static boolean doublesAreCloseAbs(double d1, double d2, double epsilon)
	{
		// Deal with infinite cases
		if (Double.isInfinite(d1)) {
			return Double.isInfinite(d2) && (d1 > 0) == (d2 > 0);
		} else if (Double.isInfinite(d2)) {
			return false;
		}
		// Compute/check error
		return (Math.abs(d1 - d2) < epsilon);
	}

	/**
	 * See if two doubles are within epsilon of each other (relative error).
	 */
	public static boolean doublesAreCloseRel(double d1, double d2, double epsilon)
	{
		// Deal with infinite cases
		if (Double.isInfinite(d1)) {
			return Double.isInfinite(d2) && (d1 > 0) == (d2 > 0);
		} else if (Double.isInfinite(d2)) {
			return false;
		}
		// Compute/check error
		d1 = Math.abs(d1);
		d2 = Math.abs(d2);
		// For two (near) zero values, return true, for just one, return false
		if (d1 < epsilonDouble)
			return (d2 < epsilonDouble);
		return (Math.abs(d1 - d2) / d1 < epsilon);
	}

	/**
	 * See if two doubles are within epsilon of each other (relative or absolute error).
	 * @param abs Absolute if true, relative if false
	 */
	public static boolean doublesAreClose(double d1, double d2, double epsilon, boolean abs)
	{
		if (abs) {
			return doublesAreCloseAbs(d1, d2, epsilon);
		} else {
			return doublesAreCloseRel(d1, d2, epsilon);
		}
	}

	/**
	 * See if two arrays of doubles are all within epsilon of each other (relative or absolute error).
	 */
	public static boolean doublesAreClose(double d1[], double d2[], double epsilon, boolean abs)
	{
		int i, n;
		n = Math.min(d1.length, d2.length);
		if (abs) {
			for (i = 0; i < n; i++) {
				if (!PrismUtils.doublesAreCloseAbs(d1[i], d2[i], epsilon))
					return false;
			}
		} else {
			for (i = 0; i < n; i++) {
				if (!PrismUtils.doublesAreCloseRel(d1[i], d2[i], epsilon))
					return false;
			}
		}
		return true;
	}
	
	/**
	 * See if two arrays of doubles are all within epsilon of each other (relative or absolute error).
	 */
	public static boolean hashMapsAreClose(HashMap map1, HashMap map2, double epsilon, boolean abs)
	{
		if (map1.size()!=map2.size()){
			return false;
		}
		Set<Entry> entries= map1.entrySet();
		for (Entry entry : entries)
		{
			double d1=(double) entry.getValue();
			double d2;
			if(map2.get(entry.getKey())!=null){
				d2=(double) map2.get(entry.getKey());
				if (abs) {				
					if (!PrismUtils.doublesAreCloseAbs(d1, d2, epsilon))
						return false;				
				} else {				
					if (!PrismUtils.doublesAreCloseRel(d1, d2, epsilon))
						return false;				
				}
			}else {
				return false;
			}
		}		
		return true;
	}


	/**
	 * See if, for all the entries given by the {@code indizes}
	 * iterator, two arrays of doubles are all within epsilon of each other (relative or absolute error).
	 * <br>
	 * Considers Inf == Inf and -Inf == -Inf.
	 */
	public static boolean doublesAreClose(double d1[], double d2[], PrimitiveIterator.OfInt indizes, double epsilon, boolean abs)
	{
		if (abs) {
			while (indizes.hasNext()) {
				int i = indizes.nextInt();
				if (!PrismUtils.doublesAreCloseAbs(d1[i], d2[i], epsilon))
					return false;
			}
		} else {
			while (indizes.hasNext()) {
				int i = indizes.nextInt();
				if (!PrismUtils.doublesAreCloseRel(d1[i], d2[i], epsilon))
					return false;
			}
		}
		return true;
	}

	/**
	 * Measure supremum norm, either absolute or relative,
	 * return the maximum difference.
	 */
	public static double measureSupNorm(double[] d1, double[] d2, boolean abs)
	{
		int n = d1.length;
		assert( n == d2.length);

		double value = 0;
		if (abs) {
			for (int i=0; i < n; i++) {
				double diff = measureSupNormAbs(d1[i], d2[i]);
				if (diff > value)
					value = diff;
			}
		} else {
			for (int i=0; i < n; i++) {
				double diff = measureSupNormRel(d1[i], d2[i]);
				if (diff > value)
					value = diff;
			}
		}
		return value;
	}

	/**
	 * Measure supremum norm for two values, absolute.
	 */
	public static double measureSupNormAbs(double d1, double d2)
	{
		if (Double.isInfinite(d1) && d1==d2)
			return 0;
		return Math.abs(d1 - d2);
	}

	/**
	 * Measure supremum norm for two values, relative,
	 * with the first value used as the divisor.
	 */
	public static double measureSupNormRel(double d1, double d2)
	{
		if (d1==d2)
			return 0;
		return (Math.abs(d1 - d2) / d1);
	}

	/**
	 * Measure supremum norm, for all the entries given by the {@code indizes}
	 * iterator, for an interval iteration.
	 */
	public static double measureSupNormInterval(double[] lower, double[] upper, boolean abs, PrimitiveIterator.OfInt indizes)
	{
		int n = lower.length;
		assert( n== upper.length);

		double value = 0;
		while (indizes.hasNext()) {
			int i = indizes.nextInt();
			double diff = measureSupNormInterval(lower[i], upper[i], abs);
			if (diff > value)
				value = diff;
		}
		return value;
	}

	/**
	 * Measure supremum norm, either absolute or relative, for an interval iteration,
	 * return the maximum difference.
	 */
	public static double measureSupNormInterval(double[] lower, double[] upper, boolean abs)
	{
		int n = lower.length;
		assert( n== upper.length);

		double value = 0;
		for (int i=0; i < n; i++) {
			double diff = measureSupNormInterval(lower[i], upper[i], abs);
			if (diff > value)
				value = diff;
		}
		return value;
	}

	/**
	 * Measure supremum norm for two values, for interval iteration.
	 */
	public static double measureSupNormInterval(double lower, double upper, boolean abs)
	{
		// Deal with infinite cases
		if (Double.isInfinite(lower)) {
			if (Double.isInfinite(upper) && (lower > 0) == (upper > 0)) {
				return 0;
			} else {
				return Double.POSITIVE_INFINITY;
			}
		} else if (Double.isInfinite(upper)) {
			return Double.POSITIVE_INFINITY;
		}

		if (lower == upper)
			return 0.0;

		// Compute/check error
		lower = Math.abs(lower);
		upper = Math.abs(upper);
		double result = upper - lower;
		result = Math.abs(result);
		if (!abs) {
			result = result / lower;
		}
		return result;
	}

	/**
	 * See if two doubles are (nearly) equal.
	 */
	public static boolean doublesAreEqual(double d1, double d2)
	{
		return doublesAreCloseAbs(d1, d2, epsilonDouble);
	}

	/**
	 * Return the maximum finite value in a double array, looking at
	 * those entries with indices given bit the integer iterator.
	 */
	public static double findMaxFinite(double[] soln, PrimitiveIterator.OfInt indices)
	{
		double max_v = Double.NEGATIVE_INFINITY;
		while (indices.hasNext()) {
			int i = indices.nextInt();

			double v = soln[i];
			if (v < Double.POSITIVE_INFINITY) {
				max_v = Double.max(v, max_v);
			}
		}
		return max_v;
	}

	/**
	 * Ensure monotonicity from below for interval iteration solution vectors.
	 * Compares the old and new values and overwrites the new value with the old
	 * value if that is larger.
	 * @param old_values old solution vector
	 * @param new_values new solution vector
	 */
	public static void ensureMonotonicityFromBelow(double[] old_values, double[] new_values)
	{
		assert(old_values.length == new_values.length);

		for (int i = 0, n = old_values.length; i < n; i++) {
			double old_value = old_values[i];
			double new_value = new_values[i];
			// from below: do max
			if (old_value > new_value) {
				new_values[i] = old_value;
			}
		}
	}

	/**
	 * Ensure monotonicity from above for interval iteration solution vectors.
	 * Compares the old and new values and overwrites the new value with the old
	 * value if that is smaller.
	 * @param old_values old solution vector
	 * @param new_values new solution vector
	 */
	public static void ensureMonotonicityFromAbove(double[] old_values, double[] new_values)
	{
		assert(old_values.length == new_values.length);

		for (int i = 0, n = old_values.length; i < n; i++) {
			double old_value = old_values[i];
			double new_value = new_values[i];
			// from above: do min
			if (old_value < new_value) {
				new_values[i] = old_value;
			}
		}
	}

	/**
	 * Check for monotonicity: If the new_values are not element-wise less-than-equal the older values
	 * (for from_above == true), then throws an exception. If from_above == false, the logic is reversed,
	 * i.e., if the new_value is not greater-than-equal, an exception is thrown.
	 * @param old_values the old values
	 * @param new_values the new values
	 * @param from_above the direction
	 */
	public static void checkMonotonicity(double[] old_values, double[] new_values, boolean from_above) throws PrismException
	{
		assert(old_values.length == new_values.length);

		for (int i = 0, n = old_values.length; i < n; i++) {
			double old_value = old_values[i];
			double new_value = new_values[i];
			if (from_above && old_value < new_value) {
				throw new PrismException("Monotonicity violated (from above): old value " + old_value + " < new value " + new_value);
			}
			if (!from_above && old_value > new_value) {
				throw new PrismException("Monotonicity violated (from below): old value " + old_value + " > new value " + new_value);
			}
		}
	}

	/**
	 * Select midpoint from two interval iteration solution vectors.
	 * Stores the result in soln_below.
	 * @param soln_below solution vector from below
	 * @param soln_above solution vector from above
	 */
	public static void selectMidpoint(double[] soln_below, double[] soln_above)
	{
		assert(soln_below.length == soln_above.length);

		for (int i = 0, n = soln_below.length; i < n; i++) {
			double below = soln_below[i];
			double above = soln_above[i];

			if (below != above) {
				// use below + ( above - below ) / 2 instead of (below+above)/2 for better numerical
				// stability
				double d = below + (above - below)/2.0;
				if (d >= below && d <= above) {
					// only store result if between below and above
					// guard against rounding problems,
					// fallback is to simply return below as is
					soln_below[i] = d;
				}
			}
		}
	}

	/**
	 * Format a large integer, represented by a double, as a string. 
	 */
	public static String bigIntToString(double d)
	{
		if (d <= Long.MAX_VALUE) {
			return "" + Math.round(d);
		} else {
			return "" + d;
		}
	}

	/**
	 * Modify a filename f, appending a counter i just before the filetype extension. 
	 */
	public static String addCounterSuffixToFilename(String f, int i)
	{
		return addSuffixToFilename(f, "" + i);
	}

	/**
	 * Modify a filename f, appending a string s just before the filetype extension. 
	 */
	public static String addSuffixToFilename(String f, String s)
	{
		int j = f.lastIndexOf(".");
		if (j != -1) {
			return f.substring(0, j) + s + f.substring(j);
		} else {
			return f + s;
		}
	}

	/**
	 * Format a fraction as a percentage to 1 decimal place.
	 */
	public static String formatPercent1dp(double frac)
	{
		return formatterPercent1dp.format(frac);
	}

	private static DecimalFormat formatterPercent1dp = new DecimalFormat("#0.0%", DecimalFormatSymbols.getInstance(Locale.UK));

	/**
	 * Format a double to 2 decimal places.
	 */
	public static String formatDouble2dp(double d)
	{
		return formatterDouble2dp.format(d);
	}

	private static DecimalFormat formatterDouble2dp = new DecimalFormat("#0.00", DecimalFormatSymbols.getInstance(Locale.UK));

	/**
	 * Format a double, as would be done by printf's %.12g
	 */
	public static String formatDouble(double d)
	{
		// Use UK locale to avoid . being changed to , in some countries.
		// To match C's printf, we have to tweak the Java version,
		// strip trailing zeros after the .
		String result = String.format(Locale.UK, "%.12g", d);
		// if there are only zeros after the . (e.g., .000000), strip them including the . 
		result = result.replaceFirst("\\.0+(e|$)", "$1");
		// handle .xxxx0000
		// we first match .xxx until there are only zeros before the end (or e)
		// as we match reluctantly (using the *?), all trailing zeros are captured
		// by the 0+ part
		result = result.replaceFirst("(\\.[0-9]*?)0+(e|$)", "$1$2");
		return result;
	}

	/**
	 * Format a double, as would be done by printf's %.(prec)g
	 */
	public static String formatDouble(int prec, double d)
	{
		// Use UK locale to avoid . being changed to , in some countries.
		// To match C's printf, we have to tweak the Java version,
		// strip trailing zeros after the .
		String result = String.format(Locale.UK, "%." + prec + "g", d);
		// if there are only zeros after the . (e.g., .000000), strip them including the . 
		result = result.replaceFirst("\\.0+(e|$)", "$1");
		// handle .xxxx0000
		// we first match .xxx until there are only zeros before the end (or e)
		// as we match reluctantly (using the *?), all trailing zeros are captured
		// by the 0+ part
		result = result.replaceFirst("(\\.[0-9]*?)0+(e|$)", "$1$2");
		return result;
	}

	/**
	 * Create a string for a list of objects, with a specified separator,
	 * e.g. ["a","b","c"], "," -&gt; "a,b,c"
	 */
	public static String joinString(List<?> objs, String separator)
	{
		String s = "";
		boolean first = true;
		for (Object obj : objs) {
			if (first) {
				first = false;
			} else {
				s += separator; 
			}
			s += obj.toString();
		}
		return s;
	}
	
	/**
	 * Create a string for an array of objects, with a specified separator,
	 * e.g. ["a","b","c"], "," -&gt; "a,b,c"
	 */
	public static String joinString(Object[] objs, String separator)
	{
		String s = "";
		boolean first = true;
		for (Object obj : objs) {
			if (first) {
				first = false;
			} else {
				s += separator; 
			}
			s += obj.toString();
		}
		return s;
	}
	
	public static void normalizeArray(double[] originalArr)
	{
		double sum = 0;
		int arrSize = originalArr.length;
		for (int i = 0; i < arrSize; i++) {
			sum += originalArr[i];
		}
		if (sum == 0) {
			for (int i = 0; i < arrSize; i++) {
				originalArr[i] = 1.0 / arrSize;
			}
		} else {
			for (int i = 0; i < arrSize; i++) {
				originalArr[i] /= sum;
			}
		}
	}
		
	public static boolean isTargetBlief(double[] belief, BitSet target)
	{
		 double prob=0;
		 for (int i = target.nextSetBit(0); i >= 0; i = target.nextSetBit(i+1)) 
		 {
			 prob+=belief[i];
		 }
		 if(Math.abs(prob-1.0)<1.0e-6)
		 {
			 return true;
		 }
		 return false;
	}
	
	 /** Picks a random item from an array of probabilities,
    normalized and summed as follows:  For example,
    if four probabilities are {0.3, 0.2, 0.1, 0.4}, then
    they should get normalized and summed by the outside owners
    as: {0.3, 0.5, 0.6, 1.0}.  If probabilities.length < checkboundary,
    then a linear search is used, else a binary search is used.    
*/

	public static int pickFromDistribution(final double[] dist,final double prob)
	{
	    if (prob<0.0 || prob>1.0)
	        throw new ArithmeticException("Invalid probability for pickFromDistribution (must be 0.0<=x<=1.0)");
	    if (dist.length==1) // quick 
	        return 0;
	    double[] probabilities = organizeDistribution(dist, false);
        // binary search
        int top = probabilities.length-1;
        int bottom = 0;
        int cur;

        while(top!=bottom)
        {
            cur = (top + bottom) / 2; // integer division

            if (probabilities[cur] > prob)
                if (cur==0 || probabilities[cur-1] <= prob)
                    return exemptZeroes(probabilities,cur);
                else // step down
                    top = cur;
            else if (cur==probabilities.length-1) // oops
                return exemptZeroes(probabilities,cur);
            else if (bottom==cur) // step up
                bottom++;  // (8 + 9)/2 = 8
            else
                bottom = cur;  // (8 + 10) / 2 = 9
        }
        
        return exemptZeroes(probabilities,bottom);  // oops
	    
	}
	
	
    public static int pickFromDistribution(final ArrayList<Entry<Integer, Double>> cdf, final double prob)
    {
    	// binary search
        int top = cdf.size()-1;
        int bottom = 0;
        int cur;
       
        while(top!=bottom)
        {
            cur = (top + bottom) / 2; // integer division

            if (cdf.get(cur).getValue() > prob)
                if (cur==0 || cdf.get(cur-1).getValue() <= prob){
                	return cdf.get(cur).getKey();
                }	                        
                else // step down
                    top = cur;
            else if (cur==cdf.size()-1)
            	return cdf.get(cur).getKey();
            else if (bottom==cur) // step up
                bottom++;  // (8 + 9)/2 = 8
            else
                bottom = cur;  // (8 + 10) / 2 = 9
        }
    	return cdf.get(bottom).getKey();
    }
    
    // allows us to have zero-probability values
    private static final int exemptZeroes(final double[] probabilities, int index)
    {
        //System.out.println(index);
        if (probabilities[index]==0.0) // I need to scan forward because I'm in a left-trail
            // scan forward
            { while(index < probabilities.length-1 && probabilities[index]==0.0) index++; }
        else
            // scan backwards
            { while(index > 0 && probabilities[index]==probabilities[index-1]) index--; }
        return index;
    }
    
    /** Normalizes probabilities, then converts them into continuing
    sums.  This prepares them for being usable in pickFromDistribution.
    If the probabilities are all 0, then selection is uniform, unless allowAllZeros
    is false, in which case an ArithmeticException is thrown.  If any of them are negative,
    or if the distribution is empty, then an ArithmeticException is thrown.
    For example, 
    {0.6, 0.4, 0.2, 0.8} -> {0.3, 0.2, 0.1, 0.4} -> {0.3, 0.5, 0.6, 1.0} */
	public static double[] organizeDistribution(final double[] probabilities, final boolean allowAllZeros)
	{
		double[] dist = probabilities.clone();
	    // first normalize
	    double sum=0.0;
	
	    if (dist.length == 0)
	        throw new ArithmeticException("Distribution has no elements");
	
	    for(int x=0;x<dist.length;x++)
	        {
	        if (dist[x]<0.0)
	            throw new ArithmeticException("Distribution has negative probabilities");
	        sum += dist[x];
	        }
	
	    if (sum==0.0)
	        if (!allowAllZeros)
	            throw new ArithmeticException("Distribution has all zero probabilities");
	        else
	            {
	            for(int x=0;x<dist.length;x++)
	                dist[x] = 1.0;
	            sum = dist.length;
	            }
	    
	    for(int x=0;x<dist.length;x++)
	        dist[x] /= sum;
	
	    // now sum
	    sum=0.0;
	    for(int x=0;x<dist.length;x++)
	        {
	        sum += dist[x];
	        dist[x] = sum;
	        }
	    
	    // now we need to work backwards setting 0 values
	    int x;
	    for(x=dist.length-1; x > 0; x--)
	        if (dist[x]==dist[x-1])  // we're 0.0
	            dist[x] = 1.0;
	        else break; 
	    dist[x] = 1.0;
	
	    return dist;
	}
		
	
	public static int pickRandomSetBit(BitSet bits)
	{
		int numSetBit = bits.cardinality();
		int[] setBits = new int[numSetBit];
		for(int i=bits.nextSetBit(0),k=0; (i>=0)&&(k<numSetBit); i=bits.nextSetBit(i+1), k++)
		{
			setBits[k]=i;
		}
		int rand=(int) Math.floor(numSetBit*Math.random());
		
		return setBits[rand];
	}
	
	public static double[] UnboxDoubleArrayList(ArrayList<Double> input)
	{
		double[] result= new double[input.size()];
		for(int i=0;i<input.size();i++)
		{
			result[i]=input.get(i);
		}		
		return result;
	}
	
	 /**
     * Calculate the factorial of n.
     *
     * @param n the number to calculate the factorial of.
     * @return n! - the factorial of n.
     */
    public static long fact(int n) 
    {
		// Base Case: 
		//    If n <= 1 then n! = 1.
		if (n <= 1) {
		    return 1;
		}
		// Recursive Case:  
		//    If n > 1 then n! = n * (n-1)!
		else {
		    return n * fact(n-1);
		}
    }
    
    public static String bitSet2BinaryString(BitSet bs)
    {
    	StringBuffer sb = new StringBuffer();

        for (int i = 0; i < bs.size(); i++) {
          sb.append((bs.get(i)) ? '1' : '0');
        }

        return sb.toString();
    }
	
	/**
	 * Check for any cycles in an 2D boolean array representing a graph.
	 * Useful for checking for cyclic dependencies in connected definitions.
	 * Returns the lowest index of a node contained in a cycle, if there is one, -1 if not.  
	 * @param matrix Square matrix of connections: {@code matr[i][j] == true} iff
	 * there is a connection from {@code i} to {@code j}.
	 */
	public static int findCycle(boolean matrix[][])
	{
		int n = matrix.length;
		int firstCycle = -1;
		// Go through nodes
		for (int i = 0; i < n; i++) {
			// See if there is a cycle yet
			for (int j = 0; j < n; j++) {
				if (matrix[j][j]) {
					firstCycle = j;
					break;
				}
			}
			// If so, stop
			if (firstCycle != -1)
				break;
			// Extend dependencies
			for (int j = 0; j < n; j++) {
				for (int k = 0; k < n; k++) {
					if (matrix[j][k]) {
						for (int l = 0; l < n; l++) {
							matrix[j][l] |= matrix[k][l];
						}
					}
				}
			}
		}
		return firstCycle;
	}

	/**
	 * Convert a string representing an amount of memory (e.g. 125k, 50m, 4g) to the value in KB.
	 * If the letter prefix is omitted, we assume it is "k" (i.e. KB).
	 */
	public static long convertMemoryStringtoKB(String mem) throws PrismException
	{
		Pattern p = Pattern.compile("([0-9]+)([kmg]?)");
		Matcher m = p.matcher(mem);
		if (!m.matches()) {
			throw new PrismException("Invalid amount of memory \"" + mem + "\"");
		}
		long num;
		try {
			num = Long.parseLong(m.group(1));
		} catch (NumberFormatException e) {
			throw new PrismException("Invalid amount of memory \"" + mem + "\"");
		}
		switch (m.group(2)) {
		case "":
		case "k":
			return num;
		case "m":
			return num * 1024;
		case "g":
			return num * (1024 * 1024);
		default:
			// Shouldn't happen
			throw new PrismException("Invalid amount of memory \"" + mem + "\"");
		}
	}

	/**
	 * Convert a string representing an amount of time (e.g. 5s, 5m, 5h, 5d, 5w) to the value
	 * in seconds.
	 * If the unit is omitted, we assume it is seconds.
	 */
	public static int convertTimeStringtoSeconds(String time) throws PrismException
	{
		Pattern p = Pattern.compile("([0-9]+)([smhdw]?)");
		Matcher m = p.matcher(time);
		if (!m.matches()) {
			throw new PrismException("Invalid time value \"" + time + "\"");
		}
		int value;
		try {
			value = Integer.parseInt(m.group(1));
		} catch (NumberFormatException e) {
			throw new PrismException("Invalid time value \"" + time + "\"");
		}
		switch (m.group(2)) {
		case "":
		case "s":  // seconds
			return value;
		case "m":  // minutes
			return value * 60;
		case "h":  // hours
			return value * (60 * 60);
		case "d":  // days
			return value * (24 * 60 * 60);
		case "w":  // weeks
			return value * (7 * 24 * 60 * 60);
		default:
			// Shouldn't happen
			throw new PrismException("Invalid time value \"" + time + "\"");
		}
	}

	/**
	 * Convert a number of bytes to a string representing the amount of memory (e.g. 125k, 50m, 4g).
	 */
	public static String convertBytesToMemoryString(long bytes) throws PrismException
	{
		String units[] = { "b", "k", "m", "g" };
		for (int i = 3; i > 0; i--) {
			long pow = 1 << (i * 10);
			if (bytes >= pow) {
				return (bytes % pow == 0 ? (bytes / pow) : String.format(Locale.UK, "%.1f", ((double) bytes) / pow)) + units[i];
			}
		}
		return bytes + units[0];
		
		/*for (String s : new String[] { "1g", "1500m", "2g", "1000m", "1024m", "1" }) {
			System.out.println(s + " => " + PrismUtils.convertMemoryStringtoKB(s) * 1024 + " => " + PrismUtils.convertBytesToMemoryString(PrismUtils.convertMemoryStringtoKB(s) * 1024));
		}*/
	}
	
	/**
	 * Utility method to create a new PrintStream for a file, but any errors are converted to PrismExceptions
	 */
	public static PrintStream newPrintStream(String filename) throws PrismException
	{
		try {
			return new PrintStream(filename);
		} catch (FileNotFoundException e) {
			throw new PrismException("File \"" + filename + "\" could not opened for output");
		}
	}
}

//------------------------------------------------------------------------------
