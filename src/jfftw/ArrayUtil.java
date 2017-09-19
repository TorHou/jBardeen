package jfftw;

//import Complex;

/**     Class containing some useful static method to manipulate
 *      double arrays containing complex data used in jfftw.
 *      All array access is by array index with no setter/getter
 *      overhead. This is not very OOP, but vastly improves
 *      efficiency.
 * 
 *      @author Will Hossack, 2008
 *
 *    This file is part of jfftw.
 *
 *    Jfftw is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    Jfftw is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with Jfftw.  If not, see <http://www.gnu.org/licenses/>.
 *
 * 
 */
public class ArrayUtil extends Object {

    /**   Static method to take two double arrays normally holding
     *    real and imaginary parts and interleave them to a single
     *    array. If imaginary array is null or length zero it is
     *    taken as zero.
     *    @param real the real array
     *    @param imag the imaginary array (can be null, or length zero).
     *    @return <code>double[]</code> interleaved array.
     *    
     */
    public static double[] interleave(double[] real, double[] imag) {
	boolean isImaginary;
	if (imag == null || imag.length == 0) {
	    isImaginary = false;
	}
	else {
	    isImaginary = true;
	}

	double[] data = new double[2*real.length]; // Array of double length
	for(int i = 0; i < real.length; i++) {
	    int j = 2*i;
	    data[j] = real[i];                     // Copy real
	    if (isImaginary) data[j+1] = imag[i];  // Copy imag if exists
	}
	return data;
    }

    /**   Static method to take split array normally holding
     *    real and imaginary parts and interleave them to a single
     *    array. 
     *    @param data the split array
     *    @return <code>double[]</code> interleaved array.
     *    
     */
    public static double[] interleave(double[][] data) {
	return interleave(data[0],data[1]);
    }


    /**   Static method to split an interleaved array into two split
     *    array, returned as a double[2][] array with [0][i] and
     *    [1][i] holding the real / imaginary parts.
     *    @param data the interleaved array.
     *    @return <code>double[][]</code> of size [2][data.length/2]
     */
    public static double[][] split(double[] data) {
	int size = data.length/2;
	double[][] sp = new double[2][size];

	for(int i = 0; i < size; i++) {
	    int j = 2*i;
	    sp[0][i] = data[j];
	    sp[1][i] = data[j + 1];
	}
	return sp;
    }

    /**   Static method to form a split array from a real and
     *    imaginary array. If imaginary array is null or of length
     *    zero, it is assumed to be zero.
     *    <p>
     *    Note the data is not copied to new space.
     *    @param real the real array
     *    @param imag the imaginary array
     *    @return <code>double[][]</code> array in split format.
     */
    public static double[][] split(double[] real, double[] imag) {
	double[][] data = new double[2][];
	data[0] = real;
	if (imag == null || imag.length == 0) {
	    data[1] = new double[real.length];
	}
	else {
	    data[1] = imag;
	}
	return data;
    }

    /**   Static method to get the ith Complex element from an
     *    interleaved array where real part is in 2*i and imaginary
     *    in 2*i + 1.
     *    @param data the interleaved array
     *    @param i the element index
     */
    public static Complex getComplex(double[] data,int i) {
	int j = 2*i;
	return new Complex(data[j],data[j + 1]);
    }

    /**    Static to set the ith component of a interleaved complex
     *     array with a Complex.
     *     @param data the interleaved array
     *     @param i the element index
     *     @param z the value
     */
    public static void setComplex(double[] data, int i, Complex z) {
	int j = 2*i;
	data[j] = z.r;
	data[j + 1] = z.i;
    }


    /**    Static to set the ith component of a interleaved complex
     *     array with a two doubles
     *     @param data the interleaved array
     *     @param i the element index
     *     @param a the real part 
     *     @param b the imaginary part.
     */
    public static void setComplex(double[] data, int i, double a, double b) {
	int j = 2*i;
	data[j] = a;
	data[j + 1] = b;
    }


    
    /**   Static method to get the ith Complex element from a
     *    split array consisting of real and imaginary parts.
     *    @param real the real array.
     *    @param imag the imaginary array.
     *    @param i the element index.
     */
    public static Complex getComplex(double[] real, double[] imag, int i) {
	return new Complex(real[i],imag[i]);
    }


    /**    Static to set the ith component of split array 
     *     array with a Complex.
     *     @param real the real array
     *     @param imag the imaginary array
     *     @param i the element index
     *     @param z the value
     */
    public static void setComplex(double[] real, double[] imag, int i,
				  Complex z) {
	real[i] = z.r;
	imag[i] = z.i;
    }

    
    /**    Static to set the ith component of split array 
     *     array with two doubles
     *     @param real the real array
     *     @param imag the imaginary array
     *     @param i the element index
     *     @param a the real part
     *     @param b the imaginary part
     */
    public static void setComplex(double[] real, double[] imag, int i,
				  double a, double b) {
	real[i] = a;
	imag[i] = b;
    }


    /**   Static method to get the ith Complex element from a
     *    split array consisting of real and imaginary parts in
     *    a [2][] array
     *    @param split the [2][] array
     *    @param i the element index.
     */
    public static Complex getComplex(double[][] split, int i) {
	return getComplex(split[0],split[1],i);
    }


    /**    Static to set the ith component of split array 
     *     array with a Complex.
     *     @param split the split array of real and imaginary parts
     *     @param i the element index
     *     @param z the value
     */
    public static void setComplex(double[][] split, int i, Complex z) {
	setComplex(split[0],split[1],i,z);
    }

    /**    Static to set the ith component of split array 
     *     array with a Complex.
     *     @param split the split array of real and imaginary parts
     *     @param i the element index
     *     @param a the real part
     *     @param b the imaginary part.
     */
    public static void setComplex(double[][] split, int i, 
				  double a, double b) {
	setComplex(split[0],split[1],i,a,b);
    }


    /**     Conjugate specified element of an interleaved array
     *      @param i the index
     */
    public static void conjugate(double[] data, int i){
	int ii = 2*i + 1;                  // Imaginary part
	data[ii] = -data[ii];              
    }
  

    /**    Get the modulus square of specified element of interleaved array
     *     @param i the index
     */
    public static double modulusSqr(double[] data, int i) {
	int j = 2*i;
	double x = data[j];
	double y = data[j + 1];
	return x*x + y*y;
    }

    /**    Get the modulus of specified element of interleaved array
     *     @param i the index
     */
    public static double modulus(double[] data, int i) {
	return Math.sqrt(modulusSqr(data,i));
    }

    /**     Multiply specified element of interleaved array by a double
     *      @param data the interleaved array
     *      @param i the index
     *      @param a the double.
     */
    public static void mult(double[] data, int i, double a) {
	int j = 2*i;
	data[j] *= a;
	data[j+1] *= a;
    }

    /**     Multiply specified element of interleaved array by a complex
     *      specified as two doubles.
     *      @param data the interleaved array
     *      @param i the index
     *      @param a real part
     *      @param b imaginary part
     */
    public static void mult(double[] data, int i, double a, double b) {
	int j = 2*i;
	int jj = j + 1;
	double x = data[j];
	double y = data[jj];
	data[j] =  a*x - b*y;
	data[jj] = b*x + a*y;
    }

    /**     Multiply specified element of interleaved array by a complex.
     *      @param data the interleaved array
     *      @param i the index
     *      @param c the Complex
     */
    public static void mult(double[] data, int i, Complex c) {
	mult(data,i,c.r,c.i);
    }

    /**   Method to multiply an interleaved array by a double.
     *    @param data the interleaved array
     *    @param a the multiplier.
     */
    public static void mult(double[] data, double a) {
	for(int i = 0; i < data.length; i++) {
	    data[i] *= a;
	}
    }

    /**   Method to multiply a split array by a double.
     *    @param split the split array
     *    @param a the multiplier.
     */
    public static void mult(double[][] split, double a) {
	mult(split[0],a);       // deal with real part
	mult(split[1],a);       // deal with imaginary part
    }


    /**   Method to multiply a interleaved array by a complex specified as
     *    two doubles.
     *    @param data the interleaved array
     *    @param a real part of multiplier
     *    @param b imaginary part of multiplier
     */    
    public static void mult(double[] data, double a, double b) {
	for(int i = 0; i < data.length; i += 2) {    // Interleaved array
	    int ii = i + 1;
	    double x = data[i];
	    double y = data[ii];
	    data[i] =  a*x - b*y;                     // Real part
	    data[ii] = b*x + a*y;                     // Imag part
	}
    }
    
    /**   Method to multiply an interleaved array by a complex.
     *    @param data the interleaved array
     *    @param c the Complex multiplier
     */
    public static void mult(double[] data, Complex c) {
	mult(data,c.r,c.i);
    }


    /**   Method to multiply a split array by a complex specified as
     *    two doubles.
     *    @param split the split array
     *    @param a real part of multiplier
     *    @param b imaginary part of multiplier
     *    @throws <code>ArrayIndexOutOfBoundsException</code>  if
     *    real and imaginary arrays of different lengths.
     */    
    public static void mult(double[][] split, double a, double b) {
	if (split[0].length != split[1].length ) {
	     throw new ArrayIndexOutOfBoundsException(
          "ArrayUtil.mult: Real and imaginary arrays of different lengths");
	}

	for(int i = 0; i < split[0].length; i++) {
	    double x = split[0][i];
	    double y = split[1][i];
	    split[0][i] = a*x - b*y;             // Real Part
	    split[1][i] = b*x + a*y;             // Imag part
	}
    }

    /**   Method to multiply a split array by a complex.
     *    @param split the split array
     *    @param c the Complex multiplier
     */
    public static void mult(double[][] split, Complex c) {
	mult(split,c.r,c.i);
    }


    /**   Method to multiply two complex interleaved data arrays
     *    placing the result in the first. Both must be the same size.
     *    @param first the first interleaved array
     *    @param second the second interleaved array (not changed)
     *    @throws <code>ArrayIndexOutOfBoundsException</code> if
     *     the two arrays are not the same length.
     */
    public static void mult(double[] first, double [] second) {
	if (first.length != second.length) {
	    throw new ArrayIndexOutOfBoundsException(
	      "ArrayUtil.mult: Arrays of different lengths passed");
	}
	
	for(int i = 0; i < first.length; i += 2) {
	    int ii = i + 1;
	    double x = first[i];           // First real 
	    double y = first[ii];          // First imag
	    double a = second[i];          // Second real
	    double b = second[ii];         // Second imag

	    first[i] =  a*x - b*y;         // Real part
	    first[ii] = b*x + a*y;         // Imag part
	}
    }



    /**   Method to multiply a interleaves complex arrays with the
     *    complex conjugate of an interleaved array, 
     *    placing the result in the first. Both must be the same size.
     *    @param first the first interleaved array
     *    @param second the second interleaved array (not changed)
     *    @throws <code>ArrayIndexOutOfBoundsException</code> if
     *     the two arrays are not the same length.
     */
    public static void multConjugate(double[] first, double [] second) {
	if (first.length != second.length) {
	    throw new ArrayIndexOutOfBoundsException(
	      "ArrayUtil.multConjugate: Arrays of different lengths passed");
	}
	
	for(int i = 0; i < first.length; i += 2) {
	    int ii = i + 1;
	    double x = first[i];           // First real 
	    double y = first[ii];          // First imag
	    double a = second[i];          // Second real
	    double b = second[ii];         // Second imag

	    first[i] =  a*x + b*y;         // Conjugate multiply
	    first[ii] = a*y - b*x;
	}
    }

    /**    Method to multiply two split arrays putting the result in the
     *     first. Both must be of the same length.
     *     @param first the first split array
     *     @param second the second split array.
     *     @throws <code>ArrayIndexOutOfBoundsException</code> if
     *     the two arrays are not the same length.
     *
     */     
    public static void mult(double[][] first, double [][]second) {
	
	//           Do sanity check that all arrays are the same length
	if (first[0].length != first[1].length || second[0].length != second[1].length ||
	    first[0].length != second[0].length ) {
	    throw new ArrayIndexOutOfBoundsException(
	      "ArrayUtil.multConjugate: Arrays of different lengths passed");
	}

	
	for(int i = 0; i < first[0].length; i++) {
	    double x = first[0][i];           // First real
	    double y = first[1][i];           // First Imag
	    double a = second[0][i];          // Second real
	    double b = second[1][i];          // Second imag

	    first[0][i] = a*x - b*y;
	    first[1][i] = b*x + a*y;
	}
    }


    /**    Method to multiply first array by complex conjugate of the second.
     *     Both must be of the same length and the second array is not altered.
     *     @param first the first split array
     *     @param second the second split array.
     *     @throws <code>ArrayIndexOutOfBoundsException</code> if
     *     the two arrays are not the same length.
     *
     */     
    public static void multConjugate(double[][] first, double [][]second) {
	
	//           Do sanity check that all arrays are the same length
	if (first[0].length != first[1].length || second[0].length != second[1].length ||
	    first[0].length != second[0].length ) {
	    throw new ArrayIndexOutOfBoundsException(
	      "ArrayUtil.multConjugate: Arrays of different lengths passed");
	}

	
	for(int i = 0; i < first[0].length; i++) {
	    double x = first[0][i];           // First real
	    double y = first[1][i];           // First Imag
	    double a = second[0][i];          // Second real
	    double b = second[1][i];          // Second imag

	    first[0][i] = a*x + b*y;          // Conjugate multiply
	    first[1][i] = a*y - b*x;
	}
    }


    /**   Method to take form the conjugate of an interleaved array
     *    by taking the negative of imaginary parts.
     *    @param data the interleaved array.
     */
    public static void conjugate(double[] data) {
	for(int i = 1; i < data.length; i += 2) {
	    data[i] = -data[i];
	}
    }

    /**   Method to take form the conjugate of a split array 
     *    by taking the negative of the imaginary parts.
     *    @param split the split array.
     */
    public static void conjugate(double[][] split) {
	double[] d = split[1];        // Use 0ne-dimensional access
	for(int i = 0; i < d.length; i++) {
	    d[i] = -d[i];
	}
    }


    /**   Method to get the modulus squared of the largest
     *    elements of an interleaved Complex array.
     *    @param data the interleaved Complex array
     *    @return <code>double</code> maximum modulus squared 
     *    
     */
    public static double maxModSqr(double[] data ) {
	double x, y, m, max = 0.0;
	for(int i = 0; i < data.length; i += 2) {
	    x = data[i];                      // Real part
	    y = data[i+1];                    // Imaginary part
	    m = x*x + y*y;
	    if (m > max) max = m;
	}
	return max;
    }

    /**   Method to get the modulus squared of the largest
     *    elements of a split Complex array with real and
     *    imaginary parts held in separate arrays.
     *    @param real the real array
     *    @param imag the imaginary array
     *    @return <code>double</code> maximum modulus squared
     */ 
    public static double maxModSqr(double[] real, double[] imag) {
	double x, y, m, max = 0.0;
	for(int i = 0; i < real.length; i++) {
	    x = real[i];
	    y = imag[i];
	    m = x*x + y*y;
	    if (m > max) max = m;
	}
	return max;
    }

    /**   Method to get the modulus squared of the largest
     *    elements of a split Complex array
     *    @param split the split array
     *    @return <code>double</code> maximum modulus squared
     */
    public static double maxModSqr(double[][] split) {
	return maxModSqr(split[0],split[1]);
    }


    /**    Method to get power (sum of modulus squared) in an
     *     Complex interleaved array.
     *     @param data the interleaved Complex array
     *     @return <code>double</code> the power.
     */
    public static double power(double[] data) {
	double x,y, p = 0;
	for(int i = 0; i < data.length; i += 2) {
	    x = data[i];
	    y = data[i + 1];
	    p += x*x + y*y;
	}
	return p;
    }

    /**    Method to get power (sum of modulus squared) in an
     *     Complex split array held in real and imaginary
     *     arrays.
     *     @param real the real array
     *     @param imag the imaginary array
     *     @return <code>double</code> the power
     */
    public double power(double[] real, double[] imag) {
	double x,y, p = 0;
	for(int i = 0; i < real.length; i++) {
	    x = real[i];
	    y = imag[i];
	    p += x*x + y*y;
	}
	return p;
    }

    /**    Method to get power (sum of modulus squared) in an
     *     Complex split array.
     *     @param split the split array
     *     return <code>double</code> the power
     */
    public double power(double[][] split) {
	return power(split[0],split[1]);
    }
    
}
    
    

	
