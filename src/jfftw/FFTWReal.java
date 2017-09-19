package jfftw;
import java.io.File;

/**   Class to implement simple one-off Real to Complex and Complex to
 *    RealFourier transforms using fftw as a native library.
 *    <p>
 *    A new plan is created for each FFT. 
 *    This is sub-optimal for a large number of
 *    FFTs, but reasonable for ``one-off''s.
 *    @author Will Hossack, 2008
 *    @version 1.3
 *
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
public class FFTWReal extends FFTW {

    /**      Default constructor, does nothing but load the sharable
     *       library.
     */
    public FFTWReal() {
    }


    /**      Constructor to form a <code>FFTWReal</code> from a 
     *       any <code>FFTW</code>. The PlanFlag is also copied.
     *       @param fft the FFTW 
     *       @since 1.1
     */
    public FFTWReal(FFTW fft) {
	super(fft);
    }

    /**      Constructor to optionally load system wisdom file if available.
     *       @param systemWisdom if true loads system wisdom file.
     */
    public FFTWReal(boolean systemWisdom) {
	super(systemWisdom);
    }

    /**      Constructor to load specified wisdom file by name.
     *       @param wisdomFile
     */
    public FFTWReal(String wisdomFile) {
	super(wisdomFile);
    }


    /**     Method to take a one dimensional Forward FFT of a real data 
     *      array of <code>n</code>  elements with output to specifed array. 
     *      The complex output array contains <code>(n/2 + 1)</code> 
     *      complex pairs.
     *      <p>
     *      This method implements out-of-place transforms only, and the
     *      input real array is not modified.
     *      @param realArray the input real array (not changed)
     *      @param complexArray the output array, 
     *      must be of length <code>2*(n/2 + 1)</code>. If
     *      <code>null</code> OR the wrong length a new array will be created.
     *      @return <code>double[]</code> holding the FFT in 
     *      <code>n/2+1</code> complex pairs.
     */
    public double[] oneDimensionalForward(double realArray[], 
					  double complexArray[] ) {

	//   Basic sanity check, (to prevent JRE crashes)

	if (realArray == null || realArray.length == 0) {
	    throw new NullPointerException("FFTWReal.oneDimensionalForward: called with null or zero length data array.");
	}

	//      Check length of complex array and make new space 
	//      if either null or the wrong length.
	int complexLength = 2*(realArray.length/2 + 1);
	if (complexArray == null || complexArray.length != complexLength) {
	    complexArray = new double[complexLength];
	}
	
	if (loaded) {
	    nativeOneDimensionalForward(realArray,complexArray,getPlanFlag());
	}
	else {
	    errorStream.println("FFTWReal.oneDimensionalForward: library not loaded.");
	}
	return complexArray;
    }



    /**     Method to take a one dimensional Forward FFT of a real data 
     *      array of <code>n</code>  elements with output to specifed array. 
     *      The complex output array contains <code>(n/2 + 1)</code> 
     *      complex pairs.
     *      <p>
     *      This method implements out-of-place transforms only, and the
     *      input real array is not modified.
     *      @param realArray the input real array (not changed)
     *      @return <code>double[]</code> holding the FFT in 
     *      <code>n/2+1</code> complex pairs in a new array.
     */    
    public double[] oneDimensionalForward(double realArray[]) {
	return oneDimensionalForward(realArray,null);
    }


    /**       Private native implement the oneDimensionalForward
     */
    private native void nativeOneDimensionalForward(double in[],
						    double out[],
						    int flag);

    /**     Method to take a one dimensional Backward FFT of
     *      a complex hermition array to give Real output.
     *      Th Complex input is held  as double array with alternative 
     *      real/imaginary parts with a total of
     *       <code>n/2+1</code> Complex pairs.
     *      <p>
     *      The output is a double array of <code>n</code> elements.
     *      @param complexArray array of <code>n/2+1</code> complex parts
     *      @param realArray output array of length n. 
     *      If <code>null</code> or wrong length a new
     *      array will be created.
     *      @return <code>double[]</code>  real array of <code>n</code>
     *      real values
     */
    public double[] oneDimensionalBackward(double complexArray[], 
					   double realArray[]) {
	
	//   Basic sanity check, stop JRE crashes
	if (complexArray == null || complexArray.length == 0) {
	    throw new NullPointerException("FFTWReal.oneDimensionalBackward: called with null or zero length data array.");
	}

	//      Check length of real output array and make 
	//      new one of null or wrong length.
	int realLength = 2*(complexArray.length/2 - 1);
	if (realArray == null || realArray.length != realLength) {
	    realArray  = new double[realLength];
	}

	if (loaded) {
	    nativeOneDimensionalBackward(complexArray,realArray,getPlanFlag());
	}
	else {
	   errorStream.println("FFTWReal.oneDimensionalBackward: library not loaded.");
	} 
	return realArray;
    }


    /**     Method to take a one dimensional Backward FFT of
     *      a complex hermition array to give Real output.
     *      Th Complex input is held  as double array with alternative 
     *      real/imaginary parts with a total of
     *       <code>n/2+1</code> Complex pairs.
     *      <p>
     *      The output is a double array of <code>n</code> elements.
     *      @param complexArray array of <code>n/2+1</code> complex parts     
     *      @return <code>double[]</code>  real array of <code>n</code>
     *      real values which will be in a new array.
     */

    public double[] oneDimensionalBackward(double complexArray[]) {
	return oneDimensionalBackward(complexArray,null);
    }
    

    /**       Private native implement the oneDimensionalForward
     */
    private native void nativeOneDimensionalBackward(double in[],
						    double out[],
						     int flag);


    /**    Method to take the two-dimensional Forward FFT of a real data
     *     held in one dimensional double array. The <code>i,j</code> 
     *     pixel of the real image is located in array 
     *     element </code>j*width + i</code>.
     *     <p>
     *     The Complex FFT is returned in a double array with real/imag 
     *     parts in even/odd elements, with the <code>k,l</code> 
     *     Complex components with:
     *      <p>
     *       Real part in location <code>2(l x wft + k)</code><br>
     *       Imag part in location <code>2(l x wft + k) + 1</code><br>
     *     where <code>wft = width/2 + 1</code>
     *     @param width the image width.
     *     @param height the image height.
     *     @param realArray the image data in <code>double[]</code> 
     *     of length <code>width*height</code>.
     *     @param complexArray the output array of length 
     *     <code>2*height*(width/2 + 1)</code>. If <code>null</code> or
     *     wrong length, then a new array will be created.
     *     @return <code>double[]</code> DFT packed into one dimensional array
     *     @throws <code>IllegalArgumentException</code> 
     *     if <code>width*height</code> does not match the array length. 
     *     
     */
    public double[] twoDimensionalForward(int width, int height, 
					  double realArray[], 
					  double complexArray[]) {

	//            Basic sanity checks to stop JRE crashes 
	if (realArray == null || realArray.length == 0) {
	    throw new NullPointerException(
 "FFTWReal.twoDimensionalForward: called with null or zero length data array.");
	}
	if (realArray.length != width*height) {
	    throw new IllegalArgumentException(
       "FFTWReal.twoDimensionalForward: called with width: " + width +
       " height: " + height + " but array length of: " + realArray.length);
	}

	//           Check the size of the output array and make a 
	//           new one if null or wrong length
	int complexLength = 2*height*(width/2 + 1);
	if (complexArray == null || complexArray.length != complexLength) {
	    complexArray = new double[complexLength];
	}


	if (loaded) {
	    nativeTwoDimensionalForward(width,height,realArray,
					complexArray,getPlanFlag());
	}
	else {
	    errorStream.println(
	       "FFTWReal.twoDimensionalForward: library not loaded.");
	} 
	return complexArray;
    }


    /**    Method to take the two-dimensional Forward FFT of a real data
     *     held in one dimensional double array. The <code>i,j</code> 
     *     pixel of the real image is located in array 
     *     element </code>j*width + i</code>.
     *     <p>
     *     The Complex FFT is returned in a double array with real/imag 
     *     parts in even/odd elements, with the <code>k,l</code> 
     *     Complex components with:
     *      <p>
     *       Real part in location <code>2(l x wft + k)</code><br>
     *       Imag part in location <code>2(l x wft + k) + 1</code><br>
     *     where <code>wft = width/2 + 1</code>
     *     @param width the image width.
     *     @param height the image height.
     *     @param realArray the image data in <code>double[]</code> 
     *     of length <code>width*height</code>.
     *     @return <code>double[]</code> DFT packed into one dimensional array
     *     @throws <code>IllegalArgumentException</code> 
     *     if <code>width*height</code> does not match the array length. 
     *     
     */
    public double[] twoDimensionalForward(int width, int height, 
					  double realArray[]) {
	return twoDimensionalForward(width, height, realArray, null);
    }
    


    /**       Private native implement the twoDimensionalForward
     */
    private native void nativeTwoDimensionalForward(int width, int height, 
						    double realArray[],
						    double complexArray[], 
						    int flag);



    /**    Method to take the two dimensional backward DFT of a hermition 
     *     complex DFT held in a double array with real/imag parts in 
     *     even/odd elements, with the <code>k,l</code> Complex components 
     *     held as follows:     
     *     <p>
     *       Real part located at <code>2(j * wft + i)</code><br>
     *       Imag part located at <code>2(j * wft + i) + 1</code><br>
     *     where <code>wft = width/2 + 1</code>.
     *     <p>
     *     The transformed real data is returned in a one-dimension double array
     *     with pixel <code>i,j</code> located at 
     *     element <code>j*width + i</code>.
     *
     *     @param width the image width
     *     @param height the image height
     *     @param complexArray the Fourier data packed in a 1-d array.
     *     @param realArray the returned real data 
     *     of length <code>width*height</code>. If null or
     *     wrong length, then a new array will be created.
     *     @return <code>double[]</code> DFT packed into one dimensional array
     *     @throws <code>IllegalArgumentException</code> if 
     *     <code>width*height</code> is not consistent with array length.
     */
    public double[] twoDimensionalBackward(int width, int height, 
					   double complexArray[], 
					   double realArray[]) {

	//            Basic sanity checks....  
	if (complexArray == null || complexArray.length == 0) {
	    throw new NullPointerException(
"FFTWReal.twoDimensionalBackward: called with null or zero length data array.");
	}
	if (complexArray.length != 2*height*(width/2 + 1)) {
	    throw new IllegalArgumentException(
	"FFTWReal.twoDimensionalBackward: called with width: " + width +
	" height: " + height + " but array length of: " + complexArray.length);
	}

	//       Check that supplied output array is of correct 
	//       length and make a new id needed
	int realLength = width*height;
	if ( realArray == null || realArray.length != realLength) { 
	    realArray = new double[realLength];
	}

	if (loaded) {
	    nativeTwoDimensionalBackward(width,height,complexArray,
					 realArray,getPlanFlag());
	}
	else {
	    errorStream.println(
	     "FFTWReal.twoDimensionalBackward: library not loaded.");
	} 
	return realArray;
    }

    /**    Method to take the two dimensional backward DFT of a hermition 
     *     complex DFT held in a double array with real/imag parts in 
     *     even/odd elements, with the <code>k,l</code> Complex components 
     *     held as follows:     
     *     <p>
     *       Real part located at <code>2(j * wft + i)</code><br>
     *       Imag part located at <code>2(j * wft + i) + 1</code><br>
     *     where <code>wft = width/2 + 1</code>.
     *     <p>
     *     The transformed real data is returned in a one-dimension double array
     *     with pixel <code>i,j</code> located at 
     *     element <code>j*width + i</code>.
     *
     *     @param width the image width
     *     @param height the image height
     *     @param complexArray the Fourier data packed in a 1-d array.
     *     @return <code>double[]</code> DFT packed into one dimensional array
     *     @throws <code>IllegalArgumentException</code> if 
     *     <code>width*height</code> is not consistent with array length.
     */

    public double[] twoDimensionalBackward(int width, int height, 
					   double complexArray[]) {
	return twoDimensionalBackward(width, height, complexArray, null);
    }
    
    /**       Private native implement the twoDimensionalBackward
     */
    private native void nativeTwoDimensionalBackward(int width, int height, 
						     double complexArray[],
						     double realArray[], 
						     int flag);



    /**    Method to take the three dimensional Forward DFT of a real cube
     *     held in one dimensional double array.
     *     the <code>i,j,k</code> pixel of the real cube is located in 
     *     array element <code>k * width * height + j * width + i</code>
     *     <p>
     *     The Complex DFT is returned in a double array with real/imag 
     *     parts in even/odd elements, with the i,j Complex components in
:    *     <p>
     *       Real part located at 
     *       <code>2(k * wft * height + j * wft + i)</code><br>
     *       Imag part located at 
     *       <code>2(k * wft * height + j * wft + i)</code> + 1<br>
     *       where <code>wft = width/2 + 1</code>.
     *
     *     @param width the cube width
     *     @param height the cube height
     *     @param depth the cube depth
     *     @param realArray the cube data in array of length 
     *     <code>width*height*depth</code>
     *     @param complexArray the output array, must be of length 
     *     <code>2*depth*height*(width/2 + 1)</code>. If null
     *     or wrong length, a new array will be allocated.
     *     @return <code>double[]</code> DFT packed into one dimensional array
     *     @throws <code>IllegalArgumentException</code> if 
     *     <code>width*height*depth</code> is not equal the array length.
     *     
     */
    public double[] threeDimensionalForward(int width, int height, int depth, 
					    double realArray[], 
					    double complexArray[] ) {

		//            Basic sanity checks.... 
	if (realArray == null || realArray.length == 0) {
	    throw new NullPointerException(
"FFTWReal.threeDimensionalForward: called with null or zero length data array.");
	}
	if (realArray.length != width*height*depth) {
	    throw new IllegalArgumentException(
    "FFTWReal.threeDimensionalForward: called with width: " + width +
    " height: " + height + " depth : " + depth + 
    " but array length of: " + realArray.length);
	}

	//      Check that complex output array is of correct length, 
	//      else make new one
	int complexLength = 2*depth*height*(width/2 + 1); 
	if ( complexArray == null || complexArray.length != complexLength ) {
	    complexArray = new double[complexLength];       // Make output array
	}


	if (loaded) {
	    nativeThreeDimensionalForward(width,height,depth,realArray,
					  complexArray,getPlanFlag());
	}
	else {
	    errorStream.println(
            "FFTWReal.threeDimensionalForward: library not loaded.");
	} 
	    
	return complexArray;
    }



    /**    Method to take the three dimensional Forward DFT of a real cube
     *     held in one dimensional double array.
     *     the <code>i,j,k</code> pixel of the real cube is located in 
     *     array element <code>k * width * height + j * width + i</code>
     *     <p>
     *     The Complex DFT is returned in a double array with real/imag 
     *     parts in even/odd elements, with the i,j Complex components in
:    *     <p>
     *       Real part located at 
     *       <code>2(k * wft * height + j * wft + i)</code><br>
     *       Imag part located at 
     *       <code>2(k * wft * height + j * wft + i)</code> + 1<br>
     *       where <code>wft = width/2 + 1</code>.
     *
     *     @param width the cube width
     *     @param height the cube height
     *     @param depth the cube depth
     *     @param realArray the cube data in array of length 
     *     <code>width*height*depth</code>
     *     @return <code>double[]</code> DFT packed into one dimensional array
     *     @throws <code>IllegalArgumentException</code> if 
     *     <code>width*height*depth</code> is not equal the array length.
     *     
     */
    public double[] threeDimensionalForward(int width, int height, int depth,  
					    double realArray[]) {
	return threeDimensionalForward(width, height, depth,  realArray, null);
    }



    /**   Private native method to implement the three-dimensional FFT
     */
    private native void nativeThreeDimensionalForward(int width, int height, 
						      int depth, 
						      double realArray[], 
						      double complexArray[],
						      int flag);


    /**    Method to take the three dimensional backward DFT of a hermition 
     *     complex DFT held in a double array with real/imag parts in 
     *     even/odd elements, with the i,j Complex components in
:    *     <p>
     *       Real part located at  
     *       <code>2(k x wft x height + j x wft + i)</code><br>
     *       Imag part located at  
     *       <code>2(k x wft x height + j x wft + i) + 1</code><br>
     *       where <code>wft = width/2 + 1</code>.
     *     <p>
     *     The transformed real data is returned in a one-dimensional 
     *     double array
     *     with pixel i,j located at element 
     *     <code>k x width x height + j x width + i</code>.
     *
     *     @param width the cube width
     *     @param height the cube height
     *     @param depth the cube depth
     *     @param complexArray Fourier data packed in a 1-d array.
     *     @param realArray real output array, must be of 
     *     length <code>width*height*depth</code>. If null
     *     or wrong length, then a new array will be allocated.
     *     @return double[] DFT packed into one dimensional array
     *     @throws <code>IllegalArgumentException </code> 
     *     if width, height, depth is not consistent with array length.
     */    
    public double[] threeDimensionalBackward(int width, int height, int depth, 
					     double complexArray[], 
					     double realArray[]) {

	//            Basic sanity checks.... 
	if (complexArray == null || complexArray.length == 0) {
	    throw new NullPointerException(
"FFTWReal.threeDimensionalBackward: called with null or zero length data array.");
	}
	if (complexArray.length != 2*depth*height*(width/2 + 1)) {
	    throw new IllegalArgumentException(
    "FFTWReal.threeDimensionalBackword: called with width: " + width +
    " height: " + height + " depth : " + depth + 
    " but array length of: " + complexArray.length);
	}



	//           Check that real output array is of currect length,
	//           else allocate new one
	int realLength = depth*height*width;
	if (realArray == null || realArray.length != realLength) {
	    realArray = new double[realLength];      // Make output array
	}

	if (loaded) {
	    nativeThreeDimensionalBackward(width,height,depth,
					   complexArray,realArray,
					   getPlanFlag());
	}
	else {
	    errorStream.println(
	"FFTWReal.threeDimensionalBackward: library not loaded.");
	} 
	
	return realArray;
    }


    /**    Method to take the three dimensional backward DFT of a hermition 
     *     complex DFT held in a double array with real/imag parts in 
     *     even/odd elements, with the i,j Complex components in
:    *     <p>
     *       Real part located at  
     *       <code>2(k x wft x height + j x wft + i)</code><br>
     *       Imag part located at  
     *       <code>2(k x wft x height + j x wft + i) + 1</code><br>
     *       where <code>wft = width/2 + 1</code>.
     *     <p>
     *     The transformed real data is returned in a one-dimensional 
     *     double array
     *     with pixel i,j located at element 
     *     <code>k x width x height + j x width + i</code>.
     *
     *     @param width the cube width
     *     @param height the cube height
     *     @param depth the cube depth
     *     @param complexArray Fourier data packed in a 1-d array.
     *     @return double[] DFT packed into one dimensional array
     *     @throws <code>IllegalArgumentException </code> 
     *     if width, height, depth is not consistent with array length.
     */    

    public double[] threeDimensionalBackward(int width, int height, int depth, 
					     double complexArray[] ) {
	return  threeDimensionalBackward(width, height, depth, 
					 complexArray, null);
    }
    

    /**   Private native method to implement the three-dimensional FFT
     */
    private native void nativeThreeDimensionalBackward(int width, int height, 
						       int depth, 
						       double complexArray[], 
						       double realArray[],
						       int flag);

}
