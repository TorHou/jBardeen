package jfftw;
import java.io.File;

/**   Class to implement simple one-off Complex to Complex Fourier
 *    transforms using jfftw as a native library.
 *    <p>
 *    A new plan is created for each FFT. 
 *    This is sub-optimal for a large number of  FFTs, but reasonable 
 *    for ``one-off''.
 *    @author Will Hossack, 2008
 *    @version 1.2
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
public class FFTWComplex extends FFTW {

    /**      Default constructor which load sharable library and sets defaults.
     */
    public FFTWComplex() {
    }
    
    /**      Constructor to form a FFTWComplex from a FFTW, 
     *       The PlanFlag is also copied.
     *       @param fft the FFTW 
     *       @since 1.1
     */
    public FFTWComplex(FFTW fft) {
	super(fft);
    }
    

    /**      Constructor to optionally load system wisdom file
     *       @param systemWisdom if true loads system wisdom file
     */
    public FFTWComplex(boolean systemWisdom) {
	super(systemWisdom);
    }

    /**      Constructor to load specified wisdom file by name
     *       @param wisdomFile
     */
    public FFTWComplex(String wisdomFile) {
	super(wisdomFile);
    }

    /**      Method to take one-dimensional Complex FFT with the data supplied
     *       in a one-dimensional interleaved <code>double</code> array.
     *       @param in the input data
     *       @param out the output data (may be the same as in)
     *       @param dirn direction +1 for forward, -1 to backward
     *       @return <code>double[]</code> the output array being <code>out</code>
     *       @throws NullPointerException of null or zero length arrays send.
     *       @throws IllegalArgumentException of in and out arrays of different length
     */
    public double[] oneDimensional(double in[], double out[], int dirn) {
	
	//            Basic sanity check
	if (in == null || in.length == 0 || out == null || out.length == 0) {
	    throw new NullPointerException(
  "FFTComplex.oneDimensional: called with null or zero length data array.");
	}
	if (in.length != out.length) {
	    throw new IllegalArgumentException(
  "FFTComplex.oneDimensional: different length input and output arrays");
	}
	
	if (loaded) {
	    nativeOneDimensional(in,out,dirn,getPlanFlag());    
	}
	else {
	    errorStream.println(
     "FFTWComplex.oneDimensional: library not loaded.");
	}
	return out;     
	
    }


    /**      Method to take one-dimensional Complex FFT with the data supplied
     *       in a one-dimensional interleaved <code>double</code> array.
     *       @param in the data to be transformed
     *       @param dirn direction +1 for forward, -1 to backward
     *       @param overwrite if true data is overwritten with FFT, else a
     *              new array is returned.
     *       @return <code>double[]</code> the FFTed array 
     *          (if overwrite = true, this will overwrite the input array).
     *       
     */
    public double[] oneDimensional(double in[], int dirn, boolean overwrite) {

	if (overwrite) {           // In-place transform
	    return oneDimensional(in,in,dirn);
	}
	else {                      // Make new array of same length
	    return oneDimensional(in,new double[in.length],dirn);
	}
    }


    /**       nativeOneDimensional method to do the work
     */
    private native void nativeOneDimensional(double in[], double out[], 
					     int dirn, int flag);


    /**   Method to take a Complex one-dimensional FFT with the real and imaging
     *    data split format in a two-dimensional array of size [2][length] 
     *    with the the i th component with real in the [0][i] 
     *    element and the imaginary element in the [1][i] element.
     *    @param in two-dimensional array, real in [0][i] and imag in [0][i]
     *    @param out two-dimensional array to return the data.
     *    @param dirn forward or inverse transform
     *    @return <code>double[][]</code> transformed data, being out
     *    @throws <code>IllegalArgumentException</code> if real and imaginary 
     *     array of different length.
     */
    public double[][] oneDimensional(double[][]in, double[][]out, int dirn) {
	//           Basic  sanity check
	if (in == null || in.length != 2 || 
	    in[0].length == 0 || in[1].length == 0) {
	    throw new NullPointerException(
   "FFTComplex.oneDimensional: called with null or zero length in array(s).");
	}

	if (out == null || out.length != 2 || 
	    out[0].length == 0 || out[1].length == 0) {
	    throw new NullPointerException(
   "FFTComplex.oneDimensional: called with null or zero length out array(s).");
	}
	
	if ( in[0].length != in[1].length || out[0].length != out[1].length ||
	     in[0].length != out[0].length ) {
	    throw new IllegalArgumentException(
   "FFTWComplex.oneDimensional: called with real array of length: " 
   + in[0].length + " and imaginary array of length: " + in[1].length);
	}

	if (loaded) {
	    //          If forward send data in normal order
	    if (dirn == FORWARD) {
		nativeOneDimensionalSplit(in[0],in[1],out[0],out[1],
					  getPlanFlag());
	    }
	    else {   // For inverse swap real/imag parts
		nativeOneDimensionalSplit(in[1],in[0],out[1],out[0],
					  getPlanFlag());
	    }
	}
	else {
	    errorStream.println(
		"FFTWComplex.oneDimensional: library not loaded.");
	}
	    
	return out;
    }

    /**   Method to take a Complex one-dimensional FFT with the real and imaging
     *    data split format in a two-dimensional array of size [2][length] 
     *    with the the i th component with real in the [0][i] 
     *    element and the imaginary element in the [1][i] element.
     *    @param in two-dimensional array, real in [0][i] and imag in [0][i]
     *    @param dirn forward or inverse transform
     *    @param overwrite if true transform will overwrite in,
     *            else a new array will be created.
     *    @return <code>double[][]</code> transformed data.
     */
    public double[][] oneDimensional(double[][]in, int dirn, 
				     boolean overwrite) {

	if (overwrite) {
	    return oneDimensional(in,in,dirn);
	}
	else {
	    return oneDimensional(in,new double[2][in[0].length],
				  dirn);
	}
    }



    /**     native OneDimensionalSplit to do the work
     */
    private native void nativeOneDimensionalSplit(double realin[], 
						  double imagin[],
						  double realout[], 
						  double imagout[],
						  int flag);
	

    /**      Method to take two-dimensional Complex FFT with the data supplied
     *       in a one-dimensional double array with real parts in even elements,
     *       and imaginary in the odd.Element i,j is located at 
     *       
     *       Real part     2*(j*width + i)
     *       Imag part     2*(j*width + i) + 1
     *
     *       @param width the width of the data
     *       @param height the height of the data
     *       @param in the input data array to be transformed
     *       @param out the output data array, may be the same as <code>in</code>
     *       @param dirn direction +1 for forward, -1 to backward
     *       @return <code>double[]</code> the FFTed array
     *       @throws <code>NullPointerException</code> if null or zero length data send.
     *       @throws <code>IllegalArgumentException</code> if width and height 
     *       does not match the both array lengths.
     *
     */
    

    public double[] twoDimensional(int width, int height, double in[],
				   double out[], int dirn) {
	//           Basic sanity checks
	if (in == null || in.length ==0 || out == null || out.length== 0) {
	    throw new NullPointerException(
"FFTComplex.twoDimensional: called with null or zero length data array.");
	}
	if (in.length != 2*width*height || out.length != in.length ) {
            throw new IllegalArgumentException(
	 "FFTWComplex.twoDimensional: called with unmatched array lengths.");
	}
	
	if (loaded) {
	    nativeTwoDimensional(width,height,in,out,dirn,getPlanFlag());
	}
	else {
	    errorStream.println("FFTWComplex.twoDimensional: library not loaded.");
	}
	    
	return out;
    }
	


    /**      Method to take two-dimensional Complex FFT with the data supplied
     *       in a one-dimensional double array with real parts in even elements,
     *       and imaginary in the odd Element i,j is located at 
     *       
     *       Real part     2*(j*width + i)
     *       Imag part     2*(j*width + i) + 1
     *
     *       @param width the width of the image data
     *       @param height the height of the image data
     *       @param in the data to be transformed
     *       @param dirn direction +1 for forward, -1 to backward
     *       @param overwrite if true data is overwritten with FFT, else a
     *              new array is returned.
     *       @return <code>double[]</code> the FFTed array 
     *       (if overwrite = true, this will
     *       be the same at the input data array.
     *
     */
    
    public double[] twoDimensional(int width, int height, double in[], 
				   int dirn, boolean overwrite) {

	if (overwrite) {                                     // In place transform
	    return twoDimensional(width,height,in,in,dirn);
	}
	else {      
	    return twoDimensional(width,height,in,new double[in.length],dirn);
	}
    }
    
    /**       nativeTwoDimensional method to do the work
     */
    private native void nativeTwoDimensional(int width, int height, 
					     double in[], double out[], 
					     int dirn, int flag);



    /**   Method to take a Complex two-dimensional FFT with the real and imaging
     *    data split format in a two-dimensional array of size [2][length] 
     *    with the the i th component with real in the [0][i] element and the 
     *    imaginary element in the [1][i] element.
     *    @param width width of the data
     *    @param height height of the data
     *    @param in two-dimensional array, real in [0][i] and imag in [0][i]
     *    @param out two-dimensional array for the output.
     *    @param dirn forward or inverse transform
     *    @return <code>double[][]</code> transformed data
     *    @throws <code>IllegalArgumentException</code> if real and imaginary 
     *    array of different length, or
     *    width*height does not match the array lengths.
     */
    public double[][] twoDimensional(int width, int height, double[][] in,
				     double[][] out,int dirn) {
	//           Basic  sanity check
	if (in == null || in.length != 2 || in[0].length == 0 || 
	    in[1].length == 0) {
	    throw new NullPointerException(
"FFTComplex.twoDimensional: called with null or zero length in array(s).");
	}
	
	if (out == null || out.length != 2 || 
	    out[0].length == 0 || out[1].length == 0) {
	    throw new NullPointerException(
"FFTComplex.twoDimensional: called with null or zero length in array(s).");
	}

	if ( in[0].length != in[1].length || in[0].length != width*height ||
	     out[0].length != out[1].length || in[0].length != out[0].length) {
	    throw new IllegalArgumentException(
	  "FFTWComplex.twoDimensional: called with width : " + width +
	  " height: " + height + " but array length(s) of " +
	  in[0].length + " and " + in[1].length);
	}
     
	if (loaded) {

	    //          If forward send data in normal order
	    if (dirn == FORWARD) {
		nativeTwoDimensionalSplit(width,height,in[0],in[1],out[0],
					  out[1],getPlanFlag());
	    }
	    else {   // For inverse swap real/imag parts
		nativeTwoDimensionalSplit(width,height,in[1],in[0],
					  out[1],out[0],getPlanFlag());
	    }
	}
	else {
	    errorStream.println("FFTWComplex.twoDimensional: library not loaded.");
	}
	
	return out;
    }
    
   
    /**   Method to take a Complex two-dimensional FFT with the real and imaging
     *    data split format in a two-dimensional array of size [2][length] 
     *    with the the i th component with real in the [0][i] element and the 
     *    imaginary element in the [1][i] element.
     *    @param width width of the data
     *    @param height height of the data
     *    @param in two-dimensional array, real in [0][i] and imag in [0][i]
     *    @param dirn forward or inverse transform
     *    @param overwrite in true it will overwrite the in array, else
     *      a new array will be created.
     *    @return <code>double[][]</code> transformed data
     */
    public double[][] twoDimensional(int width, int height, double[][] in,
				     int dirn, boolean overwrite) { 

	if (overwrite) {
	    return twoDimensional(width,height,in,in,dirn);
	}
	else {
	    return twoDimensional(width,height,in,
				  new double[2][in[0].length],dirn);
	}
    }

    /**     native TwoDimensionalSplit to do the work
     */
    private native void nativeTwoDimensionalSplit(int width, int height,
						  double realin[], double imagin[],
						  double realout[], double imagout[],
						  int flag);



    /**      Method to take three-dimensional Complex FFT with the data supplied
     *       in a one-dimensional double array with real parts in even elements,
     *       and imaginary in the odd.
     *
     *       Element (i,j,k) is located at:
     *
     *       RealPart  2*(k*width*height + j*width + i)
     *       ImagPart  2*(k*width*height * j*width + i) + 1
     *
     *       @param width the width of the image data
     *       @param height the height of the image data
     *       @param in input the data 
     *       @param out the output data, may be the same as in.
     *       @param dirn direction +1 for forward, -1 to backward
     *       @return double[] the FFTed array (if overwrite = true, this will
     *               be the same at the input data array.
     *       @throws <code>IllegalArgumentException</code> if width, height and depth does 
     *        not match the array length.
     */
    
    public double[] threeDimensional(int width, int height, int depth, 
				     double in[], double out[], int dirn) {

	//             basic sanity check
	if (in == null || in.length == 0 || out == null || out.length == 0) {
	    throw new NullPointerException(
  "FFTComplex.threeDimensional: called with null or zero length data array.");
	}

	if (in.length != 2*width*height*depth || out.length != in.length ) {
            throw new IllegalArgumentException(
     "FFTWComplex.threeDimensional: called with width: inconsistent array lengths.");
	}
	
	if (loaded) {
	    nativeThreeDimensional(width,height,depth,in,out,dirn,getPlanFlag());
	}
	else {
	   errorStream.println("FFTWComplex.threeDimensional: library not loaded.");
	} 
	return out;
    }




    /**      Method to take three-dimensional Complex FFT with the data supplied
     *       in a one-dimensional double array with real parts in even elements,
     *       and imaginary in the odd.
     *
     *       Element (i,j,k) is located at:
     *
     *       RealPart  2*(k*width*height + j*width + i)
     *       ImagPart  2*(k*width*height * j*width + i) + 1
     *
     *       @param width the width of the image data
     *       @param height the height of the image data
     *       @param in the data to be transformed
     *       @param dirn direction +1 for forward, -1 to backward
     *       @param overwrite if true data is overwritten with FFT, else a
     *              new array is returned.
     *       @return double[] the FFTed array (if overwrite = true, this will
     *               be the same at the input data array.
     *       @throws <code>IllegalArgumentException</code> if width, height and depth does 
     *       not match the array length.
     */
    
    public double[] threeDimensional(int width, int height, int depth, 
				     double in[], int dirn, boolean overwrite) {

	if (overwrite) {                                     // In place transform
	    return threeDimensional(width,height,depth,in,in,dirn);
	}
	else {      
	    return threeDimensional(width,height,depth,in,new double[in.length],dirn);
	}
    }
    
    /**       nativeTwoDimensional method to do the work
     */
    private native void nativeThreeDimensional(int width, int height, int depth, 
					       double in[], 
					       double out[], int dirn, int flag);



    /**   Method to take a Complex three-dimensional FFT with the real and imaging
     *    data split format in a two-dimensional array of size [2][length] 
     *    with the the i th component with real in the [0][i] element and the 
     *    imaginary element in the [1][i] element.
     *    @param width width of the data
     *    @param height height of the data
     *    @param depth depth of the data
     *    @param in two-dimensional array, real in [0][i] and 
     *            imag in [1][i]
     *    @param out the output array.
     *    @param dirn forward or inverse transform
     *    @return <code>double[][]</code> transformed data
     *    @throws <code>IllegalArgumentException</code> if real and imaginary 
     *     array of different length, or
     *     width*height does not match the array lengths.
     */
    public double[][] threeDimensional(int width, int height, int depth, 
				       double[][] in, 
				       double[][] out, int dirn) {
	//           Basic  sanity check
	if (in == null || in.length != 2 || in[0].length == 0 || 
	    in[1].length == 0) {
	    throw new NullPointerException(
"FFTComplex.threeDimensional: called with null or zero length data array(s).");
	}
	if (out == null || out.length != 2 || out[0].length == 0 || 
	    out[1].length == 0) {
	    throw new NullPointerException(
"FFTComplex.threeDimensional: called with null or zero length data array(s).");
	}
	if ( in[0].length != in[1].length || 
	     in[0].length != width*height*depth ||
	     out[0].length != out[1].length ||
	     in[0].length != out.length) {
	    throw new IllegalArgumentException(
      "FFTWComplex.threeDimensional: called with width : " + width +
      " height: " + height + " depth: " + depth +
      " but array length(s) of " +
      in[0].length + " and " + in[1].length);
	}
 
	if (loaded) {
	    //          If forward send data in normal order
	    if (dirn == FORWARD) {
		nativeThreeDimensionalSplit(width,height,depth,
					    in[0],in[1],out[0],out[1],
					    getPlanFlag());
	    }
	    else {   // For inverse swap real/imag parts
		nativeThreeDimensionalSplit(width,height,depth,
					    in[1],in[0],
					    out[1],out[0],getPlanFlag());
	    }
	}
	else {
	    errorStream.println("FFTWComplex.threeDimensional: library not loaded.");
	} 
	return out;
    }



    /**   Method to take a Complex three-dimensional FFT with the real and imaging
     *    data split format in a two-dimensional array of size [2][length] 
     *    with the the i th component with real in the [0][i] element and the 
     *    imaginary element in the [1][i] element.
     *    @param width width of the data
     *    @param height height of the data
     *    @param depth depth of the data
     *    @param in two-dimensional array, real in [0][i] and 
     *            imag in [1][i]
     *    @param dirn forward or inverse transform
     *    @param overwrite if true it the output will overwrite
     *       else a new array will be formed.
     *    @return <code>double[][]</code> transformed data
     *    @throws <code>IllegalArgumentException</code> if real and imaginary 
     *     array of different length, or
     *     width*height does not match the array lengths.
     */
    public double[][] threeDimensional(int width, int height, int depth, 
				       double[][] in, 
				       int dirn, boolean overwrite) {
	if (overwrite) {
	    return threeDimensional(width,height,depth,in,in,dirn);
	}
	else {
	    return threeDimensional(width,height,depth,in,
				    new double[2][in[0].length],dirn);
	}
    }


    
    /**     native TwoDimensionalSplit to do the work
     */
    private native void nativeThreeDimensionalSplit(int width, int height, int depth,
						  double realin[], double imagin[],
						  double realout[], double imagout[],
						  int flag);
    


}
    
    
