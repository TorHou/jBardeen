

/**      Demo programms to shown the use of the FFTWComplex class
 *       to take low level FFT one-dimensional tranforms 
 *       of interleaved data.
 *
 *       Note if you have interleaved data, or are willing to use
 *       interleaved data in your code, then consider using the
 *       ComplexDataArray class. They are much easier to call.
 *       @author Will Hossack, 2008
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
package jfftw.demo;
import jfftw.*;


public class FFTWComplexDemo {

    public static void main(String args[]) {

	int n = 2048;                   // Length of 1D-FFT

	//         Form in/out buffers for interleaved complex
	//         array.
	
	double[] inBuffer = new double[2*n];
	double[] outBuffer = new double[inBuffer.length];


	//         Fill the real part of the inBuffer with randoms 
	//         and imaginary with constant using
	//         the ArrayUtil method to access the ith complex
	//         element with two doubles
	for(int i = 0; i < n; i++) {
	    ArrayUtil.setComplex(inBuffer,i,Math.random() + 1.0, 1.0);
	}


	//           Form a Complex FFT object to do the work
	FFTWComplex fft = new FFTWComplex();

	//          Take out-of place forward transform.
	//          inBuffer will not be changed, and
	//          outBuffer will hold the FFT.
	outBuffer = fft.oneDimensional(inBuffer,outBuffer,FFTW.FORWARD);


	//          Take in-place backward transform of outBuffer
	//          This will overwrite outBuffer with its
	//          inverse FFT
	outBuffer = fft.oneDimensional(outBuffer,FFTW.BACKWARD,true);

	//         There is no normalisation, so
	//         after Forward FFT + Inverse FFT, each element of 
	//         outBuffer will be a  factor of n larger than expected.
	//         Use ArrayUtil method to normalise. 
	ArrayUtil.mult(outBuffer,1.0/n);


	//         inBuffer and outBuffer should now the
	//         (almost) identical. Form the real and imaginary
	//         modulus differences, this time using direct array
	//         access. 
	//         For the ith Complex element, we have
	//         Real part in array element 2*i
	//         Imaginary part in array element 2*i + 1
	double realDiff = 0.0;
	double imagDiff = 0.0;
	for(int i = 0; i < n; i++) {
	    int j = 2*i;
	    realDiff += Math.abs(inBuffer[j] - outBuffer[j]);
	    j++;
	    imagDiff += Math.abs(inBuffer[j] - outBuffer[j]);
	}

	//          Print out the differences, these will be small
	//          but NOT zero due to numerical rounding errors.
	System.out.println("Real Difference is : " + realDiff +
			   "\nImaginary Difference : " + imagDiff);


    }
}