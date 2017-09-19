package jfftw;
/**
 *         Class to handle Complex numbers and arithmetic.
 *         The complex is held a two double, and all arithmetic
 *         is in double format.
 *
 *         Much of the internal coding does not use methods for
 *         efficiency reasons.
 *
 *         @author Will Hossack, 2008
 *         @version 3.0
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
 */
public class Complex extends Object implements Cloneable {

    /**       Static to specify real part of Complex number
     */
    public static final int REAL = 0x1;

    /**       Static to specify imaginary part of Complex number
     */
    public static final int IMAG = 0x2;

    /**        Static to specify both parts of a Complex number
     */
    public static final int BOTH = 0x3;

    /**        Static to specify the modulus of a Complex number
     */
    public static final int MODULUS = 0x8;

    /**        Static to specify the phase of a Complex number
     */
    public static final int PHASE = 0x10;

       /**        Static to specify the modulus squared of a Complex number
     */
    public static final int MODULUS_SQUARED = 0x20;

    /**        Static to specify Log Power of a Complex
     */
    public static final int LOG_POWER = 0x40;


    /**           The real part
     */
    public double r;

    /**           The imaginary part
     */
    public double i;

    private static String fmt = "%g";

    /**   Constructor for Complex with two doubles
     *    @param a real part
     *    @param b imaginary part
     */
    public Complex(double a, double b) {
	set(a,b);
    }

    /**  Constructor for Complex with real part only
     *   @param a real part
     */
    public Complex(double a) {
	set(a,0.0);
    }

    /**  Default constructor for Complex with both real and imaginary zero
     */
    public Complex() {
	set(0.0,0.0);
    }

    /**  Constructor for Complex with Complex parameter
     *   @param c the Complex value
     */
    public Complex(Complex c) {
	set(c);
    }


    /**   Method to clone the current Complex
     *    @return <code>Complex</code> clone of the current Complex 
     */
    public Complex clone() {
	return new Complex(this);
    }

    /**   Method reset Complex value 
     *    @param a the real part
     *    @param b the imaginary value
     */
    public void set(double a, double b) {
	r = a;
	i = b;
    }

    /**   Method to reset the Complex value
     *    @param c the Complex value
     */
    public void set(Complex c) {
	r = c.r;
	i = c.i;
    }

    /**   Method to set the Complex value with polar parameters
     *    @param rho the radial parameter
     *    @param theta the angular parameter
     */
    public void setPolar(double rho, double theta) {
	set(rho*Math.cos(theta),rho*Math.sin(theta));
    }


    /**   Method to set Complex value to expi(theta)
     *   @param theta the angular parameter
     */
    public void setExpi(double theta) {
	setPolar(1.0,theta);
    }


    /**   Method to set a Complex to a specified modulus, but
     *    the phase set randomly between 0 to 2pi
     *    @param m the modulus
     */
    public void setRandomPhase(double m){
	double theta = 2.0*Math.PI*Math.random();
	setPolar(m,theta);
    }


    /**   Method to set the current Complex to be invalid,
     *    with both real and imaginary set to NaN
     */
    public void setInvalid() {
	set(Double.NaN,Double.NaN);
    }


    /**     Test of either components set to NaN
     *      @return <code>boolean</code> true if NaN]
     */
    public boolean isNaN() {
	if (Double.isNaN(r) || Double.isNaN(i))
	    return true;
	else
	    return false;
    }


    /**      Set the real part, imaginary part unchanged
     *       @param a the real value
     */
    public void setReal(double a) {
	r = a;
    }

    /**     Set the imaginary part, real part unchanged
     *      @param b the imaginary part
     */
    public void setImag(double b) {
	i = b;
    }

    /**      Get the real part
     *       @return <code>double</code> the real part
     */
    public double getReal() {
	return r;
    }

    /**     Get the imag part
     *      @return <code>double</code> the imaginary part
     */
    public double getImag() {
	return i;
    }

    /**     Get the modulus squared of the Complex
     *      @return <code>double</code> the modulus squared
     */ 
    public double modulusSq() {
	return r*r + i*i;
    }

    /**    Get the modulus of the Complex
     *     @return <code>double</code> the modulus
     */
    public double modulus() {
	return Math.sqrt(r*r + i*i);
    }

    /**    Get the phase from atan2(i,r) method. 
     *     @return <code>double</code> the phase in range -pi -> pi
     */
    public double phase() {
	return Math.atan2(i,r);
    }



    /**    Get the log power, being defines at log(r*r + i*i + 1.0)
     *     @return <code>double</code> the log power
     */
    public double logPower() {
	return Math.log(r*r + i+i + 1.0);
    }


    /**       Method to get double, converted from current Complex 
     *        as specified by the conversion flag.
     *        @param flag the conversion flag
     */
    public double getDouble(int flag) {
   
	/*                    Switch to return the right value
	 */
	switch (flag) {
	case REAL:
	    return r;

	case IMAG:
	    return i;
	
	case MODULUS:
	    return modulus();
	
	case MODULUS_SQUARED:
	    return modulusSq();

	case PHASE: 
	    return phase();
	
	case LOG_POWER:
	    return logPower();
	
	default :                   // Not value, so NaN
	    return Double.NaN;
	}
    }
    

    /**    Add a complex to the current returning  a new Complex.
     *     @param c the complex
     *     @return <code>Complex</code> the sum
     */
    public Complex plus(Complex c) {
	return new Complex(r + c.r,i + c.i);
    }

    /**    Add a complex to the current returning a new Complex 
     *     @param a the real part
     *     @param b the imaginary part
     *     @return <code>Complex</code> the sum
     */
    public Complex plus(double a, double b) {
	return new Complex(r + a, i + b);
    }

    /**    Add a real number to the current returning a new Complex.
     *     @param a the real part
     *     @return <code>Complex</code> the sum
     */ 
    public Complex plus(double a) {
	return new Complex(r + a, i);
    }

    /**    Subtract a complex from the current retuning a new Complex
     *     @param c the complex
     *     @return <code>Complex</code> the result
     */
    public Complex minus(Complex c) {
	return new Complex(r - c.r,i - c.i);
    }

    /**    Subtract a complex from the current returning a new Complex
     *     @param a the real part
     *     @param b the imaginary part
     *     @return <code>Complex</code> the result
     */
    public Complex minus(double a, double b) {
	return new Complex(r - a, i - b);
    }

    /**    Subtract a real number from the current retuning a new Complex
     *     @param a the real part
     *     @return <code>Complex</code> the result
     */ 
    public Complex minus(double a) {
	return new Complex(r - a, i);
    }


    /**    Subtract the current Complex FROM specified value returning
     *     a new Complex
     *     @param c the complex
     *     @return <code>Complex</code> the result
     */
    public Complex from(Complex c) {
	return new Complex(c.r - r,c.i - i);
    }

    /**    Subtract the current Complex FROM specified value retuning
     *     a new Complex
     *     @param a the real part
     *     @param b the imaginary part
     *     @return <code>Complex</code> the result
     */
    public Complex from(double a, double b) {
	return new Complex(a - r, b - i);
    }

    /**    Subtract the current Complex from a scalar returning a new
     *     Complex
     *     @param a the real part
     *     @return <code>Complex</code> the result
     */ 
    public Complex from(double a) {
	return new Complex(a - r, -i);
    }


    /**    Multiply the current by a scalar returning a new Complex
     *     @param a the scalar
     *     @return <code>Complex</code> the multiplication
     */
    public Complex mult(double a) {
	return new Complex(a*r , a*i);
    }

    /**    Multiply the current by a specified Complex, returning
     *     a new Complex.
     *     @param a the real part
     *     @param b the imaginary part
     *     @return <code>Complex</code> the multiplication
     */
    public Complex mult(double a, double b) {
	return new Complex(a*r - b*i , a*i + b*r);
    }

    /**    Multiply the current by a specified Complex, returning
     *     a new Complex.
     *     @param c the Complex
     *     @return <code>Complex</code> the multiplication
     */
    public Complex mult(Complex c) {
	return new Complex(r*c.r - i*c.i, i*c.r + r*c.i);
    }

    /**    Multiply the current by the conjugate of the
     *     specified Complex retuning a new Complex
     *     @param c the Complex
     *     @return <code>Complex</code> the multiplication
     */
    public Complex multConj(Complex c) {
	return mult(c.r, -c.i);
    }
    
    /**    Divide the current Complex by a scalar returning a new Complex
     *     @param a the scalar
     *     @return <code>Complex</code> the result
     */
    public Complex over(double a) {
	return new Complex(r/a , i/a);
    }

    /**    Divide the current by specified Complex returning new Complex
     *     @param a the real part
     *     @param b the imaginary part
     *     @return <code>Complex</code> the result
     */ 
    public Complex over(double a, double b) {
	double m = 1.0/(a*a + b*b);
	return mult(a*m,-b*m);
    }

    /**    Divide the current by specified Complex returning new Complex
     *     @param c the Complex
     *     @return <code>Complex</code> the result
     */ 
    public Complex over(Complex c) {
	return over(c.r, c.i);
    }

    /**    Divide the specified scalar BY the current Complex returning
     *     a new Complex
     *     @param a the scalar
     *     @return <code>Complex</code> The result
     */
    public Complex under(double a) {
	double s = a/modulusSq();
	return new Complex(s*r, -s*i);
    }

    /**    Divide the specified Complex  BY the current Complex
     *     @param c the Complex
     *     @return <code>Complex</code> The result
     */
    public Complex under(Complex c) {
	return c.multConj(this).over(modulusSq());
    }
	
    /**    The Complex conjugate of the current Complex
     *     @return <code>Complex</code> the conjugate
     */
    public Complex conj() {
	return new Complex(r, -i);
    }

    /**    Form a unit modulus version current Complex retaining the
     *     phase
     *     @return <code>Complex</code> Unit modulus Complex
     */
    public Complex unity() {
	double m = modulus();
	if (m != 0.0 )
	    return over(m);
	else
	    return new Complex(1.0);
    }


    /**    Form a scaled Complex from the current Complex with specified
     *     modulus, while retaining phase.
     *     @param mod the specified modulus
     *     @return <code>Complex</code> the scaled Complex.
     */
    public Complex scale(double mod) {
	double m = modulus();
	if (m != 0.0) 
	    return mult(mod/m);
	else
	    return new Complex(mod);
    }


    /**     Format the current Complex as a String
     *      @return <code>String</code> formatted Complex
     */
    public String toString() {
	return String.format("[" + this.fmt + " , " + this.fmt +"]",
			     this.r,this.i);
    }

    /**     Change the default format, by default set to "%g"
     *      @param fmt the new format String
     */
    public void setFormatString(String fmt) {
	this.fmt = fmt;
    }


    /**     Method to add to the current Complex
     *      @param c the Complex to be added
     */
    public void addTo(Complex c) {
	r += c.r;
	i += c.i;
    }

    /**     Void method to add to the current Complex
     *      @param r real part to be added
     *      @param i imaginary part of be added
     */
    public void addTo(double r, double i) {
	this.r += r;
	this.i += i;
    }

    /**    Method to mult the current Complex by specified Complex
     *     @param a real part
     *     @param b imaginary part
     */
    public void multBy(double a, double b) {
	double rt = a*r - b*i;
	i = a*i + b*r;
	r = rt;
    }

    /**    Method to mult the current Complex by specified Complex
     *     @param c the specified Complex
     */
    public void multBy(Complex c) {
	multBy(c.r,c.i);
    }

    
    /**    Method to mult the current Complex by specified real
     *     @param a the multiplier
     */
    public void multBy(double a) {
	r *= a;
	i *= a;
    }


    /**     Static add method with two parameters
     *      @param a first Complex
     *      @param b second Complex
     *      @return <code>Complex</code> the addition
     */
    public static Complex add(Complex a, Complex b) {
	return a.plus(b);
    }

    /**     Static mult method with two parameters
     *      @param a first Complex
     *      @param b second Complex
     *      @return <code>Complex</code> the multiplication
     */
    public static Complex mult(Complex a, Complex b) {
	return a.mult(b);
    }

}
