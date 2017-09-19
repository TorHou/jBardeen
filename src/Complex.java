public class Complex extends Object {

    private double x,y;
    
    /**
        Constructs the complex number z = u + i*v
        @param u Real part
        @param v Imaginary part
    */
    public Complex(double u,double v) {
        x=u;
        y=v;
    }
    
    public Complex() {
    	x=0;
    	y=0;
    }
    
    /**
        Real part of this Complex number 
        (the x-coordinate in rectangular coordinates).
        @return Re[z] where z is this Complex number.
    */
    public double real() {
        return x;
    }
    
    /**
        Imaginary part of this Complex number 
        (the y-coordinate in rectangular coordinates).
        @return Im[z] where z is this Complex number.
    */
    public double imag() {
        return y;
    }
    
    public String toString() {
    	return x + " + " + y + "i";
    }
}