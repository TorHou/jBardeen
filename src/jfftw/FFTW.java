package jfftw;
import java.io.*;

/**      Basic class for Java access to the <a href="http://www.fftw.org/">fftw3</a>
 *       library to take  Fourier Transforms or one-dimensional, 
 *       two-dimensional or three-dimensional data held 
 *       in <code>double[]</code> arrays.
 *       <p>
 *       This class set package defaults, manages the wisdom information
 *       and loads the sharable library. Actual numerical use of this package
 *       is via the extending classes <code>FFTWReal</code> and
 *       <code>FFTWComplex</code> 
 *       <p>
 *       Note this class is not declared <code>abstract</code> since
 *       direct instances of this class can be used to manipulate
 *       the static library parameters.
 *       <p>
 *       If you are a novice used, dont start here, look at the
 *       <code>DataArray</code>, at its exteding classes. These allow you
 *       to take simple FFTs without understanding the complicated bits!
 *
 *       @author Will Hossack, 2008
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
 *
 */
public class FFTW extends Object {


    /**      Version flag
     */
    public final static String version = "1.2";

    /**       Flag for Real space 
     */
    public final static int REAL = 0x1;

    /**       Flag for Fourier space 
     */
    public final static int FOURIER = 0x2;

    /**       Flag for forward FFT
     */    
    public final static int FORWARD = -1;

    /**       Flag for backward (on inverse) FFT
     */
    public final static int BACKWARD = 1;

    /**      Flag for exhaustive plan search (very slow)
     */
    public final static int EXHAUSTIVE = (1 << 3);

    /**      Flag for patient plan search (medium slow)
     */
    public final static int PATIENT = (1 << 5);

    /**     Flag for quick and dirty plan search (fast, but not optimal)
     *      <p>
     *      Note: this is default for jttfw but NOT the normal FFTW default
     */
    public final static int ESTIMATE = (1 << 6);

    
    /**    PrintStream for warning and error messages. Defaults to 
     *     <code>System.err</code>.
     */
    protected static PrintStream errorStream = System.err;
    
    /**    Integer to set to the current plan in use by the package.
     *     Defaults to FFTW.ESTIMATE. 
     */
    protected int planFlag = FFTW.ESTIMATE;

    /**     Public static boolen to determine if an excepction in thrown
     *      if the native library is missing. Setting this to <code>false</code>
     *      will allow the Java code to run but not any of the native
     *      library calls. Use this feature ONLY if you really 
     *      understand what you are doing!
     */
    public static boolean throwLoadException = true;

    /**     Boolean to determine in the sharable library has been loaded
     *      successfully. (mainly used in testing).
     */
    public static boolean loaded = false;

    private static String systemWisdomFile = "/etc/fftw/wisdom";

    /**      Default constructor, do nothing but load the sharable
     *       library. The library is statically loaded so it
     *       will only load a single copy.
     */
    public FFTW() {
	loadLibrary();	
    }

    /**      Constructor to form a new FFTW, being a copy of the
     *       specified one. The PlanFlag is copied.
     *       @param fft the FFTW 
     *       @since 1.1
     */
    public FFTW(FFTW fft) {
	loaded = fft.loaded;                     // If library loaded
	setPlanFlag(fft.getPlanFlag());          // Copy the plan flag
    }

    /**      Constructor with flag to optionally load the
     *       system wisdom file. If load fails there is a non-fatal warning.
     *       @param systemWisdom if true default system wisdom file
     */
    public FFTW(boolean systemWisdom) {
	if (systemWisdom) {
	    if (!loadWisdom()) {
		errorStream.println("FFTW warning: Failed to load system wisdom");
	    }
	}
    }
	 

    /**     Constructor to load specified wisdom file specified by name.
     *      If load fails there is a non-fatal warning.
     *      @param wisdomFileName the name of the wisdom file
     */
    public FFTW(String wisdomFileName) {
	if (!loadWisdom(wisdomFileName)) {
	    errorStream.println("FFTW warning: Failed to load wisdom " + 
			       wisdomFileName);
	}
    }


    /**     Method to get the version string
     *      @return <code>String</code>
     */
    public String getVersion() {
	return version;
    }

    /**     Method to set the plan flag, default is ESTIMATE
     *      <p>
     *      Since most jfftw call involve a one-off use of plans
     *      setting this of other options should be done with extreme
     *      care and understing of the implications.
     *      @param flag the plan flag
     */
    public void setPlanFlag(int flag) {
	planFlag = flag;
    }

    /**    Method to add to the plan flag
     *     @param addition the addition to add
     */
    public void addPlanFlag(int addition) {
	planFlag += addition;
    }

    /**     Method to get the plan flag
     *      @return int the plan flag
     */
    public int getPlanFlag() {
	return planFlag;
    }

    /**    Method to reset the default errorStream for package. The default
     *     is System.err
     *     @param es the error PrintStream
     *     @since 1.2 
     */
    public void setErrorStream(PrintStream es) {
	errorStream = es;
    }

    /**    Method to reset the system wisdom file location. Note this
     *     does not load the file.
     *     @param fileName name of default wisdom file
     */
    public void setSystemWisdom(String fileName) {
	systemWisdomFile = new String(fileName);
    }

    /**    Boolean method to load system wisdom file
     *     @return boolean true for success, false for failure
     */
    public boolean loadWisdom() {
	return loadWisdom(systemWisdomFile);
    }

    /**   Boolean method to load wisdom file from specified file by name.
     *    Note this used the FFTW library direct "C" io  and not
     *    via Java io interface. It thus can only real a file from the
     *    current local machine.
     *    @param fileName the wisdom file name
     *    @return boolean true for success, false for failure
     */
    public boolean loadWisdom(String fileName) {
	return nativeLoadWisdomFromFile(fileName);
    }

    /**  Private native method load wisdom from file
     */
    private native boolean nativeLoadWisdomFromFile(String fileName);


    /**     Boolean method to load the wisdom froma local
     *      Java <cde>String</code>
     */

    public boolean loadWisdomFromString(String wisdom) {
	return nativeLoadWisdomFromString(wisdom);
    }


    /**  Private native method load wisdom from String
     */
    private native boolean nativeLoadWisdomFromString(String wisdom);

    /**   Void method to clear the wisdom information
     */
    public void clearWisdom() {
	nativeClearWisdom();
    }

    /**  Private native to actually clear the wisdom information
     */
    private native void nativeClearWisdom();


    /** Boolean method to export the current wisdom to a file.
     *  This uses the FFTW native library call so is 
     *  via "C" io rather than the Java io interface. It can thus ONLY
     *  be written to a local file.
     *
     *  @param fileName the output file name in a String.
     *  @return <code>boolean</code> true if successful, else false
     */
    public boolean exportWisdom(String fileName) {
	return nativeExportWisdomToFile(fileName);
    }

    /**  Private native method to export the wisdom to a file
     */
    private native boolean nativeExportWisdomToFile(String fileName);

    
    /**   Method to get the wisdom information as a String 
     *    @return <code>String</code> returns the FFTW wisdom
     *           entry in a formatted String.
     */
    public String getWisdom() {
	return nativeGetWisdom();
    }

    /**  Private native method to get the wisdom as a Java string
     */
    private native String nativeGetWisdom();


    /**  Static method to read a Wisdom file into a String from a 
     *   <code>BufferedReader</code>
     *   @param in the input reader.
     *   @return String the wisdom as a String 
     */
    public static String readWisdom(BufferedReader in) {
	
	StringBuilder buffer = new StringBuilder();
	String line;
	try {
	    while((line = in.readLine()) != null) {
		buffer.append(line);               // Append the line
		buffer.append('\n');               // Append newline
	    }
	    in.close();
	} catch (Exception e) {
	    return null;
	}
	return buffer.toString();
    }

    /**  Static method to read a Wisdom file into a String from a 
     *   <code>File</code>
     *   @param file the wisdom file
     *   @return String the wisdom as a String 
     */
    public static String readWisdom(File file) {
	try {
	    BufferedReader in = new BufferedReader(new FileReader(file));
	    return  readWisdom(in);
	} catch (Exception e) {
	    errorStream.println("FFTW.readWisdon: failed to open " +
				file.getName());
	    return null;
	}
    }
    

    /**  Static method to write wisdom String to a file
     *   @param wisdom the wisdom String
     *   @param file the output file
     */
    public static boolean writeWisdom(String wisdom, File file) {
	
	try {
	    BufferedWriter out = new BufferedWriter(new FileWriter(file));
	    out.write(wisdom,0,wisdom.length());
	    out.close();
	} catch (Exception e) {
	    errorStream.println("FFTW.writeWisdom: failed to write : " +
			       file.getName());
	    return false;
	}
	return true;
    }


    /**      Force the shared library to be loaded here.
     *       Under Linux this library must be <code>libjFFTW3.so 
     *       </code>and either
     *       located in /usr/java/jdk<version>/jre/lib/i386/
     *       or under a location point at by the environmental variable
     *       LD_LIBRARY_PATH.
     */
    static void loadLibrary(){
	if (!loaded) {                                  // not loaded 
	    try {
		System.loadLibrary("jFFTW3");
		loaded = true;                          // Success
	    } catch (UnsatisfiedLinkError e) {
		if (throwLoadException) {              // Hard fail ???
		    throw new UnsatisfiedLinkError(
		     "Fail to load jFFTW3 library : " +
		     e.getMessage() + 
		     "\nCheck LD_LIBRARY_PATH");
		}
		else {
		    errorStream.println("FFTW: failed to load jFFT3 library"); 
		    loaded = false;
		}
	    }
	}
    }
       
}

    
