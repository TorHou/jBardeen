import java.io.*;
import java.util.ArrayList;


public class runIt{

	static String pathName;
	static boolean haveCoefficientsBeenRead = false;
	static Coefficients coefficients = null;
	
	public static void main(String[] args) throws Exception {
		
		if (args.length < 4){
			System.out.println("usage:");
			System.out.println("java -jar prog.jar path z-distance eFermi voltage1 [... voltageN]");
			System.exit(1);
		}
		
		// fermi CuBdD2: -1.5780, CuBdLDA: -1.7050
		
		// make sure arguments are formed well
		pathName = args[0];
		if (pathName.charAt(pathName.length()-1) != '/') pathName += '/';
		System.out.println("path: " + pathName);
		
		double zDistance = Double.parseDouble(args[1]);
		double eFermi = Double.parseDouble(args[2]);
		
		System.out.println("Distance: "+zDistance);
		System.out.println("Fermi Energy: "+eFermi+" eV");
	
		ArrayList<Double> biases = new ArrayList<Double>();
		for (int i=3;i<args.length;i++){
			biases.add(Double.parseDouble(args[i]));
		}
			
		WavecarIO_2 wcar1 = new WavecarIO_2(pathName, "WAVECAR", coefficients, zDistance, eFermi, biases);
		
		
	}
}

class Coefficients {
	public float[][][][] surfaceCoefficients;
	public float[][][][] tipCoefficients;
	public Coefficients(float[][][][] surfCoeff, float[][][][] tipCoeff){
		surfaceCoefficients = surfCoeff;
		tipCoefficients = tipCoeff;
	}
}