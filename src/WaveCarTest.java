import java.io.FileNotFoundException;
import java.io.IOException;


public class WaveCarTest {
	public WaveCarTest(String path, String wName){
		System.out.println("Name: "+ wName);
		readHeadData(path+"/"+wName+"/");
		
	}
	public void readHeadData(String wavecarName) {

		try {
			// mapping = new MappingscarIO(mappingscarName, pathName);
			LERandomAccessFile wcF = new LERandomAccessFile(wavecarName, "r");
			// wcF = new RandomAccessFile(wavecarName, "r");
			// BufferedWriter out = new BufferedWriter(new FileWriter(pathName
			// + "out.dat"));

			long wavecarRecordLength = (long) wcF.readDouble();
			//wavecarRecordLength *= 4; // sigh ... wegen byterecl flag bei der VASP kompilierung
			// wavecarRecordLength = (long) wcF2.readDouble();

			System.out.println("Record-Length: " + wavecarRecordLength);
			int nrSpin = (int) wcF.readDouble();
			System.out.println("spins" + nrSpin);
			
			double version = wcF.readDouble();
			System.out.println("version: " + version);
			
			if ((int) version != 45200 && (int) version != 45210)
				throw new WavecarCorruptException();

			wcF.seek(wavecarRecordLength);

			double nrKPoints = wcF.readDouble();
			double[] kPoints = new double[3];
			double bands = (int) wcF.readDouble();
			double enMax = (int) wcF.readDouble();
			double[][] latticeVectorsAU = new double[3][3];
			double[][] latticeVectorsAngst = new double[3][3];
			double[][] reciprocalVectorsAngst = new double[3][3];
			double[][] reciprocalVectorsAU = new double[3][3];
			
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					latticeVectorsAngst[i][j] = wcF.readDouble();
					latticeVectorsAU[i][j] = latticeVectorsAngst[i][j] / 0.529;
				}
			}
			System.out.println("latticevector "+ latticeVectorsAngst[1][1]);

			double[] productAU01 = vectorialProduct(latticeVectorsAU[0],
					latticeVectorsAU[1]);
			double[] productAngst01 = vectorialProduct(latticeVectorsAngst[0],
					latticeVectorsAngst[1]);
			double[] productAU12 = vectorialProduct(latticeVectorsAU[1],
					latticeVectorsAU[2]);
			double[] productAngst12 = vectorialProduct(latticeVectorsAngst[1],
					latticeVectorsAngst[2]);
			double[] productAU20 = vectorialProduct(latticeVectorsAU[2],
					latticeVectorsAU[0]);
			double[] productAngst20 = vectorialProduct(latticeVectorsAngst[2],
					latticeVectorsAngst[0]);
			
			
			reciprocalVectorsAngst[0][0]=2*Math.PI*productAngst12[0]/ (dotProduct(latticeVectorsAngst[0], productAngst12));
			reciprocalVectorsAngst[0][1]=2*Math.PI*productAngst12[1]/ (dotProduct(latticeVectorsAngst[0], productAngst12));
			reciprocalVectorsAngst[0][2]=2*Math.PI*productAngst12[2]/ (dotProduct(latticeVectorsAngst[0], productAngst12));
			reciprocalVectorsAU[0][0]=2*Math.PI*productAU12[0]/ (dotProduct(latticeVectorsAU[0], productAU12));
			reciprocalVectorsAU[0][1]=2*Math.PI*productAU12[1]/ (dotProduct(latticeVectorsAU[0], productAU12));
			reciprocalVectorsAU[0][2]=2*Math.PI*productAU12[2]/ (dotProduct(latticeVectorsAU[0], productAU12));
			reciprocalVectorsAngst[1][0]=2*Math.PI*productAngst20[0]/ (dotProduct(latticeVectorsAngst[0], productAngst12));
			reciprocalVectorsAngst[1][1]=2*Math.PI*productAngst20[1]/ (dotProduct(latticeVectorsAngst[0], productAngst12));
			reciprocalVectorsAngst[1][2]=2*Math.PI*productAngst20[2]/ (dotProduct(latticeVectorsAngst[0], productAngst12));
			reciprocalVectorsAU[1][0]=2*Math.PI*productAU20[0]/ (dotProduct(latticeVectorsAU[0], productAU12));
			reciprocalVectorsAU[1][1]=2*Math.PI*productAU20[1]/ (dotProduct(latticeVectorsAU[0], productAU12));
			reciprocalVectorsAU[1][2]=2*Math.PI*productAU20[2]/ (dotProduct(latticeVectorsAU[0], productAU12));
			reciprocalVectorsAngst[2][0]=2*Math.PI*productAngst01[0]/ (dotProduct(latticeVectorsAngst[0], productAngst12));
			reciprocalVectorsAngst[2][1]=2*Math.PI*productAngst01[1]/ (dotProduct(latticeVectorsAngst[0], productAngst12));
			reciprocalVectorsAngst[2][2]=2*Math.PI*productAngst01[2]/ (dotProduct(latticeVectorsAngst[0], productAngst12));
			reciprocalVectorsAU[2][0]=2*Math.PI*productAU01[0]/ (dotProduct(latticeVectorsAU[0], productAU12));
			reciprocalVectorsAU[2][1]=2*Math.PI*productAU01[1]/ (dotProduct(latticeVectorsAU[0], productAU12));
			reciprocalVectorsAU[2][2]=2*Math.PI*productAU01[2]/ (dotProduct(latticeVectorsAU[0], productAU12));

			double volumeAngst=(dotProduct(latticeVectorsAngst[0], productAngst12));
			double volumeAU=(dotProduct(latticeVectorsAU[0], productAU12));
			double surfaceAU =(Math.sqrt(Math.pow(productAU01[0],2)+Math.pow(productAU01[1],2)+Math.pow(productAU01[2],2)));
			//System.out.println("volumeAngst "+ volumeAngst);
			
			double phi12 = Math.acos(dotProduct(reciprocalVectorsAngst[0], reciprocalVectorsAngst[1])/(vectorLength(reciprocalVectorsAngst[0])*vectorLength(reciprocalVectorsAngst[1])));
			System.out.println("phi12: "+ phi12/Math.PI*180 );
			
			double[] Xproductb01 = vectorialProduct(reciprocalVectorsAngst[0], reciprocalVectorsAngst[1]);
			double sinphi012 = (dotProduct(reciprocalVectorsAngst[2], Xproductb01))/(vectorLength(Xproductb01)*vectorLength(reciprocalVectorsAngst[2]));
			

					
			
//			int nb1Max = (int) (Math.sqrt(enMax * c) / (vectorLength(reciprocalVectorsAngst[0])*Math.abs(Math.sin(phi12)))) + 1;
//			nb2Max = (int) (Math.sqrt(enMax * c) / (vectorLength(reciprocalVectorsAngst[1])*Math.abs(Math.sin(phi12)))) + 1;
//			nb3Max = (int) (Math.sqrt(enMax * c) / (vectorLength(reciprocalVectorsAngst[2])*Math.abs(sinphi012))) + 1;
//			nb1Max = (int) (Math.sqrt(enMax * c) / (vectorLength(reciprocalVectorsAngst[0])*Math.abs(Math.sin(phi12)))) ;
//			nb2Max = (int) (Math.sqrt(enMax * c) / (vectorLength(reciprocalVectorsAngst[1])*Math.abs(Math.sin(phi12)))) ;
//			nb3Max = (int) (Math.sqrt(enMax * c) / (vectorLength(reciprocalVectorsAngst[2])*Math.abs(sinphi012))) ;
//			ng1 = getAllowedGridSize(nb1Max);
//			ng2 = getAllowedGridSize(nb2Max);
//			ng3 = getAllowedGridSize(nb3Max);
//			//ng1 = 112; ng2 = 112; ng3 = 120; //CuPc
//			
//			System.out.println(this.toStringSmall());
//			

			
			
//			decay = new Decay(ng1, ng2, reciprocalVectorsAU);
//			if (debug > 0)
//				printMatrixToFile(decay.matrix, decay.max, "heyo.png");

			if (latticeVectorsAU[1][2] > 0.05 || latticeVectorsAU[1][2] > 0.05) {
				System.err
						.println("3rd axis has to be perpendicular to the other axis. WAVECAR shows differently.");
				System.exit(1);
			}
		} catch (FileNotFoundException e) {
			System.out.println(e.toString());
		} catch (IOException e) {
			System.out.println(e.toString());
		} catch (WavecarCorruptException e) {
			System.out.println(e.toString());
		}

	}
	
	public double[] vectorialProduct(double[] v1, double[] v2) {
		double[] product = new double[3];
		product[0] = v1[1] * v2[2] - v1[2] * v2[1];
		product[1] = v1[2] * v2[0] - v1[0] * v2[2];
		product[2] = v1[0] * v2[1] - v1[1] * v2[0];
		return product;
	}
	
	public double dotProduct(double[] v1, double[] v2) {
		double product;
		product = v1[0] * v2[0]+v1[1] * v2[1]+v1[2] * v2[2];
		return product;
	}
	
	public double vectorLength(double[] v1) {
		return Math.sqrt(v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]);
	}

}
