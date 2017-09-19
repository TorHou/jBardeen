import java.awt.image.BufferedImage;
import java.awt.image.WritableRaster;
import java.io.*;

import javax.imageio.ImageIO;


public class MappingscarIO {
	public String mappingscarName;
	public double eFermiAU;
	public int kPoints;
	public int spin;
	public int ngx;
	public int ngy;
	public int ngz;
	public int surface;
	public int numberOfSpecies = 0;
	public LERandomAccessFile mappingscarFile;
	double [] gX;
	double [] gY;
	int[] indexX;
	int[] indexY;
	int[] indexZ;
	private int planeWaves;
	int debug=1;
	String pathName;
	
	public MappingscarIO(String name, String pathName){
		mappingscarName=name;
		this.pathName=pathName;
		//latticeVectors = new double[3][3];
		try {
			mappingscarFile = new LERandomAccessFile(name, "r");
			
			//inS = new FileInputStream(name);
			//isR = new InputStreamReader(inS);
			//dinS = new LittleEndianDataInputStream(inS);
			mappingscarFile.skipBytes(4);
			
			double f1 = mappingscarFile.readDouble();
			eFermiAU = f1/27.2116;
			mappingscarFile.skipBytes(8);
			spin = mappingscarFile.readInt();
			kPoints = mappingscarFile.readInt();
			ngx = mappingscarFile.readInt();
			ngy = mappingscarFile.readInt();
			ngz = mappingscarFile.readInt();
			mappingscarFile.skipBytes(8);
			
			
			
		} catch (FileNotFoundException e){
			System.out.println(e.toString());
		} catch (IOException e){
			System.out.println(e.toString());
		} 
		
	}
	
	public void getIndices (int planeWaves){
		try {
			this.planeWaves = planeWaves;
			gX = new double [ngx*ngy];
			gY = new double [ngx*ngy];
			indexX = new int [planeWaves];
			indexY = new int [planeWaves];
			indexZ = new int [planeWaves];
			
			int index=0;
			for (int i=0; i<ngx; i++){
				for (int j=0; j<ngy; j++){
					gX[index] = mappingscarFile.readDouble()*0.529;
					gY[index] = mappingscarFile.readDouble()*0.529;
					
					if (debug > 0) {
						System.out.println("gX: "+gX[index]+"\tgY: "+gY[index]);
					}
					mappingscarFile.skipBytes(8);
					index++;
				}
			}
			
			
			for (int i=0; i<planeWaves; i++){
				indexX[i] = mappingscarFile.readInt();
				indexY[i] = mappingscarFile.readInt();
				indexZ[i] = mappingscarFile.readInt();
				mappingscarFile.skipBytes(8);
				if (debug > 0) {
					if (i<100) System.out.println(i+": "+indexX[i]+" "+indexY[i]+" "+indexZ[i]);
				}
				//System.out.println(i+" "+indexX[i]+" "+indexY[i]+" "+indexZ[i]);
	
			}
			
		} catch (FileNotFoundException e){
			System.out.println(e.toString());
		} catch (IOException e){
			System.out.println(e.toString());
		}
		for (int ind=1; ind<10; ind++) printIndexMap(ind);
	}
	
	public void close(){
		try {
			mappingscarFile.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public Complex[] getColumn(int x, int y){
		Complex[] output = new Complex[0];
		return 	output;
	}
	
	
	public void printIndexMap(int z){
		BufferedImage indexImage = new BufferedImage(ngx, ngy, BufferedImage.TYPE_BYTE_GRAY);
		WritableRaster raster = indexImage.getRaster();
		
		for (int i=0; i<planeWaves; i++){
			
			if (indexZ[i]-1==z) {
				//System.out.println("heyo: "+i);
				//raster.setSample((indexX[i]-1+ngx/2)%(ngx-1), (indexY[i]-1+ngy/2)%(ngy-1), 0, 200);
				raster.setSample(indexX[i]-1, indexY[i]-1, 0, 200);
			}
		}
		try{
			File outputfile = new File(pathName+"indices"+z+".png");
		    ImageIO.write(indexImage, "png", outputfile);
		} catch (IOException e){
			e.printStackTrace();
		}
			
	}
	
	public void printGX(String fileName){
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(fileName));
			for (int i=0; i<ngx; i++){
				for (int j=0; j<ngy; j++){
					out.write(gX[ngy*i+j]/0.529+"\n");
					//out.write(gX[ngx*j+j]+"\n");
				}
			}
			if (out != null) {
				out.close();
			}
		} catch (FileNotFoundException e){
			System.out.println(e.toString());
		} catch (IOException e){
			System.out.println(e.toString());
		}
	}
	
	public void printGY(String fileName){
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(fileName));
			for (int i=0; i<ngx; i++){
				for (int j=0; j<ngy; j++){
					out.write(gY[ngy*i+j]/0.529+"\n");
					//out.write(gY[ngx*j+j]+"\n");
				}
			}
			if (out != null) {
				out.close();
			}
		} catch (FileNotFoundException e){
			System.out.println(e.toString());
		} catch (IOException e){
			System.out.println(e.toString());
		}
	}

	public void printIndices(String fileName, int planeWaves){
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(fileName));
			for (int i=0; i<planeWaves; i++){
				out.write("X: "+indexX[i]+"\tY: "+indexY[i]+"\tZ: "+indexZ[i]+"\n");
			}
			if (out != null) {
				out.close();
			}
		} catch (FileNotFoundException e){
			e.printStackTrace();
		} catch (IOException e){
			e.printStackTrace();
		}
	}

	
	public double[] vectorialProduct(double[] v1, double[] v2){
		double[] product=new double[3];
		product[0]= v1[1]*v2[2] - v1[2]*v2[1];
		product[1]= v1[2]*v2[0] - v1[0]*v2[2];
		product[2]= v1[0]*v2[1] - v1[1]*v2[0];			
		return product;
	}
	
	public double dotProduct(double[] v1, double[] v2){
		double product=0.0;
		product = v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];			
		return product;
	}
	
	public String toString(){
		String output = "Name: "+mappingscarName+"\n"; 
		output += "Fermi energy: "+eFermiAU+"\n";
		output += "Spin: "+spin+"\tKPoints: "+kPoints+"\n";
		output += "ngx: "+ngx+"\tngy: "+ngy+"\tngz: "+ngz+"\n";
		
		//for (int i=0; i<numberOfSpecies; i++) output += numberOfAtomsPerSpecies[i]+"\t";
		output += "\n";
		//output += "Lattice Vector 1: "+latticeVectors[0][0]+"\t"+latticeVectors[0][1]+"\t"+latticeVectors[0][2]+"\n";
		//output += "Lattice Vector 2: "+latticeVectors[1][0]+"\t"+latticeVectors[1][1]+"\t"+latticeVectors[1][2]+"\n";
		//output += "Lattice Vector 3: "+latticeVectors[2][0]+"\t"+latticeVectors[2][1]+"\t"+latticeVectors[2][2]+"\n";
		//output += "volume: "+volume+"\tsurface: "+surface+"\n";
		
		output += "\n";
		//output += "Lattice Vector 1: "+latticeVectors[0][0]+"\t"+latticeVectors[0][1]+"\t"+latticeVectors[0][2]+"\n";
		
		return output;
	}
	
}	

@SuppressWarnings("serial")
class MappingscarCorruptException extends Exception{
	public MappingscarCorruptException(){
		super();
	}
}