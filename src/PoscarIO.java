import java.io.*;


public class PoscarIO {
	
	public double latticeParameter;
	public String poscarName;
	public double [][] latticeVectors;
	public int [] numberOfAtomsPerSpecies;
	public double volume;
	public double surface;
	
	public int numberOfSpecies = 0;
	private int numberOfAtoms;
	
	
	public PoscarIO(String path, String pname){
		poscarName = path + pname;
		
		latticeVectors = new double[3][3];
		
		PrintStream o = new PrintStream(System.out);
		
		
		try {
			FileInputStream inS = new FileInputStream(poscarName);
			InputStreamReader isR = new InputStreamReader(inS);
			BufferedReader file = new BufferedReader(isR);
			file.readLine();	// comment line
			String s1 = file.readLine().replaceAll("^\\s+","");
			Double factor = Double.valueOf(s1);	
			o.println(factor);
			for (int i=0; i<3; i++) {
				s1 = file.readLine().replaceAll("^\\s+","");
				String [] temp = s1.split("\\s+");
				for (int j=0; j<3; j++) latticeVectors[i][j] = Double.valueOf(temp[j])*factor;
				o.println(latticeVectors[i][0]+" "+latticeVectors[i][1]+" "+latticeVectors[i][2]+" ");
			}
			o.println(file.readLine());
			//System.out.println(file.readLine());
			
			s1 = file.readLine().replaceAll("^\\s+","");
			String [] temp = s1.split("\\s+");
			numberOfSpecies = temp.length;
			numberOfAtomsPerSpecies = new int[numberOfSpecies];
			numberOfAtoms=0;
			
			for (int i=0; i<temp.length; i++){
				numberOfAtomsPerSpecies[i] = Integer.valueOf(temp[i]);
				numberOfAtoms += numberOfAtomsPerSpecies[i];
				//System.out.println(numberOfAtomsPerSpecies[i]);
			}
			double [][] coordinatesDirect = new double[numberOfAtoms][3]; 
			double [][] coordinatesAngst = new double[numberOfAtoms][3];
			
			for(int i=0; i<5; i++) {
				s1 = file.readLine().replaceAll("^\\s+","");
				temp = s1.split("\\s+");
				if (isNumeric(temp[0])) break;
				//else	o.println("no number");
			}
			
			for (int j=0; j<3; j++) coordinatesDirect[0][j] = Double.valueOf(temp[j]);
			for (int i=1; i<numberOfAtoms; i++){
				s1 = file.readLine().replaceAll("^\\s+","");
				temp = s1.split("\\s+");
				for (int j=0; j<3; j++) coordinatesDirect[i][j] = Double.valueOf(temp[j]);
				
				coordinatesAngst[i] = convertDirectToAngst(coordinatesDirect[i]);
				o.println(coordinatesDirect[i][0]+" "+coordinatesDirect[i][1]+" "+coordinatesDirect[i][2]+" ");
				o.println(coordinatesAngst[i][0]+" "+coordinatesAngst[i][1]+" "+coordinatesAngst[i][2]+" ");
			}
			
			/*double[] vectorProduct1 = new double[3];
			vectorProduct1 = vectorialProduct(latticeVectors[0], latticeVectors[1]);
			//System.out.println(vectorProduct1[0]+"\t"+vectorProduct1[1]+"\t"+vectorProduct1[2]);
			
			volume = dotProduct(vectorProduct1,latticeVectors[2]);
			
			for (int i=0; i<3; i++){
				surface += vectorProduct1[i]*vectorProduct1[i];
			}
			surface = Math.sqrt(surface);
			*/
			//System.out.println(s1);
			//System.out.println(s2);
			//System.out.println("1st field: "+temp[0]);
			//System.out.println("2nd field: "+temp[1]);
			
			
			//parObject.wavecarRecordLength = (int) f1;
			
			
				
			//parObject.latticeVectors=parObject.latticeVectors/0.529;
			if (inS != null) {
				inS.close();
			}
			
		} catch (FileNotFoundException e){
			System.out.println(e.toString());
		} catch (IOException e){
			System.out.println(e.toString());
		} 
		
	}
	
	private double[] convertDirectToAngst(double[] ds) {
		double [] result = new double[3];
		for (int i=0; i<3; i++) 
			result[i] = dotProduct(ds, latticeVectors[i]);
		return result;
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
		String output = "Name: "+poscarName+"\n"; 
		output += "lattice parameter: "+latticeParameter+"\n";
		output += "species: "+numberOfSpecies+"\n";
		for (int i=0; i<numberOfSpecies; i++) output += numberOfAtomsPerSpecies[i]+"\t";
		output += "\n";
		output += "Lattice Vector 1: "+latticeVectors[0][0]+"\t"+latticeVectors[0][1]+"\t"+latticeVectors[0][2]+"\n";
		output += "Lattice Vector 2: "+latticeVectors[1][0]+"\t"+latticeVectors[1][1]+"\t"+latticeVectors[1][2]+"\n";
		output += "Lattice Vector 3: "+latticeVectors[2][0]+"\t"+latticeVectors[2][1]+"\t"+latticeVectors[2][2]+"\n";
		output += "volume: "+volume+"\tsurface: "+surface+"\n";
		
		return output;
	}

	public static boolean isNumeric(String str) {  
		try  
		{  
			double d = Double.parseDouble(str);  
		}  
		catch(NumberFormatException nfe)  
		{  
			return false;  
		}  
		return true;  
	}
	
}	



class PoscarCorruptException extends Exception{
	public PoscarCorruptException(){
		super();
	}
}