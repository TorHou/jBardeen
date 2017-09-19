import java.awt.image.BufferedImage;
import java.awt.image.WritableRaster;
import java.io.*;
import java.text.DecimalFormat;
import java.util.ArrayList;

import javax.imageio.ImageIO;

import jfftw.*;

class WavecarIO_2 {
	int debug = 1;

	LERandomAccessFile wcF = null;
	RandomAccessFile wcF2 = null;

	int[] allowedGridSizes = { 2, 3, 5, 11, 12, 15, 16, 18, 20, 22, 24, 25, 27,
			30, 32, 33, 36, 40, 42, 44, 45, 48, 50, 54, 56, 60, 64, 70, 72, 80,
			84, 90, 96, 98, 100, 108, 112, 120, 126, 140, 144, 150 };
	int[][] indices;
	int[] indexX;
	int[] indexY;
	int[] indexZ;

	double c = 0.26246582250210965422;

	double eFermiAU;
	double zHighestDirect = 0.51;
	double zTipDirect = 0.99;
	int zSamplingPoints = 50;
	double phi = 4.5 / 27.2116;

	double zMaxAngst = 12;
	double zMaxAU = zMaxAngst / 0.529;
	double zMinAngst;

	double voltage;
	double voltageAU;
	double zDistance;

	int zStepsNr = 50;

	double sigma = 0.03 / 27.2116;

	double version;

	public double nrKPoints;
	public int bands;
	public int enMax;
	public int nrSpin;
	public long wavecarRecordLength;
	public String wavecarName;
	public String mappingscarName;
	public double[][] latticeVectorsAU;
	public double[][] latticeVectorsAngst;
	public double[][] reciprocalVectorsAU;
	public double[][] reciprocalVectorsAngst;
	public double surfaceAU;
	public double volumeAU;
	public double volumeAngst;
	public double[] kPoints;
	public int nrPlaneWaves;
	public double[][] eigenEnergiesAU;
	public double[][] eigenEnergiesEV;
	public float[] cReal;
	public float[] cImag;

	public float[][][][] coefficientMatrix;
	public float[][][][] coefficientMatrixReversed;

	public float[][][][] tipCoefficientsPerBand;
	public float[][][][] surfaceCoefficientsPerBand;
	public float[][][] surfaceCoefficients;
	public float[][][] tipCoefficients;

	public double[][][] matrixSurface;
	public double[][][] matrixTip;
	public double[][] tipTersoffHamann;
	public double[][] surfaceTersoffHamann;
	public double[][] tersoffHamannMatrixSurface;
	public double[][] tersoffHamannMatrixTip;
	public double[][] bardeenMatrix;

	public double[][][] tersoffHamannCurrent;
	public double[][][] bardeenCurrent;
	public double[][][] tersoffHamannSurface;
	public double[][][] tersoffHamannTip;
	public double[] zValuesAngst;
	ZIndices zIndices;

	DecimalFormat df2 = new DecimalFormat("#0.00");
	DecimalFormat df3 = new DecimalFormat("#0.000");
	DecimalFormat df4 = new DecimalFormat("#0.0000");

	public double[] ferTotal;
	public MappingscarIO mapping;
	DecimalFormat decimalForm = new DecimalFormat("#.####");
	private int nb1Max;
	private int nb2Max;
	private int nb3Max;
	private int ng1;
	private int ng2;
	private int ng3;
	private int recordCounter = 2;
	private String pathName;

	Coefficients coefficients;
	Decay decay = null;

	private int counter =0;
	ArrayList<Double> biasesAU = null;

	private String currentCoefficients ="";
	
	
	public WavecarIO_2(String path, String wName, 
			Coefficients coeff, double zDist, double eFermi, ArrayList<Double> biases) {
		pathName = path;
		biasesAU = new ArrayList<Double>();
		for (Double bias:biases){
			biasesAU.add(bias/27.2116);
			System.out.println("bias: "+bias);
		}
		
		//voltagesAU = volt / 27.2116;
		eFermiAU = eFermi / 27.2116;
		zDistance = zDist;
		wavecarName = path + wName;
		coefficients = coeff;
	
		String mname = "";
		mappingscarName = mname;
		readHeadData(wavecarName, mname);
		zMinAngst = (zTipDirect - zHighestDirect)
				* this.latticeVectorsAngst[2][2];

		// calculate zIndices for the matching plane
		zIndices = new ZIndices();

		
		for (Double biasAU:biasesAU){
			recordCounter = 2;
			System.out.println("biasAU "+biasAU);
			calculateCurrents(biasAU.doubleValue());
		}

	
	}

	public void printToFile(double[] z, double[][][] currents, String fileName) {
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(pathName
					+ fileName));
			for (int ix = 0; ix < ng1; ix++) {
				for (int iy = 0; iy < ng2; iy++) {
					for (int iz = 0; iz < zSamplingPoints; iz++) {
						out.write(zValuesAngst[iz] + "\t"
								+ currents[ix][iy][iz] + "\n");
					}
				}
			}
			if (out != null) {
				out.close();
			}
		} catch (FileNotFoundException e) {
			System.out.println(e.toString());
		} catch (IOException e) {
			System.out.println(e.toString());
		}

	}

	public int calculatePlanewaves() {
		int countPlaneWaves = 0;
		BufferedImage indexImage = new BufferedImage(nb1Max * 2 + 1,
				nb2Max * 2 + 1, BufferedImage.TYPE_BYTE_GRAY);
		// WritableRaster raster = indexImage.getRaster();
		indices = new int[3][2 * 2 * 2 * nb1Max * nb2Max * nb3Max];

		for (int ig3 = 0; ig3 <= 2 * nb3Max; ig3++) {
			// for (int ig3=0;ig3<2;ig3++){
			double ig3p = ig3;
			int ig3s = ig3;
			if (ig3 > nb3Max) {
				ig3p = ig3 - 2 * nb3Max - 1;
				ig3s = ig3 - 2 * nb3Max + ng3 - 1;
			}
			// System.out.println("ig3: "+ig3);
			for (int ig2 = 0; ig2 <= 2 * nb2Max; ig2++) {
				double ig2p = ig2;
				int ig2s = ig2;
				if (ig2 > nb2Max) {
					ig2p = ig2 - 2 * nb2Max - 1;
					ig2s = ig2 - 2 * nb2Max + ng2 - 1;
				}
				for (int ig1 = 0; ig1 <= 2 * nb1Max; ig1++) {
					double ig1p = ig1;
					int ig1s = ig1;
					if (ig1 > nb1Max) {
						ig1p = ig1 - 2 * nb1Max - 1;
						ig1s = ig1 - 2 * nb1Max + ng1 - 1;
					}
					double[] sum = new double[3];
					for (int i = 0; i < 3; i++) {
						sum[i] = (ig1p + kPoints[0])
								* reciprocalVectorsAngst[0][i];
						sum[i] += (ig2p + kPoints[1])
								* reciprocalVectorsAngst[1][i];
						sum[i] += (ig3p + kPoints[2])
								* reciprocalVectorsAngst[2][i];
					}
					double eTotal = sum[0] * sum[0] / c + sum[1] * sum[1] / c
							+ sum[2] * sum[2] / c;
					if (eTotal < enMax) {
						indices[0][countPlaneWaves] = ig1s;
						indices[1][countPlaneWaves] = ig2s;
						indices[2][countPlaneWaves] = ig3s;
						countPlaneWaves++;
						/*
						 * if (ig3==0) { raster.setSample(ig1, ig2, 0, 200); if
						 * (ig2==0)
						 * System.out.println("ig1: "+ig1+" ig2: "+ig2+" ig3: "
						 * +ig3); if (ig2==0)
						 * System.out.println("ig1s: "+ig1s+" ig2s: "
						 * +ig2s+" ig3s: "+ig3s); }
						 */
					}

				}
			}
		}
		// System.out.println("Number of plane waves counted: "+nrPlaneWaves);
		try {
			File outputfile = new File(pathName + "indices" + 0 + ".png");
			ImageIO.write(indexImage, "png", outputfile);
		} catch (IOException e) {
			e.printStackTrace();
		}

		return countPlaneWaves;
	}

	public int calculatePlanewavesZHalf() {
		int countPlaneWaves = 0;

		BufferedImage indexImage = new BufferedImage(nb1Max * 2 + 1,
				nb2Max * 2 + 1, BufferedImage.TYPE_BYTE_GRAY);
		// WritableRaster raster = indexImage.getRaster();
		indices = new int[3][2 * 2 * nb1Max * nb2Max * nb3Max];

		for (int ig3 = 0; ig3 <= nb3Max; ig3++) {
			double ig3p = ig3;
			int ig3s = ig3;
			for (int ig2 = 0; ig2 <= 2 * (nb2Max); ig2++) {
				double ig2p = ig2;
				int ig2s = ig2;
				if (ig2 > nb2Max) {
					ig2p = ig2 - 2 * nb2Max - 1;
					ig2s = ig2 - 2 * nb2Max + ng2 - 1;
				}
				for (int ig1 = 0; ig1 <= 2 * (nb1Max); ig1++) {
					double ig1p = ig1;
					int ig1s = ig1;
					if (ig1 > nb1Max) {
						ig1p = ig1 - 2 * nb1Max - 1;
						ig1s = ig1 - 2 * nb1Max + ng1 - 1;
					}
					if (ig3p == 0 && ig2p < 0)
						break;
					if (ig3p == 0 && ig2p == 0 && ig1p < 0)
						break;

					double[] sum = new double[3];
					for (int i = 0; i < 3; i++) {
						sum[i] = (ig1p + kPoints[0])
								* reciprocalVectorsAngst[0][i];
						sum[i] += (ig2p + kPoints[1])
								* reciprocalVectorsAngst[1][i];
						sum[i] += (ig3p + kPoints[2])
								* reciprocalVectorsAngst[2][i];
					}
					double eTotal = sum[0] * sum[0] / c + sum[1] * sum[1] / c
							+ sum[2] * sum[2] / c;

					if (eTotal < enMax) {

						indices[0][countPlaneWaves] = ig1s;
						indices[1][countPlaneWaves] = ig2s;
						indices[2][countPlaneWaves] = ig3s;
						countPlaneWaves++;
						/*
						 * if (ig3==0) { raster.setSample(ig1, ig2, 0, 200); if
						 * (ig2==0)
						 * System.out.println("ig1: "+ig1+" ig2: "+ig2+" ig3: "
						 * +ig3); if (ig2==0)
						 * System.out.println("ig1s: "+ig1s+" ig2s: "
						 * +ig2s+" ig3s: "+ig3s); }
						 */
					}

				}
			}
		}
		// System.out.println("Number of plane waves counted: "+nrPlaneWaves);
		try {
			File outputfile = new File(pathName + "indices" + 0 + ".png");
			ImageIO.write(indexImage, "png", outputfile);
		} catch (IOException e) {
			e.printStackTrace();
		}

		return countPlaneWaves;
	}

	public void calculateTHConductance(int recordCounter) {
		for (int kp = 0; kp < (int) nrKPoints; kp++) {
			for (int sp = 0; sp < 1; sp++) {

				readEnergies(recordCounter);
				recordCounter++;

				int cPlaneWaves = calculatePlanewaves();
				checkForCorrectNumberOfPlanewaves(cPlaneWaves);

				initializeArraysNew();

				System.out.println("k-Point: " + (kp + 1)
						+ ", Spin-component: " + (sp + 1));
				String coefficientsFileName = pathName + "coefficients-k"
						+ (kp + 1) + "spin" + (sp + 1) + "surfaceZ"
						+ zIndices.surface + "tipZ" + zIndices.tip;
				if (counter==1){}
				else if ((new File(coefficientsFileName)).isFile()){
					readCoefficientsFromFile(coefficientsFileName);
					counter=1;
				}
					
				else
					calculateCoefficients(coefficientsFileName);

				System.out.print("calculating currents ");

				int progressPercentageThreshold = 10;

				double zStepAU = zMaxAU / (zStepsNr - 1);

				for (int iz = 0; iz < zSamplingPoints; iz++) {
					double zAU = zStepAU * iz;
					zValuesAngst[iz] = zAU * 0.529 + zMinAngst;

					/* ******************************************************
					 * Bardeen
					 * ******************************************************
					 */
					for (int tb = 0; tb < bands; tb++) {
						for (int sb = 0; sb < bands; sb++) {
							double spinFactor = 2;
							double constant_I = 2 * Math.sqrt(Math.PI) / sigma
									* Math.pow(surfaceAU, 2)
									/ Math.pow(volumeAU, 2) * spinFactor;
							if (voltage < 0)
								constant_I = -constant_I;

							if ((voltage < 0
									&& eigenEnergiesAU[tb][0] - eFermiAU > 0 && eigenEnergiesAU[sb][0]
									- eFermiAU < 0)
									|| (voltage > 0
											&& eigenEnergiesAU[tb][0]
													- eFermiAU < 0 && eigenEnergiesAU[sb][0]
											- eFermiAU > 0)) {
								double deltaE = Math
										.pow(((eigenEnergiesAU[tb][0] - eigenEnergiesAU[sb][0]) + voltageAU)
												/ sigma, 2);
								if (deltaE < 3.0) {
									// System.out.println("deltaE<3");
									for (int ix = 0; ix < ng1; ix++) {
										for (int iy = 0; iy < ng2; iy++) {
											surfaceCoefficients[ix][iy][0] = surfaceCoefficientsPerBand[ix][iy][0][sb];
											surfaceCoefficients[ix][iy][1] = surfaceCoefficientsPerBand[ix][iy][1][sb];
											tipCoefficients[ix][iy][0] = tipCoefficientsPerBand[ix][iy][0][tb];
											tipCoefficients[ix][iy][1] = tipCoefficientsPerBand[ix][iy][1][tb];
										}
									}
									for (int ix = 0; ix < ng1; ix++) {
										for (int iy = 0; iy < ng2; iy++) {
											matrixSurface[ix][iy][0] = tipCoefficients[ix][iy][0]
													* surfaceCoefficients[ix][iy][0]
													* decay.matrix[ix][iy]
													* Math.exp(-decay.matrix[ix][iy]
															* zAU);
											matrixSurface[ix][iy][1] = -tipCoefficients[ix][iy][1]
													* surfaceCoefficients[ix][iy][1]
													* decay.matrix[ix][iy]
													* Math.exp(-decay.matrix[ix][iy]
															* zAU);
											// if (matrixSurface[ix][iy][0]>0)
											// System.out.println("ix: "+ix+" iy: "+iy+" matrix: "+matrixSurface[ix][iy][0]);

										}
									}

									double[] temp = new double[2 * ng2 * ng1];
									double[] ftout = new double[2 * ng2 * ng1];

									for (int ix = 0; ix < ng1; ix++) {
										for (int iy = 0; iy < ng2; iy++) {
											temp[ng2 * ix * 2 + iy * 2] = matrixSurface[ix][iy][0];
											temp[ng2 * ix * 2 + iy * 2 + 1] = matrixSurface[ix][iy][1];
										}
									}
									FFTWComplex fft = new FFTWComplex();
									fft.twoDimensional(ng2, ng1, temp, ftout,
											FFTW.BACKWARD);

									for (int ix = 0; ix < ng1; ix++) {
										for (int iy = 0; iy < ng2; iy++) {
											bardeenMatrix[ix][iy] = Math
													.pow(ftout[ng2 * ix * 2
															+ iy * 2], 2)
													+ Math.pow(ftout[ng2 * ix
															* 2 + iy * 2 + 1],
															2);

											bardeenCurrent[ix][iy][iz] += bardeenMatrix[ix][iy]
													* constant_I
													* Math.exp(-deltaE)
													/ (2 * nrKPoints - 1);
											// if (bardeenCurrent[ix][iy][iz]>0)
											// System.out.println("cu: ix: "+ix+" iy: "+iy+" "+bardeenMatrix[ix][iy]);

										}
									}

								}// if deltaE
							}// if voltage

						}// bands
					}// bands

					/* ******************************************************
					 * Tersoff Hamann
					 * ******************************************************
					 */
					for (int b = 0; b < bands; b++) {
						if (Math.abs(eigenEnergiesAU[b][0]) < Math
								.abs(voltageAU + eFermiAU)
								&& (voltageAU < 0
										&& eigenEnergiesAU[b][0] - eFermiAU < 0 || voltageAU > 0
										&& eigenEnergiesAU[b][0] - eFermiAU > 0)) {
							// if (iz==0) System.out.println("band: "+b);
							for (int ix = 0; ix < ng1; ix++) {
								for (int iy = 0; iy < ng2; iy++) {
									surfaceCoefficients[ix][iy][0] = surfaceCoefficientsPerBand[ix][iy][0][b];
									surfaceCoefficients[ix][iy][1] = surfaceCoefficientsPerBand[ix][iy][1][b];
									tipCoefficients[ix][iy][0] = tipCoefficientsPerBand[ix][iy][0][b];
									tipCoefficients[ix][iy][1] = tipCoefficientsPerBand[ix][iy][1][b];
								}
							}

							for (int ix = 0; ix < ng1; ix++) {
								for (int iy = 0; iy < ng2; iy++) {
									// System.out.println(b+" "+ix+" "+iy);
									matrixSurface[ix][iy][0] = surfaceCoefficients[ix][iy][0]
											* Math.exp(-decay.matrix[ix][iy]
													* zAU);
									matrixSurface[ix][iy][1] = surfaceCoefficients[ix][iy][1]
											* Math.exp(-decay.matrix[ix][iy]
													* zAU);
									matrixTip[ix][iy][0] = tipCoefficients[ix][iy][0]
											* Math.exp(-decay.matrix[ix][iy]
													* zAU);
									matrixTip[ix][iy][1] = tipCoefficients[ix][iy][1]
											* Math.exp(-decay.matrix[ix][iy]
													* zAU);

								}
							}

							double[] temp = new double[2 * ng2 * ng1];
							double[] temp2 = new double[2 * ng2 * ng1];
							double[] ftout = new double[2 * ng2 * ng1];
							double[] ftout2 = new double[2 * ng2 * ng1];

							for (int ix = 0; ix < ng1; ix++) {
								for (int iy = 0; iy < ng2; iy++) {
									temp[ng2 * ix * 2 + iy * 2] = matrixSurface[ix][iy][0];
									temp[ng2 * ix * 2 + iy * 2 + 1] = matrixSurface[ix][iy][1];
									temp2[ng2 * ix * 2 + iy * 2] = matrixTip[ix][iy][0];
									temp2[ng2 * ix * 2 + iy * 2 + 1] = matrixTip[ix][iy][1];

								}
							}
							FFTWComplex fft = new FFTWComplex();
							fft.twoDimensional(ng2, ng1, temp, ftout,
									FFTW.BACKWARD);
							fft.twoDimensional(ng2, ng1, temp2, ftout2,
									FFTW.BACKWARD);

							double spinFactor = 2;
							double constant_I = -2 * Math.sqrt(Math.PI) / sigma
									* Math.pow(surfaceAU, 2)
									/ Math.pow(volumeAU, 2) * spinFactor;
							for (int ix = 0; ix < ng1; ix++) {
								for (int iy = 0; iy < ng2; iy++) {
									tersoffHamannMatrixSurface[ix][iy] = Math
											.pow(ftout[ng2 * ix * 2 + iy * 2],
													2)
											+ Math.pow(ftout[ng2 * ix * 2 + iy
													* 2 + 1], 2);
									tersoffHamannMatrixTip[ix][iy] = Math.pow(
											ftout2[ng2 * ix * 2 + iy * 2], 2)
											+ Math.pow(ftout2[ng2 * ix * 2 + iy
													* 2 + 1], 2);
									tersoffHamannCurrent[ix][iy][iz] += tersoffHamannMatrixSurface[ix][iy]
											/ (2 * nrKPoints - 1)
											/ volumeAU
											* Math.exp(-Math
													.pow((eigenEnergiesAU[b][0]
															- voltageAU - eFermiAU)
															/ sigma, 2))
											* constant_I;
									tersoffHamannSurface[ix][iy][iz] += tersoffHamannMatrixSurface[ix][iy]
											/ (2 * nrKPoints - 1) / volumeAU;
									tersoffHamannTip[ix][iy][iz] += tersoffHamannMatrixTip[ix][iy]
											/ (2 * nrKPoints - 1) / volumeAU;

								}
							}
							if (debug > 2 && iz == 1) {
								System.out.println("Current [0][0][0]: "
										+ tersoffHamannCurrent[0][0][0]);
								System.out.println("Current [2][2][0]: "
										+ tersoffHamannCurrent[2][2][0]);
								System.out.println("surface [0][0][1]: "
										+ tersoffHamannSurface[0][0][1]);
								System.out.println("surface [2][2][0]: "
										+ tersoffHamannSurface[2][2][0]);
							}
						}

					} // bands
					double percentageDone = (double) iz
							/ (double) zSamplingPoints * 100;

					if (percentageDone > progressPercentageThreshold) {
						System.out.print((int) percentageDone / 10);
						progressPercentageThreshold += 10;
					}

				} // zstep
				System.out.println(" done.");

				// give siesta File output
				String bardeenFileName = pathName + "bardeen-k" + (kp + 1)
						+ "spin" + (sp + 1) + "surfaceZ" + zIndices.surface
						+ "tipZ" + zIndices.tip + "voltage"
						+ df2.format(voltageAU*27.2116) + ".siesta";
				//writeSiestaFile(bardeenCurrent, bardeenFileName);

				String surfaceTHFileName = pathName + "surfaceTH-k" + (kp + 1)
						+ "spin" + (sp + 1) + "surfaceZ" + zIndices.surface
						+ "tipZ" + zIndices.tip + "voltage"
						+ df2.format(voltageAU*27.2116) + ".siesta";
				writeSiestaFile(tersoffHamannSurface, surfaceTHFileName);

				String tipTHFileName = pathName + "tipTH-k" + (kp + 1) + "spin"
						+ (sp + 1) + "surfaceZ" + zIndices.surface + "tipZ"
						+ zIndices.tip + "voltage" + df2.format(voltageAU*27.2116)
						+ ".siesta";
				//writeSiestaFile(tersoffHamannTip, tipTHFileName);

			} // spin
		} // kpoints

		if (mapping != null)
			mapping.close();
		try {
			wcF.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		// printToFile(zValuesAngst, tersoffHamannCurrent, "thc.dat");
		//printToFile(zValuesAngst, tersoffHamannSurface, "ths.dat");
		//printToFile(zValuesAngst, bardeenCurrent, "thb.dat");

	}

	private void calculateCoefficients(String coefficientsFileName) {
		System.out.print("calculating coefficients ");

		int progressPercentageThreshold = 0;
		try {
			DataOutputStream outCoefficients = new DataOutputStream(
					new FileOutputStream(coefficientsFileName));

			for (int b = 0; b < bands; b++) {
				readCoefficients(recordCounter);
				recordCounter++;

				// double[] temp2 = new double[2 * mapping.ngz];
				// double[] ftout = new double[2 * mapping.ngz];
				double[] temp2 = new double[2 * ng3];
				double[] ftout = new double[2 * ng3];

				for (int ix = 0; ix < ng1; ix++) {
					for (int iy = 0; iy < ng2; iy++) {
						for (int iz = 0; iz < ng3; iz++) {
							temp2[iz * 2] = coefficientMatrix[ix][iy][iz][0];
							temp2[iz * 2 + 1] = coefficientMatrix[ix][iy][iz][1];
						}
						FFTWComplex fft = new FFTWComplex();
						ftout = fft.oneDimensional(temp2, ftout, FFTW.BACKWARD);

						for (int iz = 0; iz < ng3; iz++) {
							coefficientMatrixReversed[ix][iy][iz][0] = (float) ftout[iz * 2];
							coefficientMatrixReversed[ix][iy][iz][1] = (float) ftout[iz * 2 + 1];
						}
					}
				}

				// System.out.println("Surface Z: "+surfaceZ);
				// System.out.println("Tip Z: "+tipZ);

				// System.out.println(coefficientMatrixReversed[0][0][surfaceZ]);

				// for (int ix = 0; ix < mapping.ngx; ix++) {
				// for (int iy = 0; iy < mapping.ngy; iy++) {
				for (int ix = 0; ix < ng1; ix++) {
					for (int iy = 0; iy < ng2; iy++) {
						tipCoefficientsPerBand[ix][iy][0][b] = coefficientMatrixReversed[ix][iy][zIndices.tip][0];
						tipCoefficientsPerBand[ix][iy][1][b] = coefficientMatrixReversed[ix][iy][zIndices.tip][1];
						surfaceCoefficientsPerBand[ix][iy][0][b] = 1
								/ (float) Math.sqrt(2)
								* coefficientMatrixReversed[ix][iy][zIndices.surface][0];
						surfaceCoefficientsPerBand[ix][iy][1][b] = 1
								/ (float) Math.sqrt(2)
								* coefficientMatrixReversed[ix][iy][zIndices.surface][1];
						outCoefficients
								.writeFloat(tipCoefficientsPerBand[ix][iy][0][b]);
						outCoefficients
								.writeFloat(tipCoefficientsPerBand[ix][iy][1][b]);
						outCoefficients
								.writeFloat(surfaceCoefficientsPerBand[ix][iy][0][b]);
						outCoefficients
								.writeFloat(surfaceCoefficientsPerBand[ix][iy][1][b]);

					}
				}

				double percentageDone = (double) b / (double) bands * 100;

				if (percentageDone > progressPercentageThreshold) {
					System.out.print((int) percentageDone / 10);
					progressPercentageThreshold += 10;
				}

			} // bands
			if (outCoefficients != null) {
				outCoefficients.close();
			}
		} catch (FileNotFoundException e) {
			System.out.println(e.toString());
		} catch (IOException e) {
			System.out.println(e.toString());
		}
		System.out.println(" done.");
		// writeCoefficientMatrix(0, 0);
		// System.out.println(0.85-2/latticeVectors[2][2]/0.529);

	}

	private void readCoefficientsFromFile(String coefficientsFileName) {
		System.out.println("Reading coefficients from file: "
				+ coefficientsFileName);
		System.out
				.println("If you want to calculate the coefficients from scratch delete the file: "
						+ coefficientsFileName);
		System.out.print("Reading coefficients from file ");

		int progressPercentageThreshold = 10;

		try {
			DataInputStream coefficientsFile = new DataInputStream(
					new FileInputStream(coefficientsFileName));
			// System.out.println("ok !");
			for (int b = 0; b < bands; b++) {
				recordCounter++;
				for (int ix = 0; ix < ng1; ix++) {
					for (int iy = 0; iy < ng2; iy++) {
						tipCoefficientsPerBand[ix][iy][0][b] = coefficientsFile
								.readFloat();
						tipCoefficientsPerBand[ix][iy][1][b] = coefficientsFile
								.readFloat();
						surfaceCoefficientsPerBand[ix][iy][0][b] = coefficientsFile
								.readFloat();
						surfaceCoefficientsPerBand[ix][iy][1][b] = coefficientsFile
								.readFloat();
						if (ix == 0) {
							// System.out.println(tipCoefficientsPerBand[ix][iy][0][b]);
							// System.out.println(tipCoefficientsPerBand[ix][iy][1][b]);
						}
					}
				}
				double percentageDone = (double) b / (double) bands * 100;
				if (percentageDone > progressPercentageThreshold) {
					System.out.print((int) percentageDone / 10);
					progressPercentageThreshold += 10;
				}
			}
			if (coefficientsFile != null)
				coefficientsFile.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		System.out.println(" done.");

	}

	private void checkForCorrectNumberOfPlanewaves(int calculatedNrPW) {
		if (calculatedNrPW != nrPlaneWaves) {
			System.err.println("calculated Planewaves: " + calculatedNrPW);
			System.err.println("read Planewaves: " + nrPlaneWaves);
			System.exit(1);
		}
	}

	public void calculateCurrents(double voltageAU) {
		for (int kp = 0; kp < (int) nrKPoints; kp++) {
			for (int sp = 0; sp < 1; sp++) {

				readEnergies(recordCounter);
				// System.out.println(this.toString());
				recordCounter++;

				// test if gamma-only version or other problems with nr. of plane waves
				// and calculate indices of plane waves !
				int cPlaneWaves = calculatePlanewaves();
				if (cPlaneWaves != nrPlaneWaves) {
					System.err.println("calculated Planewaves: " + cPlaneWaves);
					System.err.println("read Planewaves: " + nrPlaneWaves);
					System.exit(1);
				}
				nrPlaneWaves = cPlaneWaves;

				int tipZ = (int) ((zTipDirect - zDistance
						/ latticeVectorsAngst[2][2]) * ng3) - 1;
				int surfaceZ = (int) ((zHighestDirect + zDistance
						/ latticeVectorsAngst[2][2]) * ng3) - 1;
				
				System.out.println("surface z: "+surfaceZ);
				System.out.println("k-Point: " + (kp + 1)
						+ ", Spin-component: " + (sp + 1));
				String coefficientsFileName = pathName + "coefficients-k"
						+ (kp + 1) + "spin" + (sp + 1) + "surfaceZ" + surfaceZ
						+ "tipZ" + tipZ;
				if (currentCoefficients.compareTo(coefficientsFileName)==0){
					System.out.println("Correct coefficients already loaded.");
				}
				else if ((new File(coefficientsFileName)).isFile()) {
					currentCoefficients = coefficientsFileName;
					initializeArraysNew();
					System.out.println("Reading coefficients from file: "
							+ coefficientsFileName);
					System.out
							.println("If you want to calculate the coefficients from scratch delete the file: "
									+ coefficientsFileName);
					System.out.print("Reading coefficients from file ");
					int progressPercentageThreshold = 10;
					try {
						DataInputStream coefficientsFile = new DataInputStream(
								new FileInputStream(coefficientsFileName));
						// System.out.println("ok !");
						for (int b = 0; b < bands; b++) {
							recordCounter++;
							for (int ix = 0; ix < ng1; ix++) {
								for (int iy = 0; iy < ng2; iy++) {
									tipCoefficientsPerBand[ix][iy][0][b] = coefficientsFile
											.readFloat();
									tipCoefficientsPerBand[ix][iy][1][b] = coefficientsFile
											.readFloat();
									surfaceCoefficientsPerBand[ix][iy][0][b] = coefficientsFile
											.readFloat();
									surfaceCoefficientsPerBand[ix][iy][1][b] = coefficientsFile
											.readFloat();
									if (ix == 0) {
										// System.out.println(tipCoefficientsPerBand[ix][iy][0][b]);
										// System.out.println(tipCoefficientsPerBand[ix][iy][1][b]);
									}
								}
							}
							double percentageDone = (double) b / (double) bands
									* 100;
							if (percentageDone > progressPercentageThreshold) {
								System.out.print((int) percentageDone / 10);
								progressPercentageThreshold += 10;
							}
						}
						if (coefficientsFile != null)
							coefficientsFile.close();
					} catch (FileNotFoundException e) {
						e.printStackTrace();
					} catch (IOException e) {
						e.printStackTrace();
					}
					System.out.println(" done.");
				} else {
					initializeArraysNew();
					System.out.print("calculating coefficients ");

					int progressPercentageThreshold = 0;
					try {
						DataOutputStream outCoefficients = new DataOutputStream(
								new FileOutputStream(coefficientsFileName));

						for (int b = 0; b < bands; b++) {
							readCoefficients(recordCounter);
							recordCounter++;

							// double[] temp2 = new double[2 * mapping.ngz];
							// double[] ftout = new double[2 * mapping.ngz];
							double[] temp2 = new double[2 * ng3];
							double[] ftout = new double[2 * ng3];

							for (int ix = 0; ix < ng1; ix++) {
								for (int iy = 0; iy < ng2; iy++) {
									for (int iz = 0; iz < ng3; iz++) {
										temp2[iz * 2] = coefficientMatrix[ix][iy][iz][0];
										temp2[iz * 2 + 1] = coefficientMatrix[ix][iy][iz][1];
									}
									FFTWComplex fft = new FFTWComplex();
									ftout = fft.oneDimensional(temp2, ftout,
											FFTW.BACKWARD);

									for (int iz = 0; iz < ng3; iz++) {
										coefficientMatrixReversed[ix][iy][iz][0] = (float) ftout[iz * 2];
										coefficientMatrixReversed[ix][iy][iz][1] = (float) ftout[iz * 2 + 1];
									}
								}
							}

							// System.out.println("Surface Z: "+surfaceZ);
							// System.out.println("Tip Z: "+tipZ);

							// System.out.println(coefficientMatrixReversed[0][0][surfaceZ]);

							// for (int ix = 0; ix < mapping.ngx; ix++) {
							// for (int iy = 0; iy < mapping.ngy; iy++) {
							for (int ix = 0; ix < ng1; ix++) {
								for (int iy = 0; iy < ng2; iy++) {
									tipCoefficientsPerBand[ix][iy][0][b] = coefficientMatrixReversed[ix][iy][tipZ][0];
									tipCoefficientsPerBand[ix][iy][1][b] = coefficientMatrixReversed[ix][iy][tipZ][1];
									surfaceCoefficientsPerBand[ix][iy][0][b] = 1
											/ (float) Math.sqrt(2)
											* coefficientMatrixReversed[ix][iy][surfaceZ][0];
									surfaceCoefficientsPerBand[ix][iy][1][b] = 1
											/ (float) Math.sqrt(2)
											* coefficientMatrixReversed[ix][iy][surfaceZ][1];
									outCoefficients
											.writeFloat(tipCoefficientsPerBand[ix][iy][0][b]);
									outCoefficients
											.writeFloat(tipCoefficientsPerBand[ix][iy][1][b]);
									outCoefficients
											.writeFloat(surfaceCoefficientsPerBand[ix][iy][0][b]);
									outCoefficients
											.writeFloat(surfaceCoefficientsPerBand[ix][iy][1][b]);

								}
							}

							double percentageDone = (double) b / (double) bands
									* 100;

							if (percentageDone > progressPercentageThreshold) {
								System.out.print((int) percentageDone / 10);
								progressPercentageThreshold += 10;
							}

						} // bands
						if (outCoefficients != null) {
							outCoefficients.close();
						}
					} catch (FileNotFoundException e) {
						System.out.println(e.toString());
					} catch (IOException e) {
						System.out.println(e.toString());
					}
					System.out.println(" done.");
					// writeCoefficientMatrix(0, 0);
					// System.out.println(0.85-2/latticeVectors[2][2]/0.529);
				}

				System.out.print("calculating currents ");

				int progressPercentageThreshold = 10;

				double zStepAU = zMaxAU / (zStepsNr - 1);

				for (int ix = 0; ix < ng1; ix++) {
					for (int iy = 0; iy < ng2; iy++) {
						for (int iz = 0; iz < zSamplingPoints; iz++) {

							
							tersoffHamannCurrent[ix][iy][iz] = 0;
							tersoffHamannSurface[ix][iy][iz] = 0;
							tersoffHamannTip[ix][iy][iz] = 0;
						}
					}
				}

				
				for (int iz = 0; iz < zSamplingPoints; iz++) {
					double zAU = zStepAU * iz;
					zValuesAngst[iz] = zAU * 0.529 + zMinAngst;

					/* ******************************************************
					 * Bardeen
					 * ******************************************************
					 */
					//calculateBardeen(iz,zAU);
					
					/* ******************************************************
					 * Tersoff Hamann
					 * ******************************************************
					 */
					
					for (int b = 0; b < bands; b++) {
						
						if ((voltageAU < 0 && eigenEnergiesAU[b][0] < eFermiAU && eigenEnergiesAU[b][0] > eFermiAU+voltageAU)
							|| (voltageAU > 0 && eigenEnergiesAU[b][0] > eFermiAU && eigenEnergiesAU[b][0] < eFermiAU+voltageAU)
						) {
							// if (iz==0) System.out.println("band: "+b);
							for (int ix = 0; ix < ng1; ix++) {
								for (int iy = 0; iy < ng2; iy++) {
									surfaceCoefficients[ix][iy][0] = surfaceCoefficientsPerBand[ix][iy][0][b];
									surfaceCoefficients[ix][iy][1] = surfaceCoefficientsPerBand[ix][iy][1][b];
									tipCoefficients[ix][iy][0] = tipCoefficientsPerBand[ix][iy][0][b];
									tipCoefficients[ix][iy][1] = tipCoefficientsPerBand[ix][iy][1][b];
								}
							}

							for (int ix = 0; ix < ng1; ix++) {
								for (int iy = 0; iy < ng2; iy++) {
									// System.out.println(b+" "+ix+" "+iy);
									matrixSurface[ix][iy][0] = surfaceCoefficients[ix][iy][0]
											* Math.exp(-decay.matrix[ix][iy]
													* zAU);
									matrixSurface[ix][iy][1] = surfaceCoefficients[ix][iy][1]
											* Math.exp(-decay.matrix[ix][iy]
													* zAU);
									matrixTip[ix][iy][0] = tipCoefficients[ix][iy][0]
											* Math.exp(-decay.matrix[ix][iy]
													* zAU);
									matrixTip[ix][iy][1] = tipCoefficients[ix][iy][1]
											* Math.exp(-decay.matrix[ix][iy]
													* zAU);

								}
							}

							double[] temp = new double[2 * ng2 * ng1];
							double[] temp2 = new double[2 * ng2 * ng1];
							double[] ftout = new double[2 * ng2 * ng1];
							double[] ftout2 = new double[2 * ng2 * ng1];

							for (int ix = 0; ix < ng1; ix++) {
								for (int iy = 0; iy < ng2; iy++) {
									temp[ng2 * ix * 2 + iy * 2] = matrixSurface[ix][iy][0];
									temp[ng2 * ix * 2 + iy * 2 + 1] = matrixSurface[ix][iy][1];
									temp2[ng2 * ix * 2 + iy * 2] = matrixTip[ix][iy][0];
									temp2[ng2 * ix * 2 + iy * 2 + 1] = matrixTip[ix][iy][1];

								}
							}
							FFTWComplex fft = new FFTWComplex();
							fft.twoDimensional(ng2, ng1, temp, ftout,
									FFTW.BACKWARD);
							fft.twoDimensional(ng2, ng1, temp2, ftout2,
									FFTW.BACKWARD);

							double spinFactor = 2;
							double constant_I = -2 * Math.sqrt(Math.PI) / sigma
									* Math.pow(surfaceAU, 2)
									/ Math.pow(volumeAU, 2) * spinFactor;
							if (voltageAU>0) constant_I = - constant_I; 
							for (int ix = 0; ix < ng1; ix++) {
								for (int iy = 0; iy < ng2; iy++) {
									tersoffHamannMatrixSurface[ix][iy] = Math
											.pow(ftout[ng2 * ix * 2 + iy * 2],
													2)
											+ Math.pow(ftout[ng2 * ix * 2 + iy
													* 2 + 1], 2);
									tersoffHamannMatrixTip[ix][iy] = Math.pow(
											ftout2[ng2 * ix * 2 + iy * 2], 2)
											+ Math.pow(ftout2[ng2 * ix * 2 + iy
													* 2 + 1], 2);
									tersoffHamannCurrent[ix][iy][iz] += tersoffHamannMatrixSurface[ix][iy]
											/ (2 * nrKPoints - 1)
											/ volumeAU
											* Math.exp(-Math
													.pow((eigenEnergiesAU[b][0]
															- voltageAU - eFermiAU)
															/ sigma, 2))
											* constant_I;
									tersoffHamannSurface[ix][iy][iz] += tersoffHamannMatrixSurface[ix][iy]
											/ (2 * nrKPoints - 1) / volumeAU;
									tersoffHamannTip[ix][iy][iz] += tersoffHamannMatrixTip[ix][iy]
											/ (2 * nrKPoints - 1) / volumeAU;

								}
							}
							if (debug > 2 && iz == 1) {
								System.out.println("Current [0][0][0]: "
										+ tersoffHamannCurrent[0][0][0]);
								System.out.println("Current [2][2][0]: "
										+ tersoffHamannCurrent[2][2][0]);
								System.out.println("surface [0][0][1]: "
										+ tersoffHamannSurface[0][0][1]);
								System.out.println("surface [2][2][0]: "
										+ tersoffHamannSurface[2][2][0]);
							}
						}

					} // bands
					double percentageDone = (double) iz
							/ (double) zSamplingPoints * 100;

					if (percentageDone > progressPercentageThreshold) {
						System.out.print((int) percentageDone / 10);
						progressPercentageThreshold += 10;
					}

				} // zstep
				System.out.println(" done.");

				// give siesta File output
				String bardeenFileName = pathName + "bardeen-k" + (kp + 1)
						+ "spin" + (sp + 1) + "surfaceZ" + surfaceZ + "tipZ"
						+ tipZ + "voltage" + df2.format(voltageAU*27.2116) + ".siesta";
				//writeSiestaFile(bardeenCurrent, bardeenFileName);

				String surfaceTHFileName = pathName + "surfaceTH-k" + (kp + 1)
						+ "spin" + (sp + 1) + "surfaceZ" + surfaceZ + "tipZ"
						+ tipZ + "voltage" + df2.format(voltageAU*27.2116) + ".siesta";
				writeSiestaFile(tersoffHamannSurface, surfaceTHFileName);

				String tipTHFileName = pathName + "tipTH-k" + (kp + 1) + "spin"
						+ (sp + 1) + "surfaceZ" + surfaceZ + "tipZ" + tipZ
						+ "voltage" + df2.format(voltageAU*27.2116) + ".siesta";
				//writeSiestaFile(tersoffHamannTip, tipTHFileName);
				
				
				currentMapToPng(tersoffHamannSurface, Math.pow(10, -7), "s7_voltage"
						+ df2.format(voltageAU*27.2116) + "surfaceZ" + surfaceZ
						+ ".png", 3);
				
				/*printCurrentMap(tersoffHamannSurface, 3.0, 
						"height3_voltage" + df2.format(voltageAU*27.2116) + "surfaceZ" + surfaceZ
						+ ".dat");*/

			} // spin
		} // kpoints

		if (mapping != null)
			mapping.close();
		try {
			wcF.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		// printToFile(zValuesAngst, tersoffHamannCurrent, "thc.dat");
		printToFile(zValuesAngst, tersoffHamannSurface, "ths.dat");
		printToFile(zValuesAngst, bardeenCurrent, "thb.dat");

	}

	private void writeSiestaFile(double[][][] current2, String filePath) {
		try {
			LERandomAccessFile outSiesta = new LERandomAccessFile(filePath,
					"rw");
			outSiesta.writeInt(72);
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					if (i == 2 && j == 2) {
						outSiesta.writeDouble(zMaxAngst);
					} else
						outSiesta
								.writeDouble((float) latticeVectorsAngst[i][j]);
				}
			}
			outSiesta.writeInt(72);

			outSiesta.writeInt(16);
			outSiesta.writeInt(ng1);
			outSiesta.writeInt(ng2);
			outSiesta.writeInt(zSamplingPoints);
			outSiesta.writeInt(1);
			outSiesta.writeInt(16);

			for (int iz = 0; iz < zSamplingPoints; iz++) {
				for (int iy = 0; iy < ng2; iy++) {
					outSiesta.writeInt(ng1 * 4);
					for (int ix = 0; ix < ng1; ix++) {
						outSiesta.writeFloat((float) Math
								.abs(current2[ix][iy][iz]));
					}
					outSiesta.writeInt(ng1 * 4);
				}
			}

			if (outSiesta != null) {
				outSiesta.close();
			}
		} catch (FileNotFoundException e) {
			System.out.println(e.toString());
		} catch (IOException e) {
			System.out.println(e.toString());
		}

	}

	public void readTersoffHamannValuesFromTextFile(double[][][] map,
			String fileName) {
		try {
			BufferedReader in = new BufferedReader(new FileReader(fileName));

			/*
			 * for (int ix = 0; ix < mapping.ngx; ix++) { for (int iy = 0; iy <
			 * mapping.ngy; iy++) { for (int iz = 0; iz < zSamplingPoints; iz++)
			 * { String[] t1 = in.readLine().split("\t"); //
			 * System.out.println(t1); zValuesAngst[iz] =
			 * Double.parseDouble(t1[0]); map[ix][iy][iz] =
			 * Double.parseDouble(t1[1]);
			 * 
			 * } } }
			 */
			for (int ix = 0; ix < ng1; ix++) {
				for (int iy = 0; iy < ng2; iy++) {
					for (int iz = 0; iz < zSamplingPoints; iz++) {
						String[] t1 = in.readLine().split("\t");
						// System.out.println(t1);
						zValuesAngst[iz] = Double.parseDouble(t1[0]);
						map[ix][iy][iz] = Double.parseDouble(t1[1]);

					}
				}
			}
			if (in != null) {
				in.close();
			}
		} catch (FileNotFoundException e) {
			System.out.println(e.toString());
		} catch (IOException e) {
			System.out.println(e.toString());

		}

	}

	public double[] findMinimumAndMaximumOfCurrents(double[][][] currents) {
		double currentMax = 0;
		double currentMin = 1000;
		// for (int ix = 0; ix < mapping.ngx; ix++) {
		for (int ix = 0; ix < ng1; ix++) {
			// for (int iy = 0; iy < mapping.ngy; iy++) {
			for (int iy = 0; iy < ng2; iy++) {
				for (int iz = 0; iz < zSamplingPoints; iz++) {
					if (currents[ix][iy][iz] > currentMax)
						currentMax = currents[ix][iy][iz];
					if (currents[ix][iy][iz] < currentMin)
						currentMin = currents[ix][iy][iz];

				}
			}
		}
		return new double[] { currentMin, currentMax };
	}

	public void currentMapToPng(double[][][] currentMap, double isoCurrent,
			String outName, int scalingFactor) {

		BufferedImage isoCurrentImage = new BufferedImage(ng2 * scalingFactor,
				ng1 * scalingFactor, BufferedImage.TYPE_BYTE_GRAY);
		WritableRaster raster = isoCurrentImage.getRaster();

		double[][] heightMap = new double[currentMap.length][currentMap[0].length];
		double max = zMinAngst;
		double min = 20;
		int xmax = -1;
		int ymax = -1;

		for (int ix = 0; ix < ng1; ix++) {
			for (int iy = 0; iy < ng2; iy++) {
				heightMap[ix][iy] = zMinAngst;
				for (int iz = 0; iz < zSamplingPoints - 1; iz++) {
					// if (Math.abs(currentMap[ix][iy][iz]) <
					// Math.abs(isoCurrent)){
					// break;
					// }
					double current = Math.abs(currentMap[ix][iy][iz]);
					double current1 = Math.abs(currentMap[ix][iy][iz + 1]);

					if (current >= Math.abs(isoCurrent)
							&& current1 <= Math.abs(isoCurrent)) {

						heightMap[ix][iy] = zValuesAngst[iz]
								+ (current - isoCurrent)
								* (zValuesAngst[iz + 1] - zValuesAngst[iz])
								/ (current - current1);
						if (heightMap[ix][iy] > max) {
							max = heightMap[ix][iy];
							xmax = ix;
							ymax = iy;
						}
						if (heightMap[ix][iy] < min)
							min = heightMap[ix][iy];
						break;
					}

				}
			}
		}

		/*
		 * if (isoCurrent<0){ for (int ix = 0; ix < ng1; ix++) { for (int iy =
		 * 0; iy < ng2; iy++) { heightMap[ix][iy] = min + max -
		 * heightMap[ix][iy]; } } }
		 */

		System.out.println("current: " + isoCurrent + "\tmax: " + max
				+ "\tmin: " + min + "\txmax: " + xmax + "\tymax: " + ymax);
		for (int iy = 0; iy < ng2; iy++) {
			for (int ix = 0; ix < ng1; ix++) {
				int value = (int) Math.round((double) (heightMap[ix][iy] - min)
						/ (double) (max - min) * 255.0);
				//if (ix == 60 && iy == 70)
					//value = 0;
				for (int sx = 1; sx <= scalingFactor; sx++) {
					for (int sy = 1; sy <= scalingFactor; sy++) {
						raster.setSample(scalingFactor * iy + sy - 1,
								scalingFactor * ix + sx - 1, 0, value);
					}
				}

			}

		}

		try {
			File outputfile = new File(pathName + outName);
			ImageIO.write(isoCurrentImage, "png", outputfile);
		} catch (IOException e) {
			e.printStackTrace();
		}

	}
	
	public void printCurrentMap(double[][][] currentMap, double height,
			String outName) {
		int hZ = (int) ((zHighestDirect + height
				/ latticeVectorsAngst[2][2]) * ng3) - 1;
		System.out.println("hZ: " + hZ);
		
		
		
		
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(pathName
					+ outName));			
			for (int iy = 0; iy < ng2; iy++) {
				for (int ix = 0; ix < ng1; ix++) {
					out.write(""+currentMap[ix][iy][hZ]+" ");
				}
				out.newLine();
			}
		} catch (IOException e) {
			e.printStackTrace();
		}

	}
	
	
	public void calculateBardeen(int iz, double zAU){
		for (int tb = 0; tb < bands; tb++) {
			for (int sb = 0; sb < bands; sb++) {
				double spinFactor = 2;
				double constant_I = 2 * Math.sqrt(Math.PI) / sigma
						* Math.pow(surfaceAU, 2)
						/ Math.pow(volumeAU, 2) * spinFactor;
				if (voltage < 0)
					constant_I = -constant_I;

				if ((voltage < 0
						&& eigenEnergiesAU[tb][0] - eFermiAU > 0 && eigenEnergiesAU[sb][0]
						- eFermiAU < 0)
						|| (voltage > 0
								&& eigenEnergiesAU[tb][0]
										- eFermiAU < 0 && eigenEnergiesAU[sb][0]
								- eFermiAU > 0)) {
					double deltaE = Math
							.pow(((eigenEnergiesAU[tb][0] - eigenEnergiesAU[sb][0]) + voltageAU)
									/ sigma, 2);
					if (deltaE < 3.0) {
						// System.out.println("deltaE<3");
						for (int ix = 0; ix < ng1; ix++) {
							for (int iy = 0; iy < ng2; iy++) {
								surfaceCoefficients[ix][iy][0] = surfaceCoefficientsPerBand[ix][iy][0][sb];
								surfaceCoefficients[ix][iy][1] = surfaceCoefficientsPerBand[ix][iy][1][sb];
								tipCoefficients[ix][iy][0] = tipCoefficientsPerBand[ix][iy][0][tb];
								tipCoefficients[ix][iy][1] = tipCoefficientsPerBand[ix][iy][1][tb];
							}
						}
						for (int ix = 0; ix < ng1; ix++) {
							for (int iy = 0; iy < ng2; iy++) {
								matrixSurface[ix][iy][0] = tipCoefficients[ix][iy][0]
										* surfaceCoefficients[ix][iy][0]
										* decay.matrix[ix][iy]
										* Math.exp(-decay.matrix[ix][iy]
												* zAU);
								matrixSurface[ix][iy][1] = -tipCoefficients[ix][iy][1]
										* surfaceCoefficients[ix][iy][1]
										* decay.matrix[ix][iy]
										* Math.exp(-decay.matrix[ix][iy]
												* zAU);
								// if (matrixSurface[ix][iy][0]>0)
								// System.out.println("ix: "+ix+" iy: "+iy+" matrix: "+matrixSurface[ix][iy][0]);

							}
						}

						double[] temp = new double[2 * ng2 * ng1];
						double[] ftout = new double[2 * ng2 * ng1];

						for (int ix = 0; ix < ng1; ix++) {
							for (int iy = 0; iy < ng2; iy++) {
								temp[ng2 * ix * 2 + iy * 2] = matrixSurface[ix][iy][0];
								temp[ng2 * ix * 2 + iy * 2 + 1] = matrixSurface[ix][iy][1];
							}
						}
						FFTWComplex fft = new FFTWComplex();
						fft.twoDimensional(ng2, ng1, temp, ftout,
								FFTW.BACKWARD);

						for (int ix = 0; ix < ng1; ix++) {
							for (int iy = 0; iy < ng2; iy++) {
								bardeenMatrix[ix][iy] = Math
										.pow(ftout[ng2 * ix * 2
												+ iy * 2], 2)
										+ Math.pow(ftout[ng2 * ix
												* 2 + iy * 2 + 1],
												2);

								bardeenCurrent[ix][iy][iz] += bardeenMatrix[ix][iy]
										* constant_I
										* Math.exp(-deltaE)
										/ (2 * nrKPoints - 1);
								// if (bardeenCurrent[ix][iy][iz]>0)
								// System.out.println("cu: ix: "+ix+" iy: "+iy+" "+bardeenMatrix[ix][iy]);

							}
						}

					}// if deltaE
				}// if voltage

			}// bands
		}// bands

		
	}

	public void readHeadData(String wavecarName, String mappingscarName) {

		try {
			// mapping = new MappingscarIO(mappingscarName, pathName);
			wcF = new LERandomAccessFile(wavecarName, "r");
			// wcF = new RandomAccessFile(wavecarName, "r");
			// BufferedWriter out = new BufferedWriter(new FileWriter(pathName
			// + "out.dat"));

			wavecarRecordLength = (long) wcF.readDouble();
			// wavecarRecordLength *= 4; // sigh ... wegen byterecl flag bei der
			// VASP kompilierung
			// wavecarRecordLength = (long) wcF2.readDouble();

			
			nrSpin = (int) wcF.readDouble();
			
			version = wcF.readDouble();
			
			if ((int) version != 45200 && (int) version != 45210)
				throw new WavecarCorruptException();

			wcF.seek(wavecarRecordLength);

			nrKPoints = wcF.readDouble();
			kPoints = new double[3];
			bands = (int) wcF.readDouble();
			enMax = (int) wcF.readDouble();
			latticeVectorsAU = new double[3][3];
			latticeVectorsAngst = new double[3][3];
			reciprocalVectorsAngst = new double[3][3];
			reciprocalVectorsAU = new double[3][3];

			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					latticeVectorsAngst[i][j] = wcF.readDouble();
					latticeVectorsAU[i][j] = latticeVectorsAngst[i][j] / 0.529;
				}
			}
			
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

			reciprocalVectorsAngst[0][0] = 2 * Math.PI * productAngst12[0]
					/ (dotProduct(latticeVectorsAngst[0], productAngst12));
			reciprocalVectorsAngst[0][1] = 2 * Math.PI * productAngst12[1]
					/ (dotProduct(latticeVectorsAngst[0], productAngst12));
			reciprocalVectorsAngst[0][2] = 2 * Math.PI * productAngst12[2]
					/ (dotProduct(latticeVectorsAngst[0], productAngst12));
			reciprocalVectorsAU[0][0] = 2 * Math.PI * productAU12[0]
					/ (dotProduct(latticeVectorsAU[0], productAU12));
			reciprocalVectorsAU[0][1] = 2 * Math.PI * productAU12[1]
					/ (dotProduct(latticeVectorsAU[0], productAU12));
			reciprocalVectorsAU[0][2] = 2 * Math.PI * productAU12[2]
					/ (dotProduct(latticeVectorsAU[0], productAU12));
			reciprocalVectorsAngst[1][0] = 2 * Math.PI * productAngst20[0]
					/ (dotProduct(latticeVectorsAngst[0], productAngst12));
			reciprocalVectorsAngst[1][1] = 2 * Math.PI * productAngst20[1]
					/ (dotProduct(latticeVectorsAngst[0], productAngst12));
			reciprocalVectorsAngst[1][2] = 2 * Math.PI * productAngst20[2]
					/ (dotProduct(latticeVectorsAngst[0], productAngst12));
			reciprocalVectorsAU[1][0] = 2 * Math.PI * productAU20[0]
					/ (dotProduct(latticeVectorsAU[0], productAU12));
			reciprocalVectorsAU[1][1] = 2 * Math.PI * productAU20[1]
					/ (dotProduct(latticeVectorsAU[0], productAU12));
			reciprocalVectorsAU[1][2] = 2 * Math.PI * productAU20[2]
					/ (dotProduct(latticeVectorsAU[0], productAU12));
			reciprocalVectorsAngst[2][0] = 2 * Math.PI * productAngst01[0]
					/ (dotProduct(latticeVectorsAngst[0], productAngst12));
			reciprocalVectorsAngst[2][1] = 2 * Math.PI * productAngst01[1]
					/ (dotProduct(latticeVectorsAngst[0], productAngst12));
			reciprocalVectorsAngst[2][2] = 2 * Math.PI * productAngst01[2]
					/ (dotProduct(latticeVectorsAngst[0], productAngst12));
			reciprocalVectorsAU[2][0] = 2 * Math.PI * productAU01[0]
					/ (dotProduct(latticeVectorsAU[0], productAU12));
			reciprocalVectorsAU[2][1] = 2 * Math.PI * productAU01[1]
					/ (dotProduct(latticeVectorsAU[0], productAU12));
			reciprocalVectorsAU[2][2] = 2 * Math.PI * productAU01[2]
					/ (dotProduct(latticeVectorsAU[0], productAU12));

			volumeAngst = (dotProduct(latticeVectorsAngst[0], productAngst12));
			volumeAU = (dotProduct(latticeVectorsAU[0], productAU12));
			surfaceAU = (Math
					.sqrt(Math.pow(productAU01[0], 2)
							+ Math.pow(productAU01[1], 2)
							+ Math.pow(productAU01[2], 2)));
			
			double phi12 = Math
					.acos(dotProduct(reciprocalVectorsAngst[0],
							reciprocalVectorsAngst[1])
							/ (vectorLength(reciprocalVectorsAngst[0]) * vectorLength(reciprocalVectorsAngst[1])));
			
			double[] Xproductb01 = vectorialProduct(reciprocalVectorsAngst[0],
					reciprocalVectorsAngst[1]);
			double sinphi012 = (dotProduct(reciprocalVectorsAngst[2],
					Xproductb01))
					/ (vectorLength(Xproductb01) * vectorLength(reciprocalVectorsAngst[2]));

			nb1Max = (int) (Math.sqrt(enMax * c) / (vectorLength(reciprocalVectorsAngst[0]) * Math
					.abs(Math.sin(phi12)))) + 1;
			nb2Max = (int) (Math.sqrt(enMax * c) / (vectorLength(reciprocalVectorsAngst[1]) * Math
					.abs(Math.sin(phi12)))) + 1;
			nb3Max = (int) (Math.sqrt(enMax * c) / (vectorLength(reciprocalVectorsAngst[2]) * Math
					.abs(sinphi012))) + 1;
			nb1Max = (int) (Math.sqrt(enMax * c) / (vectorLength(reciprocalVectorsAngst[0]) * Math
					.abs(Math.sin(phi12))));
			nb2Max = (int) (Math.sqrt(enMax * c) / (vectorLength(reciprocalVectorsAngst[1]) * Math
					.abs(Math.sin(phi12))));
			nb3Max = (int) (Math.sqrt(enMax * c) / (vectorLength(reciprocalVectorsAngst[2]) * Math
					.abs(sinphi012)));
			ng1 = getAllowedGridSize(nb1Max);
			ng2 = getAllowedGridSize(nb2Max);
			ng3 = getAllowedGridSize(nb3Max);
			// ng1 = 112; ng2 = 112; ng3 = 120; //CuPc

			//System.out.println(this.toStringSmall());

			decay = new Decay(ng1, ng2, reciprocalVectorsAU);
			if (debug > 0)
				printMatrixToFile(decay.matrix, decay.max, "decay.png");

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

	public int getAllowedGridSize(int dimension) {
		int trial = dimension * 4 * 3 / 4 + 1; // corresponds to
		// trial = dimension * 4 * 3 / 2; // corresponds to
		// Precision=Normal/Low
		
		if (trial < allowedGridSizes[0]) {
			System.out.println("dimension too low, update allowedGridSizes");
			System.exit(1);
		}
		if (trial > allowedGridSizes[allowedGridSizes.length - 1]) {
			System.out.println("dimension to high, update allowedGridSizes");
			System.exit(1);
		}

		for (int i = 1; i < allowedGridSizes.length; i++) {
			if (trial <= allowedGridSizes[i]) {
				return allowedGridSizes[i];
			}
		}

		return trial; // this should never happen

	}

	public void readEnergies(long recordCounter) {
		// read energies, kpoints, occupation
		try {
			wcF.seek(recordCounter * wavecarRecordLength);

			nrPlaneWaves = (int) wcF.readDouble();
			System.out.println("nrPlaneWaves: " + nrPlaneWaves);
			// read kPoints they will not be used for the calculation
			kPoints[0] = wcF.readDouble();
			kPoints[1] = wcF.readDouble();
			kPoints[2] = wcF.readDouble();

			eigenEnergiesAU = new double[bands][2];
			eigenEnergiesEV = new double[bands][2];

			cReal = new float[nrPlaneWaves];
			cImag = new float[nrPlaneWaves];

			ferTotal = new double[bands];
			for (int b = 0; b < bands; b++) {
				eigenEnergiesEV[b][0] = wcF.readDouble();
				eigenEnergiesAU[b][0] = eigenEnergiesEV[b][0] / 27.2116;
				eigenEnergiesEV[b][1] = wcF.readDouble();
				eigenEnergiesAU[b][1] = eigenEnergiesEV[b][1] / 27.2116;
				ferTotal[b] = wcF.readDouble();
			}
		} catch (FileNotFoundException e) {
			System.out.println(e.toString());
		} catch (IOException e) {
			System.out.println(e.toString());
		}
	}

	public float[][][][] getSurfaceCoefficients() {
		return surfaceCoefficientsPerBand;
	}

	public float[][][][] getTipCoefficients() {
		return tipCoefficientsPerBand;
	}

	public void initializeSomeArrays() {
		tipCoefficientsPerBand = new float[mapping.ngx][mapping.ngy][2][bands];
		surfaceCoefficientsPerBand = new float[mapping.ngx][mapping.ngy][2][bands];
		zValuesAngst = new double[zStepsNr];
		tipCoefficients = new float[mapping.ngx][mapping.ngy][2];
		surfaceCoefficients = new float[mapping.ngx][mapping.ngy][2];

		matrixSurface = new double[mapping.ngx][mapping.ngy][2];
		tipTersoffHamann = new double[mapping.ngx][mapping.ngy];
		surfaceTersoffHamann = new double[mapping.ngx][mapping.ngy];
		tersoffHamannMatrixSurface = new double[mapping.ngx][mapping.ngy];
		bardeenMatrix = new double[mapping.ngx][mapping.ngy];
		bardeenCurrent = new double[mapping.ngx][mapping.ngy][zStepsNr];
		tersoffHamannCurrent = new double[mapping.ngx][mapping.ngy][zStepsNr];
		tersoffHamannSurface = new double[mapping.ngx][mapping.ngy][zStepsNr];

	}

	public void initializeArraysNew() {
		coefficientMatrix = new float[ng1][ng2][ng3][2];
		coefficientMatrixReversed = new float[ng1][ng2][ng3][2];
		tipCoefficientsPerBand = new float[ng1][ng2][2][bands];
		tipCoefficients = new float[ng1][ng2][2];
		surfaceCoefficientsPerBand = new float[ng1][ng2][2][bands];
		surfaceCoefficients = new float[ng1][ng2][2];
		matrixSurface = new double[ng1][ng2][2];
		matrixTip = new double[ng1][ng2][2];
		tipTersoffHamann = new double[ng1][ng2];
		surfaceTersoffHamann = new double[ng1][ng2];
		tersoffHamannMatrixSurface = new double[ng1][ng2];
		tersoffHamannMatrixTip = new double[ng1][ng2];
		tersoffHamannCurrent = new double[ng1][ng2][zStepsNr];
		tersoffHamannSurface = new double[ng1][ng2][zStepsNr];
		tersoffHamannTip = new double[ng1][ng2][zStepsNr];
		bardeenMatrix = new double[ng1][ng2];
		bardeenCurrent = new double[ng1][ng2][zStepsNr];

		zValuesAngst = new double[zStepsNr];
	}

	public void readCoefficients(long recordCounter) {
		try {
			wcF.seek(recordCounter * wavecarRecordLength);
			for (int npl = 0; npl < nrPlaneWaves; npl++) {
				float cReal = 0;
				float cImag = 0;
				if (version == 45200) {
					cReal = (float) wcF.readFloat();
					cImag = (float) wcF.readFloat();
				} else if (version == 45210) {
					cReal = (float) wcF.readDouble();
					cImag = (float) wcF.readDouble();
				}
				coefficientMatrix[indices[0][npl]][indices[1][npl]][indices[2][npl]][0] = cReal;
				coefficientMatrix[indices[0][npl]][indices[1][npl]][indices[2][npl]][1] = cImag;
				// out.write(cw[npl].real()+"\t"+cw[npl].imag()+"\n");
			}
		} catch (FileNotFoundException e) {
			System.out.println(e.toString());
		} catch (IOException e) {
			System.out.println(e.toString());
		}
	}

	public void writeCoefficientMatrix(int x, int y) {
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(pathName
					+ "out.dat"));
			BufferedWriter outReversed = new BufferedWriter(new FileWriter(
					pathName + "outReversed.dat"));

			for (int i = 0; i < mapping.ngz; i++) {
				out.write("" + coefficientMatrix[x][y][i] + "\t"
						+ coefficientMatrix[x][y][2 * i + 1]);
				out.newLine();
			}
			for (int i = 0; i < mapping.ngz; i++) {
				outReversed.write("" + coefficientMatrixReversed[x][y][2 * i]
						+ "\t" + coefficientMatrixReversed[x][y][2 * i + 1]);
				outReversed.newLine();
			}

			if (out != null) {
				out.close();
			}
			if (outReversed != null) {
				outReversed.close();
			}
		} catch (IOException e) {
			e.printStackTrace();
		} catch (NullPointerException e) {
			e.printStackTrace();
		}
	}

	public String toString() {
		String output = "Name: " + wavecarName + "\n";
		output += "version: " + version + "\n";
		output += "k-Points: " + nrKPoints + "\tbands: " + bands + "\n";
		output += "enMax: " + enMax + "\tspin: " + nrSpin + "\n";
		output += "WAVECAR Record-Length: " + wavecarRecordLength + "\n";

		output += "a1: " + printVector2(latticeVectorsAngst[0]) + "\n";
		output += "a2: " + printVector2(latticeVectorsAngst[1]) + "\n";
		output += "a3: " + printVector2(latticeVectorsAngst[2]) + "\n";
		output += "b1: " + printVector4(reciprocalVectorsAngst[0]) + "\n";
		output += "b2: " + printVector4(reciprocalVectorsAngst[1]) + "\n";
		output += "b3: " + printVector4(reciprocalVectorsAngst[2]) + "\n";
		output += "volume [Angstrom]: " + volumeAngst + "\n";
		output += "nb1Max: " + nb1Max + "\n";
		output += "nb2Max: " + nb2Max + "\n";
		output += "nb3Max: " + nb3Max + "\n";
		output += "ng1: " + ng1 + "\n";
		output += "ng2: " + ng2 + "\n";
		output += "ng3: " + ng3 + "\n";

		output += "Number of plane waves: " + nrPlaneWaves + "\tspin: "
				+ nrSpin + "\n";
		output += "K-Points: " + kPoints[0] + "\t" + kPoints[1] + "\t"
				+ latticeVectorsAU[2] + "\n";

		output += "Eigen Energies AU: " + eigenEnergiesAU[0][0] + "\t"
				+ eigenEnergiesAU[1][0] + "\t" + eigenEnergiesAU[2][0] + "\t"
				+ eigenEnergiesAU[3][0] + "\n";
		output += "Eigen Energies eV: " + eigenEnergiesEV[0][0] + "\t"
				+ eigenEnergiesEV[1][0] + "\t" + eigenEnergiesEV[2][0] + "\t"
				+ eigenEnergiesEV[3][0] + "\n";
		output += "FerTotal: " + ferTotal[0] + "\t" + ferTotal[bands - 1]
				+ "\n";
		// output += mapping.toString();

		return output;
	}

	public String toStringSmall() {
		String output = "Name: " + wavecarName + "\n";
		output += "version: " + version + "\n";
		output += "k-Points: " + nrKPoints + "\tbands: " + bands + "\n";
		output += "enMax: " + enMax + "\tspin: " + nrSpin + "\n";
		output += "WAVECAR Record-Length: " + wavecarRecordLength + "\n";

		output += "a1: " + printVector2(latticeVectorsAngst[0]) + "\n";
		output += "a2: " + printVector2(latticeVectorsAngst[1]) + "\n";
		output += "a3: " + printVector2(latticeVectorsAngst[2]) + "\n";
		output += "b1: " + printVector4(reciprocalVectorsAngst[0]) + "\n";
		output += "b2: " + printVector4(reciprocalVectorsAngst[1]) + "\n";
		output += "b3: " + printVector4(reciprocalVectorsAngst[2]) + "\n";
		output += "volume [Angstrom]: " + volumeAngst + "\n";
		output += "nb1Max: " + nb1Max + "\n";
		output += "nb2Max: " + nb2Max + "\n";
		output += "nb3Max: " + nb3Max + "\n";
		output += "ng1: " + ng1 + "\n";
		output += "ng2: " + ng2 + "\n";
		output += "ng3: " + ng3 + "\n";

		output += "Number of plane waves: " + nrPlaneWaves + "\tspin: "
				+ nrSpin + "\n";
		output += "K-Points: " + kPoints[0] + "\t" + kPoints[1] + "\t"
				+ latticeVectorsAU[2] + "\n";
		return output;
	}

	public void printMatrixToFile(double[][] matrix, double max, String fileName) {
		System.out.println(matrix.length);
		System.out.println(matrix[0].length);
		BufferedImage indexImage = new BufferedImage(matrix.length,
				matrix[0].length, BufferedImage.TYPE_BYTE_GRAY);
		WritableRaster raster = indexImage.getRaster();

		for (int ix = 0; ix < matrix.length; ix++) {
			for (int iy = 0; iy < matrix[0].length; iy++) {
				raster.setSample(ix, iy, 0,
						Math.round(matrix[ix][iy] / max * 255));
			}
		}
		try {
			File outputfile = new File(pathName + fileName);
			ImageIO.write(indexImage, "png", outputfile);
		} catch (IOException e) {
			e.printStackTrace();
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
		product = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
		return product;
	}

	public double vectorLength(double[] v1) {
		return Math.sqrt(v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]);
	}

	public String printVector2(double[] v1) {
		String out = "";
		out += df2.format(v1[0]) + "\t" + df2.format(v1[1]) + "\t"
				+ df2.format(v1[2]);
		return out;
	}

	public String printVector3(double[] v1) {
		String out = "";
		out += df3.format(v1[0]) + "\t" + df3.format(v1[1]) + "\t"
				+ df3.format(v1[2]);
		return out;
	}

	public String printVector4(double[] v1) {
		String out = "";
		out += df4.format(v1[0]) + "\t" + df4.format(v1[1]) + "\t"
				+ df4.format(v1[2]);
		return out;
	}

	class ZIndices {
		int tip;
		int surface;

		ZIndices() {
			tip = (int) ((zTipDirect - zDistance / latticeVectorsAngst[2][2]) * ng3) - 1;
			surface = (int) ((zHighestDirect + zDistance
					/ latticeVectorsAngst[2][2]) * ng3) - 1;
		}
	}

}

class WavecarCorruptException_2 extends Exception {
	private static final long serialVersionUID = -3107556003856126691L;

	public WavecarCorruptException_2() {
		super();
	}
}

class Decay_2 {
	int debug = 1;
	int ng1;
	int ng2;
	double phi = 4.7 / 27.2116;
	public double max;
	double[][] matrix;

	public Decay_2(int ng1, int ng2, double[][] vectorsAU) {
		matrix = new double[ng1][ng2];
		for (int iy = 0; iy < ng2; iy++) {
			int iyIndex = iy;
			if (iy > ng2 / 2)
				iyIndex = iy - ng2;
			for (int ix = 0; ix < ng1; ix++) {
				int ixIndex = ix;
				if (ix > ng1 / 2)
					ixIndex = ix - ng1;

				matrix[ix][iy] = Math.sqrt(2 * phi
						+ Math.pow(vectorsAU[0][0] * ixIndex, 2)
						+ Math.pow(vectorsAU[1][1] * iyIndex, 2));
				if (max < matrix[ix][iy])
					max = matrix[ix][iy];
				if (debug > 1) {
					if (ix < 2 && iy < 2) {
						System.out.println("ix: " + ix + ", iy: " + iy
								+ ", decay: " + matrix[ix][iy]);
						System.out.println("vectorsAU[0][0]: "
								+ vectorsAU[0][0] + "\tvectorsAU[0][0]: "
								+ vectorsAU[0][0]);
					}
					// if (ix<ng1/2+2 && ix>ng1/2-2 && iy<ng2/2+2 && iy>ng2/2-4
					// )
					// System.out.println("ix: "+ix+", iy: "+iy+", decay: "+matrix[ix][iy]);
				}
			}
		}
	}

}
