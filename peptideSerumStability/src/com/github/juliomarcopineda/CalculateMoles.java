package com.github.juliomarcopineda;
import java.util.*;

public class CalculateMoles {
	public static final double I = 131.17;
	public static final double L = 131.17;
	public static final double K = 146.19;
	public static final double M = 149.21;
	public static final double F = 165.19;
	public static final double T = 119.12;
	public static final double W = 204.23;
	public static final double V = 117.15;
	public static final double R = 174.20;
	public static final double H = 155.16;
	public static final double A = 89.09;
	public static final double N = 132.12;
	public static final double D = 133.10;
	public static final double C = 121.16;
	public static final double E = 147.13;
	public static final double Q = 146.15;
	public static final double G = 75.07;
	public static final double P = 115.13;
	public static final double S = 105.09;
	public static final double Y = 181.19;
	public static final double B = 244.31;
	public static final double Ac = 43.04;
	
	private List<Double> molecularWeight;
	private Map<String, Double> aminoAcidWeight;
	
	public CalculateMoles(ArrayList<Double> weight) {
		this.molecularWeight = new ArrayList<Double>();
		for (int i = 0; i < weight.size(); i++) {
			this.molecularWeight.add(weight.get(i));
		}
//		this.aminoAcidWeight();
//		for (int i = 0; i < sequence.length(); i++) {
//			String aa = sequence.substring(i, i + 1);
//			double weight = this.aminoAcidWeight.get(aa);
//			this.molecularWeight.add(i, weight);
//		}
	}
	
	public String toString() {
		return this.molecularWeight.toString();
	}
	
	public int size() {
		return this.molecularWeight.size();
	}
	
	private void aminoAcidWeight() {
		this.aminoAcidWeight = new TreeMap<String, Double>();
		String aminoAcids = "ILKMFTWVRHANDCEQGPSYB";
		ArrayList<Double> weight = new ArrayList<Double>();
		weight.add(I);
		weight.add(L);
		weight.add(K);
		weight.add(M);
		weight.add(F);
		weight.add(T);
		weight.add(W);
		weight.add(V);
		weight.add(R);
		weight.add(H);
		weight.add(A);
		weight.add(N);
		weight.add(D);
		weight.add(C);
		weight.add(E);
		weight.add(Q);
		weight.add(G);
		weight.add(P);
		weight.add(S);
		weight.add(Y);
		weight.add(B);
		for (int i = 0; i < aminoAcids.length(); i++) {
			String aa = aminoAcids.substring(i, i + 1);
			this.aminoAcidWeight.put(aa, weight.get(i));
		}
	}
	
	public ArrayList<Double> runningSum(int window) {
		if (window < 0 || window > this.molecularWeight.size()) {
			throw new IllegalArgumentException();
		}
		ArrayList<Double> sum = new ArrayList<Double>();
		for (int i = 0; i < -window + (this.molecularWeight.size() + 1); i++) {
			List<Double> windowSum = this.molecularWeight.subList(i, window + i);
			double tempSum = 0.0;
			for (int j = 0; j < windowSum.size(); j++) {
				tempSum += windowSum.get(j);
			}
			tempSum = tempSum - (18 * (windowSum.size() - 1));
			sum.add(Math.round(tempSum * 100.0) / 100.0);
		}
		return sum;
	}
	
	public ArrayList<Double> cyclicRunningSum(int window) {
		if (window < 0 || window > this.molecularWeight.size()) {
			throw new IllegalArgumentException();
		}
		ArrayList<Double> sum = new ArrayList<Double>();
		for (int i = 0; i < -window + (this.molecularWeight.size() + 1); i++) {
			List<Double> windowSum = this.molecularWeight.subList(i, window + i);
			double tempSum = 0.0;
			for (int j = 0; j < windowSum.size(); j++) {
				tempSum += windowSum.get(j);
			}
			tempSum = tempSum - (18 * (windowSum.size() - 1));
			
			sum.add(Math.round(tempSum * 100.0) / 100.0);
		}
		return sum;
	}
	
	public Map<String, Double> mapSequence(String sequence, ArrayList<Double> sum, int window) {
		Map<String, Double> fragments = new TreeMap<String, Double>();
		for (int i = 0; i < -window + (sequence.length() + 1); i++) {
			double fragmentSum = sum.get(i);
			String fragment = sequence.substring(i, window + i);
			fragments.put(fragment, fragmentSum);
		}
		return fragments;
	}
	
}
