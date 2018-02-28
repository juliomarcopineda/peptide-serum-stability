import java.io.*;
import java.util.*;

public class M2pepbiotin {
	public static final double BIOTIN = 244.31;
	public static final double a = 43.04;
	public static final int WINDOW = 8;
	public static final double THRESHOLD = 2.0;
	
	public static void main(String[] args) throws FileNotFoundException {
		printIntro();
		ArrayList<Double> weight = new ArrayList<Double>();
		Scanner console = new Scanner(new File("AcM2pep(WtoY)Biotin.txt"));
		while (console.hasNext()) {
			weight.add(console.nextDouble());
		}
		weight.add(BIOTIN);
		String sequence = "YEQDPWGVKWWYGGGSKKKB";
		CalculateMoles m2pep = new CalculateMoles(weight);
		peptideIntro(weight, sequence, m2pep);
		Map<String, Double> fragments = determineFragments(weight, sequence, m2pep);
		
		Scanner input = new Scanner(System.in);
		//compareData(fragments, input, sequence);
		
	}
	
	public static void printIntro() {
		System.out.println("Serum stability study program");
		System.out.println("Created by Julio Pineda");
		System.out.println();
		System.out.println("This program will produce all the possible fragments");
		System.out.println("of a given linear amino acid sequence with");
		System.out.println("a fragment size of " + WINDOW +".");
		System.out.println();
	}
	
	public static void peptideIntro(ArrayList<Double> weight, String sequence, 
			                        CalculateMoles m2pep) {
		System.out.println("M2pepBiotin sequence: " + sequence);
		System.out.println("Total Length: " + m2pep.size());
		System.out.println("Total molecular weight: " + m2pep.runningSum(m2pep.size()).get(0));
		System.out.println();
	}
	
	public static Map<String, Double> determineFragments(ArrayList<Double> weight, String sequence, 
										  CalculateMoles m2pep) {
		System.out.println("Fragment size: " + WINDOW);
		ArrayList<Double> sum = m2pep.runningSum(WINDOW);
		System.out.println("Number of fragments: " + sum.size());
		Map<String, Double> fragments = m2pep.mapSequence(sequence, sum, WINDOW);
		for (String s : fragments.keySet()) {
			double sumFragment = fragments.get(s);
			System.out.println(s + "-----" + sumFragment);
		}
		System.out.println();
		return fragments;
	}
	
	public static void compareData(Map<String, Double> fragments, Scanner input, String sequence) {
		System.out.print("Do you want to compare your data with generated fragments? (y/n) ");
		String answer = input.next();
		if (answer.equalsIgnoreCase("y")) {
			System.out.println("Compare with Mass Spec data with threshold of " + THRESHOLD);
			System.out.print("Molecular weight from data? ");
			int matches = 0;
			double data = input.nextDouble();
			for (String s : fragments.keySet()) {
				double theoreticalWeight = fragments.get(s);
				if (data > theoreticalWeight - THRESHOLD && data < theoreticalWeight + THRESHOLD) {
					matches++;
					System.out.println("You have found " + matches + " match(es)");
					System.out.println(s);
					System.out.println("Theoretical weight: " + theoreticalWeight);
					System.out.println("Data: " + data);
					int cleavage = sequence.indexOf(s) + 1;
					if (cleavage == 1) {	
						System.out.println("Possible cleavage site before " + cleavage + "st amino acid");
						System.out.println("and after " + (cleavage + s.length() - 1) + "th amino acid");
						System.out.println();
					} else if (cleavage == 2) {
						System.out.println("Possible cleavage site before " + cleavage + "2nd amino acid");
						System.out.println("and after " + (cleavage + s.length() - 1) + "th amino acid");
						System.out.println();
					} else {
						System.out.println("Possible cleavage site before " + cleavage + "th amino acid");
						System.out.println("and after " + (cleavage + s.length() - 1) + "th amino acid");
						System.out.println();
					}
				}
			}
			if (matches == 0) {
				System.out.println("Sorry no matches.");
			}
		} else {
			System.out.println("Have a good day!");
		}
	}
}
