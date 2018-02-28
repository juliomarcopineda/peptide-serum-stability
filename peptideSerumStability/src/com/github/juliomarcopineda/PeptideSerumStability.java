package com.github.juliomarcopineda;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;

/**
 * PeptideSerumStability is a command-line tool that generates peptide fragments from input data and detects possible peptide fragments 
 * from Mass Spectrometry Data.
 * 
 * @author Julio Pineda
 *
 */
public class PeptideSerumStability {

	public static void main(String[] args) throws IOException {
		String peptide = args[0].toUpperCase();
		
		if (peptide.isEmpty()) {
			System.out.println("Please input a peptide sequence.");
			System.exit(0);
		}
		else if (peptide.equals("TEST")) {
			// M2pep
			peptide = "YEQDPWGVKWWYGGGSKKKB";
		}
		
		printIntro(peptide);
		
		PeptideFragmentDetector detector = new PeptideFragmentDetector(peptide);
		detector.determineAllFragments();
		
		BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
		String prompt = "in"; // Go in the while loop for the first time
		
		while (!prompt.substring(0, 1).toUpperCase().equals("Q")) {
			menuPrompt();
			prompt = br.readLine();
			
			if (prompt.substring(0, 1).toUpperCase().equals("C")) {
				String proceed = "Y";
				
				while (proceed.substring(0, 1).toUpperCase().equals("Y")) {
					System.out.print("Enter data from mass spectrometry (must be decimal): ");
					double data = Double.parseDouble(br.readLine());
					System.out.println();
					detector.suggestFragments(data);
					
					System.out.print("More data? (Y/N) ");
					proceed = br.readLine();
				}
			}
			else if (prompt.substring(0, 1).toUpperCase().equals("P")) {
				String proceed = "Y";

				while (proceed.substring(0, 1).toUpperCase().equals("Y")) {
					System.out.print("What fragment size to print? (Must be integer) ");
					int size = Integer.parseInt(br.readLine());
					detector.printFragments(size);
					
					System.out.print("Print more? (Y/N) ");
					proceed = br.readLine();
				}
			}
			
			System.out.println();
		}
		
		System.out.println("Goodbye!");
		br.close();
	}

	private static void printIntro(String peptide) {
		System.out.println("-------------------------");
		System.out.println("Peptide Serum Stability");
		System.out.println("-------------------------");
		System.out.println();
		System.out.println("This program will generate all possible fragments of the given input peptide sequence");
		System.out.println("and will calculate their corresponding theoretical molecular weights.");
		System.out.println();
		System.out.println("Peptide: " + peptide);
		System.out.println();
	}
	
	private static void menuPrompt() {
		System.out.println("What do you want to do?");
		System.out.println("| (C)ompare | (P)rint | (Q)uit |");
		System.out.print("Enter: ");
	}
}
