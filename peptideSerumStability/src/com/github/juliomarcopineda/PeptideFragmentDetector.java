package com.github.juliomarcopineda;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * PeptideFragmentDetector accepts a peptide represented as a String of one-letter amino acids from 
 * N-terminus to C-terminus. Then, this class will determine all the possible fragments from the 
 * peptide sequence. This class is also able to accept mass spectrometry data to suggest possible
 * fragments that exist in one's serum.
 * 
 * @author Julio Pineda
 *
 */
public class PeptideFragmentDetector {
	private static final double THRESHOLD = 2.0;
	
	private String peptide;
	private Map<Character, Double> aminoAcidWeight;
	
	private Map<String, Double> fragments;
	private Map<String, Integer> positions;
	
	public PeptideFragmentDetector(String peptide) {
		this.peptide = peptide;
		this.aminoAcidWeight = createAminoAcideWeightMap();
	}
	
	/**
	 * Generates all the possible fragments from the input peptide. Stores the information in a map of <fragment, molecular weight>.
	 */
	public void determineAllFragments() {
		System.out.println("Finding all possible fragments from input peptide...");
		
		Map<String, Double> fragments = new HashMap<>();
		Map<String, Integer> positions = new HashMap<>();
		for (int size = 1; size <= peptide.length(); size++) {

			int i = 0;
			while ((i + size) <= this.peptide.length()) {
				int startIndex = i;
				int endIndex = i + size;
				
				String fragment = this.peptide.substring(startIndex, endIndex);
				double weight = calculateFragmentWeight(fragment);
				
				
				fragments.put(fragment, weight);
				positions.put(fragment, startIndex + 1);
				
				i++;
			}
			
			
		}
		this.fragments = fragments;
		this.positions = positions;
		
		System.out.println();
		System.out.println("Finished generating all peptide fragments.");
		System.out.println();
	}
	
	public void suggestFragments(double massSpecData) {
		System.out.println("Comparing Mass Spectrometry Data with generated fragments with threshold of " + THRESHOLD);
		
		List<String> suggestedFragments = new ArrayList<>();
		
		for (Map.Entry<String, Double> entry : this.fragments.entrySet()) {
			String fragment = entry.getKey();
			double theoreticalWeight = entry.getValue();
			
			if ((massSpecData >= theoreticalWeight - THRESHOLD) && (massSpecData <= theoreticalWeight + THRESHOLD)) {
				suggestedFragments.add(fragment);
			}
		}
		
		if (suggestedFragments.isEmpty()) {
			System.out.println("No fragments found for data input.");
			System.out.println();
		}
		else {
			System.out.println("Found " + suggestedFragments.size() + " fragments that matched for data input...");
			for (String fragment : suggestedFragments) {
				System.out.println("FRAGMENT: " + fragment);
				System.out.println("Theoretical weight: " + Math.round(this.fragments.get(fragment) * 100.0) / 100.0);
				System.out.println("Starting position: " + this.positions.get(fragment));
				System.out.println();
			}
		}
	}
	
	public void printFragments() {
		System.out.println("Printing all possible fragments");
		
		this.fragments.entrySet()
			.stream()
			.forEach(e -> {
				System.out.println(e.getKey());
			});
		
		System.out.println();
	}
	
	public void printFragments(int size) {
		System.out.println("Printing all possible fragments of size " + size);
		
		this.fragments.entrySet()
			.stream()
			.filter(e -> e.getKey().length() == size)
			.forEach(e -> {
				System.out.println(e.getKey() + "----" + Math.round(e.getValue() * 100.0) / 100.0);
			});
		
		System.out.println();
	}
	
	/**
	 * Given a peptide fragment represented as a String, returns the molecular weight associated with the fragment.
	 * 
	 * @param fragment
	 * @return
	 */
	private double calculateFragmentWeight(String fragment) {
		double molecularWeight = 0;
		
		// Iterate over each character of the fragment and determine its individual molecular weight
		for (int i = 0; i < fragment.length(); i++) {
			char aminoAcid = fragment.charAt(i);
			molecularWeight += this.aminoAcidWeight.get(aminoAcid);
		}
		
		// Subtract the molecular weight of water
		molecularWeight = molecularWeight - (18.0 * (fragment.length() - 1));
		
		return molecularWeight;
	}
	
	private Map<Character, Double> createAminoAcideWeightMap() {
		Map<Character, Double> aminoAcidWeight = new HashMap<>();
		aminoAcidWeight.put('I', 131.17);
		aminoAcidWeight.put('L', 131.17);
		aminoAcidWeight.put('K', 146.19);
		aminoAcidWeight.put('M', 149.21);
		aminoAcidWeight.put('F', 165.19);
		aminoAcidWeight.put('T', 119.12);
		aminoAcidWeight.put('W', 204.23);
		aminoAcidWeight.put('V', 117.15);
		aminoAcidWeight.put('R', 174.20);
		aminoAcidWeight.put('H', 115.16);
		aminoAcidWeight.put('A', 89.09);
		aminoAcidWeight.put('N', 132.12);
		aminoAcidWeight.put('D', 133.10);
		aminoAcidWeight.put('C', 121.16);
		aminoAcidWeight.put('E', 147.13);
		aminoAcidWeight.put('Q', 146.15);
		aminoAcidWeight.put('G', 75.07);
		aminoAcidWeight.put('P', 115.13);
		aminoAcidWeight.put('S', 105.09);
		aminoAcidWeight.put('Y', 181.19);
		aminoAcidWeight.put('B', 224.31);
		aminoAcidWeight.put('1', 43.04); // Acetylated N-terminus
		
		return aminoAcidWeight;
	}
}
