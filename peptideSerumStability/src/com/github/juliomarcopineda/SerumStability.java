package com.github.juliomarcopineda;

import java.util.List;

import com.github.juliomarcopineda.peptide.Peptide;

public class SerumStability {
	public static void main(String[] args) {
		String inputFile = "../peptide-test.txt";
		
		InputParser parser = new InputParser(inputFile);
		parser.parse();
		List<Peptide> peptides = parser.getPeptides();
		
		printInputData(peptides);
	}
	
	private static void printInputData(List<Peptide> peptides) {
		peptides.stream()
			.forEach(peptide -> {
				String seq = peptide.getSequence();
				
				System.out.println(seq);
				System.out.println(peptide.getType());
				
				if (peptide.getConnections() != null && !peptide.getConnections()
					.isEmpty()) {
					System.out.println(peptide.getConnections());
				}
				
				if (peptide.getMassSpecData() != null) {
					System.out.println(peptide.getMassSpecData());
				}
				
				if (!peptide.getGraph()
					.isEmpty()) {
					peptide.getGraph()
						.entrySet()
						.forEach(e -> {
							int source = e.getKey();
							List<Integer> targets = e.getValue();
							
							if (source < seq.length()) {
								System.out.println(seq.charAt(source) + " -> " + printTargets(targets, seq));
							}
							else {
								System.out.println(source + " -> " + printTargets(targets, seq));
							}
							
						});
				}
				
				System.out.println();
			});
		
	}
	
	private static String printTargets(List<Integer> targets, String sequence) {
		
		StringBuilder sb = new StringBuilder();
		
		for (int target : targets) {
			if (sb.length() != 0) {
				sb.append(", ");
			}
			
			if (target < sequence.length()) {
				sb.append(sequence.charAt(target));
			}
			else {
				sb.append(target);
			}
		}
		
		return sb.toString();
	}
}
