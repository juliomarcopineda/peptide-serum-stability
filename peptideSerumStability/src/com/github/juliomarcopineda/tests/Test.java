package com.github.juliomarcopineda.tests;

import java.util.Arrays;
import java.util.List;
import java.util.Map;

import com.github.juliomarcopineda.FragmentAnalyzer;
import com.github.juliomarcopineda.PeptideSerumStability;
import com.github.juliomarcopineda.peptide.Peptide;
import com.github.juliomarcopineda.peptide.PeptideType;

public class Test {
	public static void main(String[] args) {
		// Prevent equal connections
		// Prevent connections in decreasing order
		
		String peptideSequence = "CGYEQDPWGVRYWYGCKKKKB";
		
		for (int conn1 = 0; conn1 < peptideSequence.length() - 1; conn1++) {
			for (int conn2 = conn1 + 1; conn2 < peptideSequence.length(); conn2++) {
				if (conn1 == conn2) {
					continue;
				}
				
				List<Integer> connections = Arrays.asList(conn1, conn2);
				
				for (PeptideType type : PeptideType.values()) {
					
					try {
						Map<Integer, List<Integer>> graph = PeptideSerumStability.createGraphStructure(peptideSequence, connections, type);
						
						Peptide peptide = new Peptide();
						peptide.setSequence(peptideSequence);
						peptide.setConnections(connections);
						peptide.setType(type);
						peptide.setGraph(graph);
						
						FragmentAnalyzer analyzer = new FragmentAnalyzer(peptide);
						analyzer.findAllFragments()
							.measureAllFragmentWeights();
					}
					catch (Exception e) {
						System.out.println("ERROR!");
						System.out.println(type);
						System.out.println("Connection 1: " + conn1);
						System.out.println("Connection 2: " + conn2);
						
						e.printStackTrace();
						System.exit(-1);
						
					}
					
				}
			}
		}
		
		//		Map<String, Double> fragments = analyzer.getFragmentWeights();
		//				fragments.entrySet()
		//					.stream()
		//					.filter(e -> e.getKey()
		//						.contains("#"))
		//					.forEach(System.out::println);
		
		//		double threshold = 5.0;
		//		double massSpecData = 655.74;
		//		
		//		Map<String, Double> suggestedFragments = analyzer.suggestFragments(massSpecData, threshold);
		//		suggestedFragments.entrySet()
		//			.stream()
		//			.forEach(System.out::println);
	}
}
