package com.github.juliomarcopineda.experiments;

import java.util.AbstractMap;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import com.github.juliomarcopineda.FragmentAnalyzer;
import com.github.juliomarcopineda.peptide.Peptide;
import com.github.juliomarcopineda.peptide.PeptideType;

public class Test {
	public static void main(String[] args) {
		String peptideSequence = "CGYEQDPWGVRYWYGCKKKKB";
		List<Integer> connections = Arrays.asList(0, 15);
		// List<Integer> connections = new ArrayList<>();
		PeptideType type = PeptideType.DFBP;
		Map<Integer, List<Integer>> graph = createGraphStructure(peptideSequence, connections, type);
		
		Peptide peptide = new Peptide();
		peptide.setSequence(peptideSequence);
		peptide.setConnections(connections);
		peptide.setType(type);
		peptide.setGraph(graph);
		
		FragmentAnalyzer analyzer = new FragmentAnalyzer(peptide);
		analyzer.findAllFragments()
			.measureAllFragmentWeights();
		
		Map<String, Double> fragments = analyzer.getFragmentWeights();
		fragments.entrySet()
			.stream()
			.filter(e -> e.getKey()
				.charAt(0) == '#')
			.forEach(System.out::println);
		
		//		double threshold = 5.0;
		//		double massSpecData = 655.74;
		//		
		//		Map<String, Double> suggestedFragments = analyzer.suggestFragments(massSpecData, threshold);
		//		suggestedFragments.entrySet()
		//			.stream()
		//			.forEach(System.out::println);
	}
	
	private static Map<Integer, List<Integer>> createGraphStructure(String peptideSequence, List<Integer> connections, PeptideType type) {
		
		// Create the graph structure of the peptide base sequence
		Map<Integer, List<Integer>> graph = IntStream.range(0, peptideSequence.length() - 1)
			.mapToObj(source -> {
				int target = source + 1;
				
				List<Integer> targets = new ArrayList<>();
				targets.add(target);
				
				return new AbstractMap.SimpleEntry<Integer, List<Integer>>(source, targets);
			})
			.sorted(Map.Entry.comparingByKey())
			.collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue, (collision1, collision2) -> collision1, LinkedHashMap::new));
		
		// Add any cyclic connections if connections is not empty
		if (!connections.isEmpty()) {
			switch (type) {
				case AMIDE:
					for (int i = 0; i < connections.size(); i++) {
						if (i % 2 == 0) {
							graph.get(connections.get(i))
								.add(connections.get(i + 1));
						}
						else {
							graph.get(connections.get(i))
								.add(connections.get(i - 1));
						}
					}
					
					break;
				
				case DFBP:
					int dfbpIndex = peptideSequence.length();
					
					for (int connection : connections) {
						// Add connections from DFBP
						if (!graph.containsKey(dfbpIndex)) {
							List<Integer> targets = new ArrayList<>();
							targets.add(connection);
							
							graph.put(dfbpIndex, targets);
						}
						else {
							graph.get(dfbpIndex)
								.add(connection);
						}
						
						// Add connections to DFBP
						graph.get(connection)
							.add(dfbpIndex);
					}
					
					break;
				case DISULFIDE:
					// Create disulfide bridge
					int s1Index = peptideSequence.length();
					int s2Index = s1Index + 1;
					
					List<Integer> s1ToS2 = new ArrayList<>();
					s1ToS2.add(s2Index);
					
					List<Integer> s2ToS1 = new ArrayList<>();
					s2ToS1.add(s1Index);
					
					graph.put(s1Index, s1ToS2);
					graph.put(s2Index, s2ToS1);
					
					// Add connections from peptide base
					graph.get(connections.get(0))
						.add(s1Index);
					graph.get(connections.get(1))
						.add(s2Index);
					
					// Add connections from disulfide bridge
					graph.get(s1Index)
						.add(connections.get(0));
					graph.get(s2Index)
						.add(connections.get(1));
					
					break;
				case LINEAR:
					break;
			}
		}
		
		return graph;
	}
}
