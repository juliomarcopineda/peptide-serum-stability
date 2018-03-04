package com.github.juliomarcopineda;

import java.util.AbstractMap;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import com.github.juliomarcopineda.peptide.Peptide;
import com.github.juliomarcopineda.peptide.PeptideType;

public class FragmentAnalyzer {
	private Peptide peptide;
	private Map<Character, Double> weights;
	
	public FragmentAnalyzer(Peptide peptide) {
		this.peptide = peptide;
		this.weights = createAminoAcideWeightMap();
	}
	
	public FragmentAnalyzer() {
		this.weights = createAminoAcideWeightMap();
	}
	
	public void findFragment(Map<Integer, List<Integer>> graph, double data) {
		for (Map.Entry<Integer, List<Integer>> entry : graph.entrySet()) {
			int start = entry.getKey();
			
		}
	}
	
	public void walkGraph(int start, Map<Integer, List<Integer>> graph, double data) {
		List<Integer> fragmentIndex = new ArrayList<>();
		fragmentIndex.add(start);
		
		walkGraph(start, start, graph, fragmentIndex, data);
	}
	
	private void walkGraph(int root, int start, Map<Integer, List<Integer>> graph, List<Integer> fragmentIndex, double data) {
		double fragmentWeight = fragmentWeight(fragmentIndex);
		
		//		if ((fragmentWeight - data) > 10.0) { // fragment is too big
		//			System.out.println("Too big");
		//			return;
		//		}
		if (!graph.containsKey(start)) { // traversed the end of the peptide
			System.out.println("End of peptide");
			System.out.println();
			return;
		}
		//		else if (Math.abs(data - fragmentWeight) <= 2.0) { // found a fragment!
		//			System.out.println("FOUND ONE!");
		//			return;
		//		}
		else {
			List<Integer> targets = graph.get(start);
			
			//			if (targets.contains(before)) {
			//				targets.remove(targets.indexOf(before));
			//			}
			
			for (int i = 0; i < targets.size(); i++) {
				int target = targets.get(i);
				
				if (root == target) {
					System.out.println("Formed cycle: " + fragmentIndex.toString() + " " + target);
					System.out.println();
					return;
				}
				
				fragmentIndex.add(target);
				System.out.println(fragmentIndex);
				
				walkGraph(root, target, graph, fragmentIndex, data);
				fragmentIndex.remove(fragmentIndex.size() - 1);
				System.out.println(fragmentIndex);
			}
		}
	}
	
	private double fragmentWeight(List<Integer> fragmentIndex) {
		double sum = 0;
		
		for (int index : fragmentIndex) {
			
			char symbol = this.peptide.getSequence()
				.charAt(index);
			sum += this.weights.get(symbol);
		}
		
		sum = sum - (18.0 * (fragmentIndex.size() - 1));
		
		return sum;
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
	
	public static void main(String[] args) {
		String peptideSequence = "DGYEQDPWGVRYWYGKKKKKB";
		List<Integer> connections = Arrays.asList(0, 15);
		PeptideType type = PeptideType.AMIDE;
		Map<Integer, List<Integer>> graph = createGraphStructure(peptideSequence, connections, type);
		
		Peptide peptide = new Peptide();
		peptide.setSequence(peptideSequence);
		peptide.setConnections(connections);
		peptide.setType(type);
		peptide.setGraph(graph);
		
		FragmentAnalyzer analyzer = new FragmentAnalyzer(peptide);
		
		graph.entrySet()
			.forEach(e -> {
				int source = e.getKey();
				List<Integer> targets = e.getValue();
				
				System.out.println(source + " -> " + targets);
			});
		
		System.out.println();
		
		//		for (Map.Entry<Integer, List<Integer>> entry : graph.entrySet()) {
		//			int start = entry.getKey();
		//			analyzer.walkGraph(start, graph, 886.54);
		//			System.out.println();
		//			
		//			break;
		//			
		//		}
		
		int start = 1;
		analyzer.walkGraph(start, graph, 886.54);
		System.out.println();
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
