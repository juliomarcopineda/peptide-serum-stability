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
	
	private List<List<Integer>> fragments;
	
	public FragmentAnalyzer(Peptide peptide) {
		this.peptide = peptide;
		this.weights = createAminoAcideWeightMap();
		this.fragments = new ArrayList<>();
	}
	
	public List<List<Integer>> getFragments() {
		return fragments;
	}
	
	public void setFragments(List<List<Integer>> fragments) {
		this.fragments = fragments;
	}
	
	/**
	 * Initiates the process of finding all the possible fragments of the peptide.
	 */
	public void findAllFragments() {
		Map<Integer, List<Integer>> graph = this.peptide.getGraph();
		for (Map.Entry<Integer, List<Integer>> entry : graph.entrySet()) {
			int start = entry.getKey();
			
			walkGraph(start, graph);
		}
	}
	
	/**
	 * Given a starting node (amino acid) and the graph representation of the peptide, traverse all the possible
	 * paths from that starting node.
	 * 
	 * Recursive backtracking is the algorithm used to traverse every possible path from the starting node. See the private helper method to check the
	 * details of the recursive backtracking implementation.
	 * 
	 * @param start
	 * @param graph
	 */
	public void walkGraph(int start, Map<Integer, List<Integer>> graph) {
		List<Integer> fragmentIndex = new ArrayList<>();
		fragmentIndex.add(start);
		
		walkGraph(start, -1, start, graph, fragmentIndex);
	}
	
	/**
	 * Helper method for the recursive backtracking implementation of traversing the graph given a starting node.
	 * 
	 * The root parameter is the starting node of the algorithm. If the root is visited, a cycle is formed and the traversing is terminated.
	 * The before parameter keeps track of the last node visited by the algorithm. This prevents the traversing to go back prematurely and also prevent
	 * an infinite loop in the linker of cyclic peptides.
	 * The start parameter indicates the current location of traversing the graph.
	 * The fragmentIndex keeps track of the traversed nodes during the recursive backtrackign step.
	 * 
	 * This method also saves all possible fragments visited by the algorithm into the fragments field.
	 * 
	 * @param root
	 * @param before
	 * @param start
	 * @param graph
	 * @param fragmentIndex
	 */
	private void walkGraph(int root, int before, int start, Map<Integer, List<Integer>> graph, List<Integer> fragmentIndex) {
		if (!graph.containsKey(start)) { // Traversed the end of the peptide
			return;
		}
		else {
			List<Integer> targets = graph.get(start);
			
			// Iterate over the different choices during a branching path
			for (int i = 0; i < targets.size(); i++) {
				int target = targets.get(i);
				
				if (target == before) { // Prevents walk to backtrack
					continue;
				}
				
				if (root == target) { // Formed a cycle and terminates walk in this direction
					return;
				}
				
				// Keeping track of traversed nodes
				fragmentIndex.add(target);
				
				// Create a new list for each valid fragment
				List<Integer> copyFragmentIndex = new ArrayList<>();
				copyFragmentIndex.addAll(fragmentIndex);
				this.fragments.add(copyFragmentIndex);
				
				// Recursive step
				walkGraph(root, start, target, graph, fragmentIndex);
				
				// Backtracking step when finding a terminal/base case
				fragmentIndex.remove(fragmentIndex.size() - 1);
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
		String peptideSequence = "CGYEQDPWGVRYWYGCKKKKB";
		List<Integer> connections = Arrays.asList(0, 15);
		PeptideType type = PeptideType.DFBP;
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
		
		analyzer.findAllFragments();
		List<List<Integer>> fragments = analyzer.getFragments();
		fragments.stream()
			.forEach(System.out::println);
		System.out.println();
		System.out.println(fragments.size());
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
