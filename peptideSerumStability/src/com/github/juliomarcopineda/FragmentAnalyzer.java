package com.github.juliomarcopineda;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.github.juliomarcopineda.peptide.Peptide;
import com.github.juliomarcopineda.peptide.PeptideType;

public class FragmentAnalyzer {
	private Peptide peptide;
	private Map<Character, Double> weights;
	
	private List<List<Integer>> fragments;
	private Map<String, Double> fragmentWeights;
	
	public FragmentAnalyzer(Peptide peptide) {
		this.peptide = peptide;
		this.weights = createAminoAcideWeightMap();
		this.fragments = new ArrayList<>();
	}
	
	public Map<String, Double> getFragmentWeights() {
		return fragmentWeights;
	}
	
	public void setFragmentWeights(Map<String, Double> fragmentWeights) {
		this.fragmentWeights = fragmentWeights;
	}
	
	public List<List<Integer>> getFragments() {
		return fragments;
	}
	
	public void setFragments(List<List<Integer>> fragments) {
		this.fragments = fragments;
	}
	
	/**
	 * Populates the map of fragment (represented by String) and its theoretical molecular weight.
	 * 
	 * @return
	 */
	public FragmentAnalyzer measureAllFragmentWeights() {
		Map<String, Double> fragmentWeights = new HashMap<>();
		
		PeptideType type = this.peptide.getType();
		
		for (List<Integer> fragmentIndex : this.fragments) {
			String fragment = getPeptideStringRepresentation(fragmentIndex, type);
			double weight = calculateFragmentWeight(fragment);
			
			fragmentWeights.put(fragment, weight);
		}
		
		this.fragmentWeights = fragmentWeights;
		
		return this;
	}
	
	/**
	 * Given the peptide represented by its index and the peptide type, return the String representation 
	 * of the peptide.
	 * 
	 * @param peptideIndex
	 * @param type
	 * @return
	 */
	private String getPeptideStringRepresentation(List<Integer> peptideIndex, PeptideType type) {
		String peptideSequence = this.peptide.getSequence();
		
		StringBuilder stringBuilder = new StringBuilder();
		
		for (int index : peptideIndex) {
			if (index < peptideSequence.length()) {
				stringBuilder.append(peptideSequence.charAt(index));
			}
			else {
				switch (type) {
					case DFBP:
						stringBuilder.append("2");
						break;
					case DISULFIDE:
						stringBuilder.append("S");
						break;
					case LINEAR:
						break;
					case AMIDE:
						break;
					default:
						break;
				}
			}
		}
		
		return stringBuilder.toString();
	}
	
	/**
	 * Given a fragment with its String representation, calculates its theoretical molecular weight.
	 * 
	 * @param fragment
	 * @return
	 */
	private double calculateFragmentWeight(String fragment) {
		double sum = 0;
		
		for (int i = 0; i < fragment.length(); i++) {
			char symbol = fragment.charAt(i);
			sum += this.weights.get(symbol);
		}
		
		sum = sum - (18.0 * (fragment.length() - 1));
		
		return sum;
	}
	
	/**
	 * Initiates the process of finding all the possible fragments of the peptide.
	 */
	public FragmentAnalyzer findAllFragments() {
		Map<Integer, List<Integer>> graph = this.peptide.getGraph();
		for (Map.Entry<Integer, List<Integer>> entry : graph.entrySet()) {
			int start = entry.getKey();
			
			walkGraph(start, graph);
		}
		
		return this;
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
	
	/**
	 * Populate the amino acid weight map using the weights.csv file.
	 * 
	 * @return
	 */
	private Map<Character, Double> createAminoAcideWeightMap() {
		Map<Character, Double> aminoAcidWeight = new HashMap<>();
		
		try {
			BufferedReader reader = new BufferedReader(new FileReader("../weights.csv"));
			String line = reader.readLine(); // ignore header
			while ((line = reader.readLine()) != null) {
				
				String[] split = line.split(",");
				char symbol = split[0].charAt(0);
				double weight = Double.parseDouble(split[1]);
				
				aminoAcidWeight.put(symbol, weight);
			}
			
			reader.close();
		}
		catch (IOException e) {
			e.printStackTrace();
		}
		
		return aminoAcidWeight;
	}
	
}
