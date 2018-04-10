package com.github.juliomarcopineda;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.github.juliomarcopineda.peptide.Peptide;
import com.github.juliomarcopineda.peptide.PeptideType;

/**
 * FragmentAnalyzer accepts a Peptide and determines all the possible fragments that can result from a peptide serum stability study. 
 * 
 * @author Julio Pineda
 *
 */
public class FragmentAnalyzer {
	private Peptide peptide;
	private Map<Character, Double> weights;
	
	private List<List<Integer>> fragments;
	private Map<String, Double> fragmentWeights;
	
	/**
	 * Constructor that accepts a Peptide object.
	 * 
	 * @param peptide
	 */
	public FragmentAnalyzer(Peptide peptide) {
		this.peptide = peptide;
		this.weights = createAminoAcideWeightMap();
		this.fragments = new ArrayList<>();
	}
	
	public FragmentAnalyzer() {
		this.weights = createAminoAcideWeightMap();
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
	 * Given the mass spec data and a threshold, returns all the fragments that are withing this threshold.
	 * 
	 * @param massSpecData
	 * @param threshold
	 * @return
	 */
	public Map<String, Double> suggestFragments(double massSpecData, double threshold) {
		Map<String, Double> suggestedFragments = new HashMap<>();
		
		for (Map.Entry<String, Double> entry : this.fragmentWeights.entrySet()) {
			String fragment = entry.getKey();
			double theoreticalWeight = entry.getValue();
			
			double diff = Math.abs(massSpecData - theoreticalWeight);
			if (diff <= threshold) {
				suggestedFragments.put(fragment, theoreticalWeight);
			}
		}
		
		return suggestedFragments;
	}
	
	/**
	 * Populates the map of fragment (represented by String) and its theoretical molecular weight.
	 * 
	 * @return
	 */
	public FragmentAnalyzer measureAllFragmentWeights() {
		System.out.println("Calculating molecular weights of all fragments...");
		
		Map<String, Double> fragmentWeights = new HashMap<>();
		PeptideType type = this.peptide.getType();
		
		for (List<Integer> fragmentIndex : this.fragments) {
			String fragment = getPeptideStringRepresentation(fragmentIndex, type);
			double weight = calculateFragmentWeight(fragment);
			
			fragmentWeights.put(fragment, weight);
		}
		
		if (!type.equals(PeptideType.LINEAR)) {
			findBranchedAndCyclicFragments(fragmentWeights);
		}
		
		this.fragmentWeights = fragmentWeights;
		
		System.out.println("Done!");
		System.out.println();
		
		return this;
	}
	
	/**
	 * Determines all the posisble branched fragments of a Peptide and populated the map of <fragment, weight>. 
	 * Note that the current implementation can only handle a cyclic peptide with only one cycle.
	 * 
	 * @param fragmentWeights
	 */
	private void findBranchedAndCyclicFragments(Map<String, Double> fragmentWeights) {
		PeptideType type = this.peptide.getType();
		List<Integer> connections = this.peptide.getConnections();
		String peptideSequence = this.peptide.getSequence();
		
		Map<Integer, List<List<Integer>>> connectionInFragments = findConnectionsInFragments(connections);
		
		int connection1 = connections.get(0);
		int connection2 = connections.get(1);
		
		List<List<Integer>> fragmentsWith1 = connectionInFragments.get(connection1);
		List<List<Integer>> fragmentsWith2 = connectionInFragments.get(connection2);
		
		// Stop finding branched fragments if there are no fragments with connections
		if (fragmentsWith1 == null || fragmentsWith2 == null || fragmentsWith1.isEmpty() || fragmentsWith2.isEmpty()) {
			return;
		}
		
		// Start building branched fragment using linear fragment with the first connection
		for (List<Integer> fragmentWith1 : fragmentsWith1) {
			List<Integer> afterConnection1 = fragmentWith1.subList(fragmentWith1.indexOf(connection1), fragmentWith1.size());
			
			StringBuilder sb1 = new StringBuilder();
			sb1.append(getPeptideStringRepresentation(fragmentWith1, type));
			
			// Append the linkers to the fragment (DFBP or S or SS) to form branched fragments
			switch (type) {
				case DFBP:
					// Skip if the fragment contains DFBP
					if (fragmentWith1.contains(peptideSequence.length())) {
						continue;
					}
					
					// Only append DFBP if the connection is not in the beginning
					if (fragmentWith1.get(0) != connection1) {
						sb1.append("#2");
						fragmentWeights.put(sb1.toString(), calculateBranchedFragmentWeight(sb1.toString()));
					}
					
					break;
				case DISULFIDE:
					// Skip if the fragment contains S or SS
					if (fragmentWith1.contains(peptideSequence.length()) || fragmentWith1.contains(peptideSequence.length() + 1)) {
						continue;
					}
					
					// Only append S then SS if the connection is not in the beginning
					if (fragmentWith1.get(0) != connection1) {
						sb1.append("#S");
						fragmentWeights.put(sb1.toString(), calculateBranchedFragmentWeight(sb1.toString()));
						sb1.append("S");
						fragmentWeights.put(sb1.toString(), calculateBranchedFragmentWeight(sb1.toString()));
					}
					
					break;
				case AMIDE:
					break;
				case LINEAR:
					break;
				default:
					break;
			}
			
			// Start building branched fragments using the linear fragment with the second connection
			for (List<Integer> fragmentWith2 : fragmentsWith2) {
				
				// Skip if appending two fragments would form a linear peptide
				if (isLinear(fragmentWith1, connection1, fragmentWith2, connection2)) {
					continue;
				}
				
				List<Integer> beforeConnection2 = fragmentWith2.subList(0, fragmentWith2.indexOf(connection2) + 1);
				
				StringBuilder sb2 = new StringBuilder();
				sb2.append(getPeptideStringRepresentation(fragmentWith1, type));
				if (isValidBranchedFragment(afterConnection1, beforeConnection2)) {
					switch (type) {
						case DFBP:
							// Skip if the fragment contains DFBP
							if (fragmentWith2.contains(peptideSequence.length())) {
								continue;
							}
							
							sb2.append("#2#" + getPeptideStringRepresentation(fragmentWith2, type));
							fragmentWeights.put(sb2.toString(), calculateBranchedFragmentWeight(sb2.toString()));
							break;
						case DISULFIDE:
							// Skip if the fragment contains S or SS
							if (fragmentWith2.contains(peptideSequence.length()) || fragmentWith2.contains(peptideSequence.length() + 1)) {
								continue;
							}
							
							sb2.append("#SS#" + getPeptideStringRepresentation(fragmentWith2, type));
							fragmentWeights.put(sb2.toString(), calculateBranchedFragmentWeight(sb2.toString()));
							break;
						case AMIDE:
							sb2.append("#" + getPeptideStringRepresentation(fragmentWith2, type));
							fragmentWeights.put(sb2.toString(), calculateBranchedFragmentWeight(sb2.toString()));
							break;
						default:
							break;
					}
				}
			}
		}
		
		// Start forming branched fragments starting with fragments with the second connection
		for (List<Integer> fragmentWith2 : fragmentsWith2) {
			List<Integer> beforeConnection2 = fragmentWith2.subList(0, fragmentWith2.indexOf(connection2) + 1);
			
			StringBuilder sb1 = new StringBuilder();
			sb1.append(getPeptideStringRepresentation(fragmentWith2, type));
			
			// Append the linkers (DFBP or S or SS) to form possible branched fragments
			switch (type) {
				case DFBP:
					// Skip if linear fragment contains DFBP
					if (fragmentWith2.contains(peptideSequence.length())) {
						continue;
					}
					
					// Only append DFBP if connection is not in the end of the fragment
					if (fragmentWith2.get(fragmentWith2.size() - 1) != connection2) {
						sb1.append("#2");
						fragmentWeights.put(sb1.toString(), calculateBranchedFragmentWeight(sb1.toString()));
					}
					
					break;
				case DISULFIDE:
					// Skip if linear fragment contains S or SS
					if (fragmentWith2.contains(peptideSequence.length()) || fragmentWith2.contains(peptideSequence.length() + 1)) {
						continue;
					}
					
					// Only append S then SS if connection is not in the end of the fragment
					if (fragmentWith2.get(fragmentWith2.size() - 1) != connection2) {
						sb1.append("#S");
						fragmentWeights.put(sb1.toString(), calculateBranchedFragmentWeight(sb1.toString()));
						sb1.append("S");
						fragmentWeights.put(sb1.toString(), calculateBranchedFragmentWeight(sb1.toString()));
					}
					
					break;
				default:
					break;
			}
			
			// Start building branched fragments using linear fragments with the first connection
			for (List<Integer> fragmentWith1 : fragmentsWith1) {
				
				// Skip if appending two fragments would form a linear peptide
				if (isLinear(fragmentWith1, connection1, fragmentWith2, connection2)) {
					continue;
				}
				
				List<Integer> afterConnection1 = fragmentWith1.subList(fragmentWith1.indexOf(connection1), fragmentWith1.size());
				
				StringBuilder sb2 = new StringBuilder();
				sb2.append(getPeptideStringRepresentation(fragmentWith1, type));
				if (isValidBranchedFragment(afterConnection1, beforeConnection2)) {
					switch (type) {
						case DFBP:
							// Skip if fragment contains DFBP
							if (fragmentWith1.contains(peptideSequence.length())) {
								continue;
							}
							
							sb2.append("#2#" + getPeptideStringRepresentation(fragmentWith1, type));
							
							fragmentWeights.put(sb2.toString(), calculateBranchedFragmentWeight(sb2.toString()));
							break;
						case DISULFIDE:
							// Skip if fragment contains S or SS
							if (fragmentWith1.contains(peptideSequence.length()) || fragmentWith1.contains(peptideSequence.length() + 1)) {
								continue;
							}
							
							sb2.append("#SS#" + getPeptideStringRepresentation(fragmentWith1, type));
							fragmentWeights.put(sb2.toString(), calculateBranchedFragmentWeight(sb2.toString()));
							break;
						case AMIDE:
							sb2.append("#" + getPeptideStringRepresentation(fragmentWith1, type));
							fragmentWeights.put(sb2.toString(), calculateBranchedFragmentWeight(sb2.toString()));
							break;
						case LINEAR:
							break;
						default:
							break;
					}
				}
			}
		}
		
		List<List<Integer>> possibleCyclicFragments = connectionInFragments.get(-1);
		
		// Stop finding cyclic fragments if there are no possible cyclic fragments
		if (possibleCyclicFragments == null || possibleCyclicFragments.isEmpty()) {
			return;
		}
		
		for (List<Integer> possibleCyclicFragmentIndex : possibleCyclicFragments) {
			if (isCyclicFragment(possibleCyclicFragmentIndex, connections)) {
				StringBuilder sb = new StringBuilder();
				sb.append("#" + getPeptideStringRepresentation(possibleCyclicFragmentIndex, type) + "#");
				fragmentWeights.put(sb.toString(), calculateCyclicFragmentWeight(sb.toString(), type));
			}
		}
		
	}
	
	/**
	 * Determines if appending fragment1 and fragment2 results in a linear peptide.
	 * 
	 * @param fragment1
	 * @param connection1
	 * @param fragment2
	 * @param connection2
	 * @return
	 */
	private boolean isLinear(List<Integer> fragment1, int connection1, List<Integer> fragment2, int connection2) {
		if (fragment1.get(0) == connection1 && fragment2.get(fragment2.size() - 1) == connection2) {
			return true;
		}
		else if (fragment1.get(fragment1.size() - 1) == connection1 && fragment2.get(0) == connection2) {
			return true;
		}
		else {
			return false;
		}
	}
	
	/**
	 * Given the amino acids after the first connection and the amino acids before a connection, determines if the branched fragments is truly valid.
	 * A branched fragment is only valid if the two inputs do not share any elements.
	 * 
	 * Returns true if the branched fragments is valid.
	 * 
	 * @param afterConnection1
	 * @param beforeConnection2
	 * @return
	 */
	private boolean isValidBranchedFragment(List<Integer> afterConnection1, List<Integer> beforeConnection2) {
		int afterConnection1Size = afterConnection1.size();
		int beforeConnection2Size = beforeConnection2.size();
		
		if (beforeConnection2Size < afterConnection1Size) {
			List<Integer> difference = new ArrayList<>(afterConnection1);
			difference.removeAll(beforeConnection2);
			
			if (difference.size() < afterConnection1Size) {
				return false;
			}
		}
		else {
			List<Integer> difference = new ArrayList<>(beforeConnection2);
			difference.removeAll(afterConnection1);
			
			if (difference.size() < beforeConnection2Size) {
				return false;
			}
		}
		
		return true;
	}
	
	/**
	 * Determines if a possible cyclic fragment is truly a cyclic fragment.
	 * 
	 * @param possibleCyclicFragmentIndex
	 * @param connections
	 * @return
	 */
	private boolean isCyclicFragment(List<Integer> possibleCyclicFragmentIndex, List<Integer> connections) {
		
		for (int connection : connections) {
			int index = possibleCyclicFragmentIndex.indexOf(connection);
			
			if (index != connection) {
				return false;
			}
		}
		
		String peptideSequence = peptide.getSequence();
		if (possibleCyclicFragmentIndex.contains(peptideSequence.length()) || possibleCyclicFragmentIndex.contains(peptideSequence.length() + 1)) {
			return false;
		}
		
		return true;
	}
	
	/**
	 * Given the String representation of a branched fragment, calculates the molecular weight of this branched fragment.
	 * 
	 * @param branchedFragment
	 * @return
	 */
	private double calculateBranchedFragmentWeight(String branchedFragment) {
		double sum = 0;
		
		String[] split = branchedFragment.split("#");
		for (String fragment : split) {
			double weight = calculateFragmentWeight(fragment);
			sum += weight;
		}
		
		if (!this.peptide.getType()
			.equals(PeptideType.AMIDE)) {
			sum = sum - (18.0 * (split.length - 1));
		}
		
		return sum;
	}
	
	/**
	 * Caluculates the molecular weights of cyclic fragments given its string representation and the peptide type.
	 * 
	 * @param cyclicFragment
	 * @param type
	 * @return
	 */
	private double calculateCyclicFragmentWeight(String cyclicFragment, PeptideType type) {
		double sum = 0;
		
		String[] split = cyclicFragment.split("#");
		for (String fragment : split) {
			double weight = calculateFragmentWeight(fragment);
			sum += weight;
		}
		
		switch (type) {
			case DFBP:
				sum = sum + this.weights.get('2') - (18.0 * 2);
				break;
			case DISULFIDE:
				sum = sum + (this.weights.get('S') * 2 - 18.0) - (18.0 * 2);
				break;
			default:
				break;
			
		}
		
		return sum;
	}
	
	/**
	 * Returns a map of <index of connection, list of fragments>. This is a convenient map where the key is the index where the connection is in a linear
	 * fragment. If the key is -1, the list of linear fragments can be possible cyclic fragments.
	 * 
	 * @param connections
	 * @return
	 */
	private Map<Integer, List<List<Integer>>> findConnectionsInFragments(List<Integer> connections) {
		Map<Integer, List<List<Integer>>> connectionsInFragments = new HashMap<>();
		
		for (List<Integer> fragmentIndex : this.fragments) {
			int index = determineConnectionInFragment(fragmentIndex, connections);
			
			if (!connectionsInFragments.containsKey(index)) {
				List<List<Integer>> fragments = new ArrayList<>();
				fragments.add(fragmentIndex);
				connectionsInFragments.put(index, fragments);
			}
			else {
				connectionsInFragments.get(index)
					.add(fragmentIndex);
			}
		}
		
		return connectionsInFragments;
	}
	
	/**
	 * Given a linear fragment and the connection points of a peptide, returns the index where the connection is in the linear fragment. Returns -1 
	 * if the linear fragments has 2 connection fragments (possible cyclic fragment).
	 * 
	 * @param fragmentIndex
	 * @param connections
	 * @return
	 */
	private int determineConnectionInFragment(List<Integer> fragmentIndex, List<Integer> connections) {
		int result = -2;
		
		List<Integer> connectionIndices = new ArrayList<>();
		
		for (int index : fragmentIndex) {
			if (connections.contains(index)) {
				connectionIndices.add(index);
			}
		}
		
		if (connectionIndices.size() == 1) {
			result = connectionIndices.get(0);
		}
		else if (connectionIndices.size() == 2) {
			result = -1;
		}
		
		return result;
	}
	
	/**
	 * Given the peptide represented by its index and the peptide type, return the String representation 
	 * of the peptide.
	 * 
	 * @param peptideIndex
	 * @param type
	 * @return
	 */
	public String getPeptideStringRepresentation(List<Integer> peptideIndex, PeptideType type) {
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
		System.out.println("Finding all peptide fragments...");
		
		Map<Integer, List<Integer>> graph = this.peptide.getGraph();
		for (Map.Entry<Integer, List<Integer>> entry : graph.entrySet()) {
			int start = entry.getKey();
			
			walkGraph(start, graph);
		}
		
		System.out.println("Done");
		System.out.println();
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
	 * The fragmentIndex keeps track of the traversed nodes during the recursive backtracking step.
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
		if (!graph.containsKey(start)) { // Traversed the end of the peptide, the graph does not have a key for the end of the peptide
			return;
		}
		else if (hasDuplicates(fragmentIndex)) { // Illegal fragment. Fragments can't have the same amino acid in multiple positions
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
	
	private boolean hasDuplicates(List<Integer> fragmentIndex) {
		Set<Integer> setCheck = new HashSet<>(fragmentIndex);
		
		if (setCheck.size() < fragmentIndex.size()) {
			return true;
		}
		else {
			return false;
		}
	}
	
	/**
	 * Populate the amino acid weight map using the weights.csv file.
	 * 
	 * @return
	 */
	private Map<Character, Double> createAminoAcideWeightMap() {
		Map<Character, Double> aminoAcidWeight = new HashMap<>();
		
		InputStream weightsStream = FragmentAnalyzer.class.getResourceAsStream("/weights.csv");
		
		try {
			BufferedReader reader = new BufferedReader(new InputStreamReader(weightsStream));
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
