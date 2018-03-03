package com.github.juliomarcopineda;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.AbstractMap;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * This class parses through the designated input text format for the peptide serum stability analysis. The general format as follows:
 * 
 * Line 0: [peptide sequence] [peptide type] [optional: first connection index] [optional: second connection index] ... (more indices if desired).
 * Line 1: Mass spectrometry data delimited with white space.
 * 
 * More lines can be added if more sequences want to be analyzed.
 * 
 * The peptide type can have the following valid options: linear, disulfide, dfbp and amide.
 * If the peptide type is not linear, the indices afterwards must be even in number. The index is assumed to be zero-index.
 * 
 * @author Julio Pineda
 *
 */
public class InputParser {
	private String inputFile;
	
	private List<Peptide> peptides;
	
	public InputParser(String inputFile) {
		this.inputFile = inputFile;
	}
	
	public InputParser() {
		
	}
	
	/**
	 * Begins the process of parsing through the input text file.
	 * 
	 * Throws an IllegalArgumentException if the first line for each peptide contains an odd number of arguments. This suggests that the number of connection
	 * indices entered was incorrect. The number of indices must be even to have a valid input.
	 */
	public void parse() {
		// Initialize array list of peptides
		List<Peptide> peptides = new ArrayList<>();
		
		try (BufferedReader reader = new BufferedReader(new FileReader(inputFile))) {
			
			String line;
			int lineNumber = 0;
			while ((line = reader.readLine()) != null) {
				if (lineNumber % 2 == 0) {
					String[] split = line.split("\\s+");
					
					if (split.length % 2 != 0) {
						throw new IllegalArgumentException();
					}
					
					// Extract input arguments from text file
					String peptideSequence = split[0];
					PeptideType type = PeptideType.valueOf(split[1].toUpperCase());
					List<Integer> indexConnections = new ArrayList<>();
					
					if (!type.equals(PeptideType.LINEAR)) {
						
						for (int i = 2; i < split.length; i++) {
							indexConnections.add(Integer.parseInt(split[i]));
						}
					}
					
					// Build graph structure from sequence and index connections
					Map<Integer, List<Integer>> graph = createGraphStructure(peptideSequence, indexConnections, type);
					
					// Create Peptide object from data above
					Peptide peptide = new Peptide();
					peptide.setSequence(peptideSequence);
					peptide.setType(type);
					// peptide.setGraph(graph);
					
					// Add to list of peptides
					peptides.add(peptide);
				}
				else {
					
				}
				
				lineNumber++;
			}
			
		}
		catch (IOException e) {
			e.printStackTrace();
		}
		
		this.peptides = peptides;
	}
	
	/**
	 * Given the peptide sequence, indices of connections and the peptide type, creates a directed network representation of the peptide where 
	 * each amino acid (symbol) is a node, and each node has a directed edge pointing towards its connection in the sequence.
	 * 
	 * @param peptideSequence
	 * @param indexConnections
	 * @return
	 */
	private Map<Integer, List<Integer>> createGraphStructure(String peptideSequence, List<Integer> connections, PeptideType type) {
		
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
	
	public static void main(String[] args) {
		String seq = "YEQDPWGVKK";
		
		InputParser test = new InputParser();
		List<Integer> connections = Arrays.asList(2, 8);
		
		Map<Integer, List<Integer>> graph = test.createGraphStructure(seq, connections, PeptideType.DISULFIDE);
		
		System.out.println("SEQUENCE: " + seq);
		System.out.println();
		
		graph.entrySet()
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
