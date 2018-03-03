package com.github.juliomarcopineda;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

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
					Map<Character, List<Character>> graph = createGraphStructure(peptideSequence, indexConnections, type);
					
					// Create Peptide object from data above
					Peptide peptide = new Peptide();
					peptide.setSequence(peptideSequence);
					peptide.setType(type);
					peptide.setGraph(graph);
					
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
	private Map<Character, List<Character>> createGraphStructure(String peptideSequence, List<Integer> indexConnections, PeptideType type) {
		
		Map<Character, List<Character>> graph = new LinkedHashMap<>();
		
		// Create the graph structure of the peptide base sequence
		for (int i = 0; i < peptideSequence.length() - 1; i++) {
			char source = peptideSequence.charAt(i);
			char target = peptideSequence.charAt(i + 1);
			
		}
		
		return graph;
	}
}
