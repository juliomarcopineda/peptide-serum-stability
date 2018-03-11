package com.github.juliomarcopineda;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.AbstractMap;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import com.github.juliomarcopineda.peptide.Peptide;
import com.github.juliomarcopineda.peptide.PeptideType;

public class PeptideSerumStability {
	
	public static void main(String[] args) {
		
		String test = "bla";
		System.out.println(PeptideType.valueOf(test));
		
		System.exit(0);
		
		if (args[0].toLowerCase()
			.equals("input")) {
			
			if (args.length != 3) {
				System.out.println("Please add the right number of arguments for the choice \"input\"");
				System.exit(1);
			}
			
			String inputFile = args[1];
			String outputFile = args[2];
		}
		else if (args[0].toLowerCase()
			.equals("interactive")) {
			
			printIntro();
			
			Peptide peptide = new Peptide();
			
			try (BufferedReader br = new BufferedReader(new InputStreamReader(System.in))) {
				// Set the peptide sequence
				System.out.println("Please enter the peptide sequence. If peptide is not linear, do not include the linkers");
				System.out.print("Peptide sequence: ");
				
				String peptideSequence = br.readLine();
				peptide.setSequence(peptideSequence);
				
				// Set the peptide type
				System.out.println("Please enter the peptide type (linear, amide, DFBP, disulfide)");
				System.out.print("Peptide Type: ");
				
				String typeString = br.readLine();
				try {
					PeptideType type = PeptideType.valueOf(typeString.toUpperCase());
					peptide.setType(type);
					
					// Set the conenctions of the peptide
					List<Integer> connections = new ArrayList<>();
					if (!type.equals(PeptideType.LINEAR)) {
						System.out.println("Peptide type is cyclic.");
						System.out.println("Please enter the indices where the cyclic conenctions are. Seperate the index with commas");
						
						String connectionsString = br.readLine();
						String[] split = connectionsString.split(",");
						
						for (String indexString : split) {
							connections.add(Integer.parseInt(indexString));
						}
					}
					peptide.setConnections(connections);
					
					// Set the graph structure of the peptide
					Map<Integer, List<Integer>> graph = createGraphStructure(peptideSequence, connections, type);
					peptide.setGraph(graph);
					
				}
				catch (IllegalArgumentException e) {
					System.out.println("Please input a valid peptide type (linear, amide, DFBP, disulfide)");
					System.exit(1);
				}
				
				// Create FragmentAnalyzer and then start interactive session
				FragmentAnalyzer analyzer = new FragmentAnalyzer(peptide);
				interactiveSession(analyzer);
			}
			catch (IOException e) {
				e.printStackTrace();
			}
		}
	}
	
	/**
	 * Given a fragment analyzer, runs the interactive session for the user in the command line.
	 * 
	 * @param analyzer
	 */
	private static void interactiveSession(FragmentAnalyzer analyzer) {
		
	}
	
	/**
	 * Creates the graph structure reperesntation of a peptide. This graph is represented as a map from the index and then its connections.
	 * 
	 * @param peptideSequence
	 * @param connections
	 * @param type
	 * @return
	 */
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
	
	/**
	 * Intro message when using the interactive mode of the program.
	 * 
	 * @param peptide
	 */
	private static void printIntro() {
		System.out.println("-------------------------");
		System.out.println("Peptide Serum Stability");
		System.out.println("-------------------------");
		System.out.println();
		System.out.println("This program will generate all possible fragments of the given input peptide sequence");
		System.out.println("and will calculate their corresponding theoretical molecular weights.");
		System.out.println();
	}
}
