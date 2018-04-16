package com.github.juliomarcopineda.tests;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ThreadLocalRandom;

import com.github.juliomarcopineda.FragmentAnalyzer;
import com.github.juliomarcopineda.PeptideSerumStability;
import com.github.juliomarcopineda.peptide.Peptide;
import com.github.juliomarcopineda.peptide.PeptideType;

public class ConnectionsTest {
	public static void main(String[] args) {
		// Prevent equal connections
		// Prevent connections in decreasing order
		
		PrintStream out = null;
		try {
			out = new PrintStream(new FileOutputStream("../log.txt"));
		}
		catch (FileNotFoundException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		System.setOut(out);
		
		String peptideSequence = "CGYEQDPWGVRYWYGCKKKKB";
		System.out.println("TEST PEPTIDE: " + peptideSequence);
		System.out.println();
		
		for (int conn1 = 0; conn1 < peptideSequence.length() - 1; conn1++) {
			for (int conn2 = conn1 + 1; conn2 < peptideSequence.length(); conn2++) {
				if (conn1 == conn2) {
					continue;
				}
				
				System.out.println("TEST STAPLE CONNECTIONS: " + conn1 + " " + conn2);
				List<Integer> connections = Arrays.asList(conn1, conn2);
				
				for (PeptideType type : PeptideType.values()) {
					System.out.println("PEPTIDE TYPE: " + type);
					try {
						Map<Integer, List<Integer>> graph = PeptideSerumStability.createGraphStructure(peptideSequence, connections, type);
						
						Peptide peptide = new Peptide();
						if (type.equals(PeptideType.CUSTOM)) {
							double customWeight = ThreadLocalRandom.current()
								.nextDouble(0, 1000);
							
							System.out.println(customWeight);
							
							peptide.setCustomWeight(customWeight);
						}
						
						peptide.setSequence(peptideSequence);
						peptide.setConnections(connections);
						peptide.setType(type);
						peptide.setGraph(graph);
						
						FragmentAnalyzer analyzer = new FragmentAnalyzer(peptide);
						analyzer.findAllFragments()
							.measureAllFragmentWeights();
						
						System.out.println();
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
	}
}
