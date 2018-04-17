package com.github.juliomarcopineda.tests;

import java.util.Arrays;
import java.util.List;
import java.util.Map;

import com.github.juliomarcopineda.FragmentAnalyzer;
import com.github.juliomarcopineda.PeptideSerumStability;
import com.github.juliomarcopineda.peptide.Peptide;
import com.github.juliomarcopineda.peptide.PeptideType;

/**
 * This class prints out all the branched fragments for a peptide to check if the branched fragment builder 
 * is working.
 * 
 * @author Julio Pineda
 *
 */
public class BranchedFragmentTest {
	
	public static void main(String[] args) {
		Peptide peptide = new Peptide();
		
		String peptideSequence = "CGYEQDPWGVRYWYGCKKKKB";
		PeptideType type = PeptideType.CUSTOM;
		double customWeight = 383.32;
		List<Integer> connections = Arrays.asList(4, 19);
		Map<Integer, List<Integer>> graph = PeptideSerumStability.createGraphStructure(peptideSequence, connections, type);
		
		peptide.setSequence(peptideSequence);
		peptide.setType(type);
		peptide.setConnections(connections);
		peptide.setCustomWeight(customWeight);
		peptide.setGraph(graph);
		
		FragmentAnalyzer analyzer = new FragmentAnalyzer(peptide);
		analyzer.findAllFragments()
			.measureAllFragmentWeights();
		
		analyzer.getFragmentWeights()
			.keySet()
			.stream()
			.filter(s -> s.contains("#"))
			.forEach(System.out::println);
	}
	
}
