package com.github.juliomarcopineda.peptide;

import java.util.List;
import java.util.Map;

/**
 * This class represents the data structure required for a Peptide object for this project. This class contains the peptide sequence, the peptide type 
 * determined by the PeptideType enum and the graph structure representing the sequence.
 * 
 * @author Julio Pineda
 *
 */
public class Peptide {
	private String sequence;
	private PeptideType type;
	private Map<Integer, List<Integer>> graph;
	private List<Double> massSpecData;
	
	public List<Double> getMassSpecData() {
		return massSpecData;
	}
	
	public void setMassSpecData(List<Double> massSpecData) {
		this.massSpecData = massSpecData;
	}
	
	public String getSequence() {
		return sequence;
	}
	
	public void setSequence(String sequence) {
		this.sequence = sequence;
	}
	
	public PeptideType getType() {
		return type;
	}
	
	public void setType(PeptideType type) {
		this.type = type;
	}
	
	public Map<Integer, List<Integer>> getGraph() {
		return graph;
	}
	
	public void setGraph(Map<Integer, List<Integer>> graph) {
		this.graph = graph;
	}
	
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((graph == null) ? 0 : graph.hashCode());
		result = prime * result + ((sequence == null) ? 0 : sequence.hashCode());
		result = prime * result + ((type == null) ? 0 : type.hashCode());
		return result;
	}
	
	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		Peptide other = (Peptide) obj;
		if (graph == null) {
			if (other.graph != null)
				return false;
		}
		else if (!graph.equals(other.graph))
			return false;
		if (sequence == null) {
			if (other.sequence != null)
				return false;
		}
		else if (!sequence.equals(other.sequence))
			return false;
		if (type != other.type)
			return false;
		return true;
	}
}
