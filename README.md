# peptide-serum-stability
Please visit the [wiki](https://github.com/juliomarcopineda/peptide-serum-stability/wiki) for user guides and other information.
# Background
peptide-serum-stability assists in the analysis of peptide fragments in peptide serum stability studies.

Serum stability studies assess the ability of a peptide to resist degradation by serum proteases. During the study, the peptide is incubated in serum and then aliquots are taken from the serum periodically. These aliquots are analyzed using mass spectrometry (such as MALDI-ToF) and measurements from the mass spectrometry are produced.

Instead of manually calculating molecular weights of peptide fragments and matching to the mass spectrometry measurements, peptide-serum-stability performs this process programatically. The program generates all possible peptide fragments with the given peptide sequence, and then calculates the fragments' theoretical molecular weights. Then, for every mass spectrometry measurement, the program will suggest possible fragments that match up to this data input.

The current implementation of peptide-serum-stability can analyze linear and certain types of cyclic peptides. The cyclic peptides that can be analyzed are amide, disulfide and DFBP.

# Releases
If you need the runnable JAR file "stability.jar", please look at the releases page: 

https://github.com/juliomarcopineda/peptide-serum-stability/releases
