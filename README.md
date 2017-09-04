# SpectrumAI
R script to automatically inspect MS2 spectra of single-subsititution peptide identifications


SpectrumAI.R is the script to curate PSMs of single-substitution peptides.

plot_mirror_image.R is the script to draw MS2 spectrum mirror image of synthetic and endogenous peptide.

plot_annotated_spectrum.R is the script to draw annotated spectra of any PSM. (matched b,y ions to theorectial spectrum will be highlighted) 


Prerequisite R library

requires R libraries "mzR" and "MSnbase".

Install them from Bioconductor

source("https://bioconductor.org/biocLite.R")

biocLite("mzR")

biocLite("MSnbase")


How does SpectrumAI works?

![My image](https://github.com/yafeng/SpectrumAI/sequence_example.png)

Imagine a 12-amino-acid peptide is identified with single substitution at 8th residue, in order to pass SpectrumAI, it must have matched MS2 peaks (within 10 ppm fragment ion mass tolerance) from at least one of the following groups: b7&b8, y4&y5, y4&b7 or y5&b8. Second, the sum intensity of the supporting flanking MS2 ions must be larger than the median intensity of all fragmentation ions. An exception to these criteria is made when the substituted amino acid has a proline residue to its N-terminal side. Because CID/HCD fragmentation at the C-terminal side of a proline residue is thermodynamically unfavored, SpectrumAI only demands the presence of any b or y fragment ions containing substituted amino acids, in this case, b8 to b11, y5 to y11. 


Evaluation and a false single amino acid sequence variant spotted by SpectrumAI

![My image](https://github.com/yafeng/SpectrumAI/SpectrumAI.png)
