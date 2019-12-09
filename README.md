# SpectrumAI
R script to automatically inspect MS2 spectra of single-subsititution peptide identifications.

Please cite this paper when you have used the SpectrumAI in your publications.

Zhu Y, Orre LM, Johansson HJ, Huss M, Boekel J, Vesterlund M, Fernandez-Woodbridge A, Branca RMM, Lehtio J: Discovery of coding regions in the human genome by integrated proteogenomics analysis workflow. Nat Commun 2018, 9(1):903.  [PMID: 29500430](https://www.ncbi.nlm.nih.gov/pubmed/29500430)

SpectrumAI.R is the script to curate PSMs of single-substitution peptides.

plot_mirror_image.R is the script to draw MS2 spectrum mirror image of synthetic and endogenous peptide.

plot_annotated_spectrum.R is the script to draw annotated spectra of any PSM. (matched b,y ions to theorectial spectrum will be highlighted)


# Prerequisite R library

Requires R libraries `mzR` , `protViz`, and `MSnbase`.

Install them from Bioconductor

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
	install.packages("BiocManager")

## mzR will be installed automatically as MSnbase depends on it
BiocManager::install(c("MSnbase", "protViz"))
```

# How does SpectrumAI work?
![My image](https://github.com/yafeng/SpectrumAI/blob/master/image/sequence_example.png)

Assume a 12-amino-acid peptide is identified with single substitution at 8th residue, in order to pass SpectrumAI, it must have matched MS2 peaks (within 10 ppm fragment ion mass tolerance) from at least one of the following groups: b7&b8, y4&y5, y4&b7 or y5&b8. Second, the sum intensity of the supporting flanking MS2 ions must be larger than the median intensity of all fragmentation ions. An exception to these criteria is made when the substituted amino acid has a proline residue to its N-terminal side. Because CID/HCD fragmentation at the C-terminal side of a proline residue is thermodynamically unfavored, SpectrumAI only demands the presence of any b or y fragment ions containing substituted amino acids, in this case, b8 to b11, y5 to y11.


# Evaluation and a false single amino acid sequence variant spotted by SpectrumAI

![My image](https://github.com/yafeng/SpectrumAI/blob/master/image/SpectrumAI.png)

# how to use
## To curate a list of variant peptide identifications (support only MSGF+ search outputs)
1. Prepare the list of variant peptide PSMs table in tabular format. It should contain at least the following columns with exactly same names (the order can be different):

	"SpectraFile" -  The filename of spectra file in which the variant peptide PSM is identified.

	"ScanNum" - The scan number

	"Peptide" - The peptide sequence. Modifications are noted. For instance, in MSGF+ search engine output, "M+15.995GYEEAE", M+15.995 is used to indicate oxidation on Methonine.

	"sub_pos" -  position in peptide indicate which amino acid is substituted. index starts from 1.

2. Open SpectrumAI.R script in R Studio.

3. A few variables to set in advance.
```{r}
setwd()  #set your working directory
mzml_path = "" # set file path to which raw files are located
infile_name =""  # PSM table file name
outfile_name =""  #set corresponding output file name
Frag.ions.tolerance= 0.02 # 0.02 Da tolerance for MS2 fragment ions mass accuracy.
relative=FALSE  # set TRUE if ppm value is used for Frag.ions.tolerance
```

**ERROR! An unexpected error regarding mass shift of modifications on label-free peptides was found and corrected. Changes made on 2019-12-09. Commit 7e748e4. If you ran SpectrumAI on label-free peptides, please rerun it. /my sincere apology**

**Warning! SpectrumAI has not been tested in low-resolution MS data, it may not perform as expected.**

4. After you set everything, ctrl-A to select all codes and click run.

## To generate mirror plots for a list of synthetic peptides (to be updated)
"plot_mirror_image.R" is the R script to generate the mirror plot in Figure 3C.
The mirror plot is generally used to assess if the peptide spectrum match of your candidate
peptide is correct. By generating a clean MS spectrum of synthetic peptide sequence, you
expect it to resemble the original spetrum of your candidate peptide as a mirrored image
without any mismatch of major b and y ions.

To make this mirror plot for a list of candidate peptides, you need have these files
at hand.

PSM table of candidate peptides including spectra file information, scan number, peptide sequences.

PSM table of synthetic peptides including spectra file information, scan number, peptide sequences.

All mzML spectra file that are present in candidate peptides PSM table are stored in one folder.

All mzML spectra file that are present in synthetic peptides PSM table are stored in one folder.

## To generate annotated spectra of peptide identifications (to be updated)
