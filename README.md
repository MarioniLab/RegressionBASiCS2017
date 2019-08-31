[![doi 10.5281/zenodo.3382836](https://zenodo.org/badge/DOI/10.5281/zenodo.3382836.svg)](https://doi.org/10.5281/zenodo.3382836)

# RegressionBASiCS2017

This repository contains scripts for processing raw data, quality control, running the model and performing downstream analysis for [Eling *et al.* (2018) Correcting the Mean-Variance Dependency for Differential Variability Testing Using Single-Cell RNA Sequencing Data](https://www.cell.com/cell-systems/fulltext/S2405-4712(18)30278-3).

* [Preprocessing](../master/Preprocessing/) contains a Makefile to map raw data from [L&ouml;nnberg *et al.*](http://immunology.sciencemag.org/content/2/9/eaal2192) to the mouse reference genome.
  This folder also contains a file for data preparation and quality control of all dataset used in this study.

* All scripts to run the model can be found in [ChainGeneration](../master/ChainGeneration).

* To recreate the analysis and reproduce the figures, scripts in [Figures](../master/Figures) can be used.

# Correction

While revising code used to analyze the data described in Eling et al. (2018) for inclusion in a Bioconductor workflow associated to the BASiCS package, we discovered a problem with the analyses of datasets generated using the Fluidigm C1 system (Antolović et al., 2017; Lönnberg et al., 2017; Martinez-Jimenez et al., 2017). 
Specifically, we observed that the number of spike-in molecules present in the cell lysis volume had been miscomputed in our original analyses, such that all true spike-in numbers were inadvertently scaled by the same constant factor. 
Where available, these quantities are used as an input for BASiCS, and therefore, some of the outputs originally reported in Eling et al. were incorrect.

In principle, the original error in calculating the exact number of spike-in molecules per reaction only scales the arbitrary units in which gene expression is measured and thus should not alter any interpretation or downstream analysis. 
However, when we re-analyzed these data to correct Eling et al., we noted that mis-scaling led to a poor initialization of the MCMC sampler, which led to less stable estimates of the mean and dispersion for lowly expressed genes. 
In essence, more iterations were needed to achieve optimal good convergence. 
Consequently, using the correct input spike-in numbers leads to changes in the set of differentially expressed and variable genes identified.

The calculation of the ERCC spike-in molecules has been corrected in the [ChainGeneration](../master/ChainGeneration/) folder and the code to generate the corrected Figures can be found in the [Corrections](../master/Corrections/) folder.
