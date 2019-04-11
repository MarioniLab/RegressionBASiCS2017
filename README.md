# RegressionBASiCS2017

This repository contains scripts for processing raw data, quality control, running the model and performing downstream analysis for Eling *et al.*, Robust expression variability testing reveals heterogeneous T cell responses.

* [Preprocessing](../master/Preprocessing/) contains a Makefile to map raw data from [L&ouml;nnberg *et al.*](http://immunology.sciencemag.org/content/2/9/eaal2192) to the mouse reference genome.
  This folder also contains a file for data preparation and quality control of all dataset used in this study.

* All scripts to run the model can be found in [ChainGeneration](../master/ChainGeneration).

* To recreate the analysis and reproduce the figures, scripts in [Figures](../master/Figures) can be used.

The publication can be found at:

[Correcting the Mean-Variance Dependency for Differential Variability Testing Using Single-Cell RNA Sequencing Data](https://www.cell.com/cell-systems/fulltext/S2405-4712(18)30278-3)
