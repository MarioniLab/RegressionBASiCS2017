# RegressionBASiCS2017

This repository contains scripts for processing raw data, quality control, running the model and performing downstream analysis for Eling *et al.*, Robust expression variability testing reveals heterogeneous T cell responses.

* [Preprocessing](../master/Preprocessing/) contains a Makefile to map raw data from [L&ouml;nnberg *et al.*](http://immunology.sciencemag.org/content/2/9/eaal2192) to the mouse reference genome.
  This folder also contains a file for data preparation and quality control of all dataset used in this study.

* All scripts to run the model can be found in [ChainGeneration](../master/ChainGeneration).

* To recreate the analysis and reproduce the figures, scripts in [Figures](../master/Figures) can be used.

The preprint can be found at:

[Robust expression variability testing reveals heterogeneous T cell responses](https://www.biorxiv.org/content/early/2017/12/21/237214)
