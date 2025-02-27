# The stage at which neural noise arises determines whether context improves or degrades decision-making



[![Static Badge](https://img.shields.io/badge/bioRXiv-10.1101%2F2024.03.26.586597-red.svg?style=flat)](https://doi.org/10.1101/2024.03.26.586597)


__Authors: Bo Shen, Duc Nguyen, Jailyn Wilson, Paul Glimcher, Kenway Louie__

__New York University, Grossman School of Medicine__

---

## Contents

- [Overview](#overview)
- [Repo Contents](#repo-contents)
- [System Requirements](#system-requirements)
- [License](./LICENSE)
- [Citation](#citation)

# Overview

This is a repository for a manuscript under review, titled "Origins of noise in both improving and degrading decision making". This is an unpublished work. The contents may change during the revision.

# Repo Contents

- [Empirical data](./myData): Human choice data created by the current project. Please check the [README file](./myData/README.txt) for further information about the data structure and variables.
- [Data collection code](./DataCollectionCode): The Matlab code used for presenting choice information and collecting the participants' responses with Psychtoolbox 3.0.19.
- [Simulation code](./ModelSimulationCode): Matlab code for the modeling part of the paper. To replicate the simulation, please follow the main files [Fig1](./ModelSimulationCode/Fig1.m), [Fig2](./ModelSimulationCode/Fig2.m), [Fig4](./ModelSimulationCode/Fig4.m), [Supplementary Fig1ab](./ModelSimulationCode/SuppFig1ab.m), and [Supplementary Fig1d](./ModelSimulationCode/SuppFig1d.m).
- [Data analysis code](./BehavioralDataAnalysisCode.Rmd): `R` code for analyzing the empirical data. To replicate the data analysis, please change the directory of the data according to your local directory environment and follow the sections in the code.
- [Model fitting code](./ModelfitCode): Matlab code for model fitting can be found in [Modelfit.m](./ModelfitCode/ModelFit.m). Posterior check for fitting performance can be found in [Modelfit_PostCheck.m](./ModelfitCode/Modelfit_PostCheck.m).


# System Requirements

## Hardware and Software Requirements

The `R` code was developped with [R 4.3.3](https://www.r-project.org/#!). Excution of the `R` code requires only a standard computer with enough RAM. For minimal performance, this will be a computer with about 2 GB of RAM.

The `Matlab` code was developped with GPU parallel computation. For optimal performance, we recommend a computer with a graphic card. The code was developped on a Windows-11 laptop with the following specs:

RAM: 16.0 GB  
CPU: Intel(R) Core(TM) i7-10870H CPU @ 2.20GHz   
GPU: NVIDIA GeForce RTX 3070 Laptop GPU  

Matlab code was developped with Matlab 2023a (The Mathworks Inc., 2023).

### OS Requirements

The package development version is tested on *Mac OSX*, *Linux* and *Windows* operating systems. The 'R' code would be compatible with all of those operating systems when [`R`](https://www.r-project.org/#!) and [`RStudio`](https://posit.co/download/rstudio-desktop/) are properly set.

# Citation

For usage of the package and associated manuscript, please cite our [bioRXiv preprint](https://doi.org/10.1101/2024.03.26.586597).


