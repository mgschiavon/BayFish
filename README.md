# BayFish
Bayesian inference of transcription dynamics from population snapshots of smFISH

BayFish is a computational pipeline to infer kinetic parameters of gene expression from sparse single-molecule RNA fluorescence *in situ* hybridization (smFISH) data at multiple time points after induction. Given an underlying model of gene expression, BayFish uses a Markov Chanin Monte Carlo method to estimate the posterior probability of the model parameters and quantify the parameter uncertainty given the observed smFISH data.

## Directories

The `CODE_MATLAB` directory contains MATLAB code used to process the data (`DATA_*.m`), run the pipeline (`SIM_*.m`), and analyze (`FIG_*.m`) the BayFish results.

The `CODE_CPP` directory contains C++ code used to run the BayFish pipeline (`main.cpp`), as well as the required function files (`*.h`).

The `DATA` directory contains an example of data to be processed: smFISH measurements of the neuronal activity inducible gene *Npas4* in primary neurons (see Gómez-Schiavon *et al.*, 2017; https://doi.org/10.1101/109603).

## Instructions
To use the MATLAB version, please refer to CODE_MATLAB/README.md
To use the C++ version, please refer to CODE_CPP/README.md

## Referencing

If you use this code or the data associated with it please cite:

Gómez-Schiavon *et al.* (2017); https://doi.org/10.1101/109603.

Latest release: [![DOI](https://zenodo.org/badge/82009613.svg)](https://zenodo.org/badge/latestdoi/82009613)

## Copyright

(C) Copyright 2017 Mariana Gómez-Schiavon

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
