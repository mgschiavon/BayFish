# BayFish - C++ version
Bayesian inference of transcription dynamics from population snapshots of smFISH

BayFish is a computational pipeline to infer kinetic parameters of gene expression from sparse single-molecule RNA fluorescence *in situ* hybridization (smFISH) data at multiple time points after induction. Given an underlying model of gene expression, BayFish uses a Markov Chanin Monte Carlo method to estimate the posterior probability of the model parameters and quantify the parameter uncertainty given the observed smFISH data.

NOTE: To use the MATLAB version, please refer to CODE_MATLAB/README.md

## Instructions
BayFish pipeline general description. A name or flag must be assigned to the used data set (*myGene*, e.g. `Npas4`).

1. Define data.
2. Define mathematical model.
3. Calculate probability distributions.
4. Compare model & data (i.e. calculate/define posteriors).
5. Run Metropolis-Hastings & obtain "best" kinetic parameters.
6. Explore sampling effect (i.e. synthetic data).

### (1) Define data:

(1.1) 

### (2) Define mathematical model:

(2.1) 

### (3) Calculate probability distributions

(3.1) 

### (4) Compare model & data:

(4.1) 

### (5) Run Metropolis-Hastings & obtain "best" kinetic parameters:

(5.1) 

### (6) Explore sampling effect (i.e. synthetic data):

(6.1) 

## Referencing

If you use this code or the data associated with it please cite:

Gómez-Schiavon *et al.* (2017); https://doi.org/10.1101/109603.

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
