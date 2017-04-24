# BayFish - C++ version
Bayesian inference of transcription dynamics from population snapshots of smFISH

BayFish is a computational pipeline to infer kinetic parameters of gene expression from sparse single-molecule RNA fluorescence *in situ* hybridization (smFISH) data at multiple time points after induction. Given an underlying model of gene expression, BayFish uses a Markov Chanin Monte Carlo method to estimate the posterior probability of the model parameters and quantify the parameter uncertainty given the observed smFISH data.

NOTE: To use the MATLAB version, please refer to CODE_MATLAB/README.md

## Instructions
BayFish pipeline general description. A c++ compiler and the Armadillo Library (http://arma.sourceforge.net/) must be previously installed.

To run BayFish c++ version, first modify the first section of `main.cpp` to specify the model parameters:

```c++
    ////////////////////////////////////////////////////////////////////////////
    // Model parameters:
    const int T = 4;          // Number of time points.
    int myT[T] = {0,5,15,25}; // Time points; only factors of 5 are allowed.
    int N = 2;                // Number of promoter states (2 or 3).
    int maxM = 300;           // Maximum mRNA molecules.
    // Fixed biophysical parameters:
    Par pB;                   // ...in basal state.
    pB.d = 0.0462;            // Degradation rate (1/min).
    Par pS = pB;              // ...after stimulus.
    // Data:
    double a = 0;             // If N='3S', threshold to define third TS state.
    char* myDataCode = "Npas4"; // Code for data to load.
    // Metropolis Random Walk (MRW) parameters:
    int mrwI = 100000;        // Iterations.
    int mrwS = 7;             // Seed to use.
    // MRW sigma for parameter transition proposal in basal state:
    Par zigB;
    zigB.kON = 1e-5;
    zigB.kOFF = 1e-5;
    zigB.mu0 = 1e-5;
    zigB.mu = 0.01;
    // MRW sigma for parameter transition proposal in stimulus state:
    Par zigS;
    zigS.kON = 1e-5;
    zigS.kOFF = 1e-5;
    zigS.mu = 0.01;
    // MRW limits for parameter transition proposal in basal state:
    Par lB_m, lB_M;
    lB_m.kON = 1e-6;    lB_M.kON = 1e-2;
    lB_m.kOFF = 1e-4;   lB_M.kOFF = 1;
    lB_m.mu0 = 1e-5;    lB_M.mu0 = 1e-1;
    lB_m.mu = 1e-3;     lB_M.mu = 1;
    // MRW limits for parameter transition proposal in stimulus state:
    Par lS_m, lS_M;
    lS_m.kON = 1e-4;    lS_M.kON = 1;
    lS_m.kOFF = 1e-6;   lS_M.kOFF = 1e-2;
    lS_m.mu = 0.01;     lS_M.mu = 10;
    ////////////////////////////////////////////////////////////////////////////
```

1. Define data.
2. Define mathematical model.
3. Calculate probability distributions.
4. Compare model & data (i.e. calculate/define posteriors).
5. Run Metropolis-Hastings & obtain "best" kinetic parameters.
6. Explore sampling effect (i.e. synthetic data).

### (1) Define data:

The data must be in text format and named as myData_*myGene*_t*#*.txt, where *myGene* (e.g. `Npas4`) is the name or flag must assigned to the used data set, and each time point measure is a different file (e.g. `myData_Npas4_t5.txt`). Each file is the list of individual cell measurements (e.g. `[6.55,4.76,46]`) as tab-separated values; the first and second columns correspond to the measurement of active transcription sites (`TS1`, `TS2`), and the third column has the count of free mRNA molecules (`mRNA`).

```c++
    ////////////////////////////////////////////////////////////////////////////
    // Data:
    double a = 0;             // If N='3S', threshold to define third TS state.
    char* myDataCode = "Npas4"; // Code for data to load.
    ////////////////////////////////////////////////////////////////////////////
    // Load data matrix:
    myData x[T];
    for(int t = 0; t < T; t++)
        x[t].loadData(N,maxM,a,myDataCode,myT[t]);
```

### (2) Define mathematical model:

The mathematical model is defined depending on the number of promoter states (currently, the code supports either supports 2 or 3 promoter states; e.g. `2`), and the maximum number of free mRNA molecules to consider (e.g. `300`).

### (3) Calculate probability distributions



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
