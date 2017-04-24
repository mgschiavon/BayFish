# BayFish - C++ version
Bayesian inference of transcription dynamics from population snapshots of smFISH

BayFish is a computational pipeline to infer kinetic parameters of gene expression from sparse single-molecule RNA fluorescence *in situ* hybridization (smFISH) data at multiple time points after induction. Given an underlying model of gene expression, BayFish uses a Markov Chain Monte Carlo method to estimate the posterior probability of the model parameters and quantify the parameter uncertainty given the observed smFISH data.

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

Then, compile `main.cpp`:

```
g++ main.cpp -l armadillo -o RunMRW.exe
```

where `g++` is the compiler being used, `-l armadillo` specifies the Armadillo library is going to be used, and `-o RunMRW.exe` is the output/executable file. Finally, run `RunMRW.exe`. Two output files will be produce, a list of parameters per iteration (`*_Par.dat`) and a list of log-likelihood per time point per iteration (`*_logL.data`). See details in the following sections.

### Define data:

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

The mathematical model is defined depending on the number of promoter states (currently, the code supports either supports 2 or 3 promoter states; e.g. `N = 2`), and the maximum number of free mRNA molecules to consider (e.g. `maxM = 300`).

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
////////////////////////////////////////////////////////////////////////////

    // Define model structure:
    ModelStruct ms(N,maxM);
```

The `Par` class (see `Model.h`) includes the following parameters:

```c++
class Par
{
public:
    double kON, kOFF;
    double kONs, kOFFs;
    double mu0, mu, muS;
    double d;
    
    Par()
    {
        kON = 0;	// Promoter activation rate (OFF->ON)
        kOFF = 0;	// Promoter deactivation rate (ON->OFF)
        kONs = 0;	// Promoter "super" activation rate (ON->ONs)
        kOFFs = 0;	// Promoter "super" deactivation rate (ONs->ON)
        mu0 = 0;	// mRNA synthesis rate of promoter in OFF state
        mu = 0;		// mRNA synthesis rate of promoter in ON state
        muS = 0;	// mRNA synthesis rate of promoter in ONs state
        d = 0;		// mRNA degradation rate
    }
};
```

### (3) Metropolis-Hastings

All fixed parameters values must be specified in the `pB` and `pS` structures. The minimum and maximum values for the initial random parameter values, as well as the variance for the proposals must be specified for each parameter being fitted.

```c++
////////////////////////////////////////////////////////////////////////////
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

    // Create MRW structure:
    arma_rng::set_seed(mrwS); // Set seed for random number generator.
    mrwPar mrw;
    mrw.pB = mrw.ParToMat(pB);
    mrw.pS = mrw.ParToMat(pS);
    mrw.zigB = mrw.ParToMat(zigB);
    mrw.zigS = mrw.ParToMat(zigS);
	
	    // First iteration:
    mrw.initPar(mrw.ParToMat(lB_m),mrw.ParToMat(lB_M),  // Initialize parameters.
            mrw.ParToMat(lS_m),mrw.ParToMat(lS_M));
    MRWp << 1 <<' ' << "B" << ' ';
    mrw.pB.raw_print(MRWp);
    MRWp << 1 << ' ' << "S" << ' ';
    mrw.pS.raw_print(MRWp);
    
    mat L(1,T);     // Log-likelihood per time point.
    L = LxT(&ms,x,mrw.MatToPar(mrw.pB),mrw.MatToPar(mrw.pS),T,myT);
    MRWl << 1 << ' ';
    L.raw_print(MRWl);
    
    // Iterate:
    for(int i = 2; i <= mrwI; i++)
    {
        mat ptB = mrw.ptB();
        mat ptS = mrw.ptS(ptB);
        if(min(min(join_cols(ptB+(mrw.pB==0),ptS+(mrw.pS==0))))>1e-8)
        {
            mat Lt(1,T);
            Lt = LxT(&ms,x,mrw.MatToPar(ptB),mrw.MatToPar(ptS),T,myT);
            
            // If proposal is accepted, update system:
            mat r = randu(1,1);
            mat alpha = min(join_rows(ones(1,1),exp(sum(Lt,1)-sum(L,1))),1);
            if(r(0) <= alpha(0))
            {
                mrw.pB = ptB;
                mrw.pS = ptS;
                L = Lt;
            }
        }
        MRWp << i <<' ' << "B" << ' ';
        mrw.pB.raw_print(MRWp);
        MRWp << i << ' ' << "S" << ' ';
        mrw.pS.raw_print(MRWp);
        MRWl << i << ' ';
        L.raw_print(MRWl);
    }
```

### (4) Output files:

The output files are named accordingly to the data set used, the model (i.e. the number of promoter states, `N`, and the maximum mRNA number, `maxM`), and the random seed used (`s`); for example: `MRW_Npas4_N2(300)_s7_Par.dat`.

```c++
    // OUTPUT FILES
    // Parameters:
    char myOutputFile[255];
    strcpy (myOutputFile,"MRW_");
    strcat (myOutputFile,myDataCode);
    strcat (myOutputFile,"_N%d(%d)_s%d_Par.dat");
    sprintf(myOutputFile,myOutputFile,N,maxM,mrwS);
    ofstream MRWp(myOutputFile,ios::out);
    MRWp.precision(4);
    MRWp << "Iteration" << ' ' << "[B/S]" << ' ';
    MRWp << "kON" << ' ' << "kOFF" << ' ' << "kONs" << ' ' << "kOFFs" << ' ';
    MRWp << "mu0" << ' ' << "mu" << ' ' << "muS" << ' ' << "d" << endl;
    // Log-likelihoods:
    strcpy (myOutputFile,"MRW_");
    strcat (myOutputFile,myDataCode);
    strcat (myOutputFile,"_N%d(%d)_s%d_logL.dat");
    sprintf(myOutputFile,myOutputFile,N,maxM,mrwS);
    ofstream MRWl(myOutputFile,ios::out);
    MRWl.precision(6);
    MRWl << "Iteration" << ' ';
    for(int t = 0; t < (T-1); t++)
        MRWl << "logL[" << t << "]" << ' ';
    MRWl << "logL[" << (T-1) << "]" << endl;
```

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
