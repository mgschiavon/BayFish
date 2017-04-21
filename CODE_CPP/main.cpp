/*
 * (C) Copyright 2017 Mariana GÃ³mez-Schiavon
 *
 *    This file is part of BayFish.
 *
 *    BayFish is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    BayFish is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with BayFish.  If not, see <http://www.gnu.org/licenses/>.
 *
 * BayFish pipeline
 * SIMULATE: Run Metropolis-Hastings & obtain "best" kinetic parameters.
 * 
 * Created by Mariana Gómez-Schiavon 
 * June 2016 
 * 
 * Simulations : Generate model, load data, and run Metropolis Random Walk.
 * 
 */

#include <iostream>
#include <string>
#include <armadillo>
#include "Data.h"
#include "Model.h"
#include "ProbDistr.h"
#include "MRW.h"
#include <iomanip>


using namespace std;
using namespace arma;

// Armadillo documentation is available at:
// http://arma.sourceforge.net/docs.html

int main(int argc, char** argv)
{
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
    char* myDataCode = "Fos"; // Code for data to load.
    // Metropolis Random Walk (MRW) parameters:
    int mrwI = 3;             // Iterations.
    int mrwS = 1;             // Seed to use.
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
    
    // Define model structure:
    ModelStruct ms(N,maxM);
    // Load data matrix:
    myData x[T];
    for(int t = 0; t < T; t++)
        x[t].loadData(N,maxM,a,myDataCode,myT[t]);
    // Create MRW structure:
    arma_rng::set_seed(mrwS); // Set seed for random number generator.
    mrwPar mrw;
    mrw.pB = mrw.ParToMat(pB);
    mrw.pS = mrw.ParToMat(pS);
    mrw.zigB = mrw.ParToMat(zigB);
    mrw.zigS = mrw.ParToMat(zigS);
    
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
    
    MRWp.close();
    MRWl.close();
  return 0;
  }
