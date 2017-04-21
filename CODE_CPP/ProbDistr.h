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
 * PROBDISTR: Calculates probability distributions for a given transition 
 *            matrix and time
 * 
 * Created by Mariana Gómez-Schiavon 
 * June 2016 
 * 
 * ProbDistr : Stationary probability distribution and protein distribution 
 *  dynamics.
 * 
 *  mat Pss(mat A) : Given the transition matrix A, returns the stationary 
 *      probability distribution vector.
 * 
 *  double logL(Mat<int> x, mat P) : Calculate the log-likelihood of observing 
 *      the data x given the probability distribution vector P.
 * 
 *  mat LxT(ModelStruct *ms, myData *x, Par pB, Par pS, int T, int *myT) : 
 *      Iterate over time points myT[T] to estimate the log-likelihood of 
 *      observing the data x under the model ms with paramters pB in basal 
 *      state and pS after stimulus, and returns a matrix L(1,T).
 * 
 */

#ifndef PROBDISTR_H
#define PROBDISTR_H

#include <stdio.h>
#include <cmath>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <algorithm>
#include <vector>
#include <armadillo>
#include "Data.h"
#include "Model.h"

using namespace std;
using namespace arma;

mat Pss(mat A)
{
    cx_vec eigval;
    cx_mat eigvec;
    eig_gen(eigval,eigvec,A); 
    
    mat Pss = abs(eigvec.col((abs(eigval)).index_min()));
    mat n = sum(Pss);
    Pss = Pss/n(0,0);
    
    return Pss;
}

double logL(Mat<int> x, mat P)
{
    P.reshape(size(x));
    double L = sum(sum(x % trunc_log(P)));
    return L;
}

mat LxT(ModelStruct *ms, myData *x, Par pB, Par pS, int T, int *myT)
{
    mat L(1,T);
    mat Ab = ms->TransM(pB);
    mat As = ms->TransM(pS);
    mat At5 = expmat(As*5);
    mat P[T];
    for(int t = 0; t < T; t++)
    {
        if(t==0)
        {
            P[t] = Pss(Ab);
        }
        else
        {
            int temp = myT[t] - myT[t-1];
            P[t] = P[t-1];
            while(temp > 0)
            {
                P[t] = At5*P[t];
                temp -=5;
            }
        }
        L(0,t) = logL(x[t].data, P[t]);
    }
    return L;
};

#endif /* PROBDISTR_H */

