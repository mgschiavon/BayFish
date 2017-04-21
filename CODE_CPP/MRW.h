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
 * MRW: Functions and steps required for the Metropolis Random Walk algorithm
 *      to estimate the posterior probability distribution under a specific 
 *      model.
 * 
 * Created by Mariana Gómez-Schiavon 
 * June 2016 
 * 
 * MRW : My Metropolis Random Walk (MRW) algorithm.
 * 
 *  class mrwPar : Structure to follow the MRW progress.
 *      mat zigB : Variance for parameter proposals in basal state.
 *      mat zigS : Variance for parameter proposals in stimulus state.
 *      mat pB : Current parameters in basal state.
 *      mat pS : Current parameters in stimulus state.
 * 
 *      mat ParToMat(Par p) : Translates a Par structure to a vector.
 * 
 *      Par MatToPar(mat m) : Translates a vector to a Par structure.
 * 
 *      void initPar(mat lB_m, mat lB_M, mat lS_m, mat lS_M) : Initialize 
 *          parameters to be fitted in the MRW by choosing a uniformly 
 *          distributed random number between [l_m,l_M]. Notice that when 
 *          the parameter does not change with stimulus, the ptB value is 
 *          copied.
 * 
 *      mat ptB() : Calculates the next proposal parameters in basal state to be 
 *          evaluated by the Metropolis algorithm.
 * 
 *      mat ptS(mat ptB) : Calculates the next proposal parameters in stimulus 
 *          state to be evaluated by the Metropolis algorithm. Notice that when 
 *          the parameter does not change with stimulus, the ptB value is 
 *          copied.
 */

#ifndef MRW_H
#define MRW_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <armadillo>
#include "Model.h"

using namespace std;
using namespace arma;

class mrwPar
{
public:
    mat zigB;
    mat zigS;
    mat pB;
    mat pS;
    
    mrwPar() { }
    
    mat ParToMat(Par p)
    {
       mat temp(1,8);
       temp(0) = p.kON;
       temp(1) = p.kOFF;
       temp(2) = p.kONs;
       temp(3) = p.kOFFs;
       temp(4) = p.mu0;
       temp(5) = p.mu;
       temp(6) = p.muS;
       temp(7) = p.d;
       
       return temp;
    }
    
    Par MatToPar(mat m)
    {
       Par temp;
       temp.kON = m(0);
       temp.kOFF = m(1);
       temp.kONs = m(2);
       temp.kOFFs = m(3);
       temp.mu0 = m(4);
       temp.mu = m(5);
       temp.muS = m(6);
       temp.d = m(7);
       
       return temp;
    }
    
    void initPar(mat lB_m, mat lB_M, mat lS_m, mat lS_M)
    {
        pB = ((zigB==0)%pB) + ((zigB>0)%(lB_m+(randu(size(pB))%(lB_M-lB_m))));
        pS = ((zigS==0)%pB) + ((zigS>0)%(lS_m+(randu(size(pB))%(lS_M-lS_m))));
    }
    
    mat ptB()
    {
        mat pt;
        pt = pB + (randn(size(pB))%sqrt(zigB));
        return pt;
    }    
    
    mat ptS(mat ptB)
    {
        mat pt;
        pt = ((zigS==0)%ptB) + ((zigS>0)%pS) + (randn(size(pS))%sqrt(zigS));
        return pt;
    }
};

#endif /* MRW_H */

