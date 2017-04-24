/*
 * (C) Copyright 2017 Mariana Gómez-Schiavon
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
 * MODEL: Define underlying model structure
 * 
 * Created by Mariana Gómez-Schiavon 
 * June 2016 
 * 
 * Model : Model structure and transition matrix.
 * 
 *  class Par : Structure with the kinetic parameters
 *      double kON   : Promoter activation rate (OFF->ON) 
 *      double kOFF  : Promoter deactivation rate (ON->OFF) 
 *      double kONs  : Promoter "super" activation rate (ON->ONs) 
 *      double kOFFs : Promoter "super" deactivation rate (ONs->ON) 
 *      double mu0   : mRNA synthesis rate of promoter in OFF state 
 *      double mu    : mRNA synthesis rate of promoter in ON state 
 *      double muS   : mRNA synthesis rate of promoter in ONs state 
 *      double d     : mRNA degradation rate 
 * 
 *  class ModelStruct(int myN, int myMaxM) : Create instructions to construct 
 *      the transition matrix under the given model (N = myN) and maximum 
 *      number of mRNA molecules (maxM = myMaxM).
 *      Mat<int> S : Lists of possible system states (species), i.e. each row 
 *          is [#ON promoters, #mRNA molecules] if N=2, or 
 *          [#ON promoters, #ONs promoters, #mRNA molecules] if N=3.
 *      mat R : Lists the positions of the transition matrix with 
 *          propensities different from zero and the kind of reaction occuring.
 *      int N : Model family, i.e. number of states (e.g. 2).
 *      int maxM : Maximum mRNA number to consider (e.g. 300).
 * 
 *      mat TransM(Par p) : Given the input parameters and previously defined N 
 *          and maxM, returns the transition matrix for the model.
 * 
 */

#ifndef MODEL_H
#define MODEL_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <armadillo>

using namespace std;
using namespace arma;

class Par
{
public:
    double kON, kOFF;
    double kONs, kOFFs;
    double mu0, mu, muS;
    double d;
    
    Par()
    {
        kON = 0;
        kOFF = 0;
        kONs = 0;
        kOFFs = 0;
        mu0 = 0;
        mu = 0;
        muS = 0;
        d = 0;
    }
};

class ModelStruct
{
public:
    Mat<int> S;
    mat R;
    int N;
    int maxM;
    
    ModelStruct(int myN, int myMaxM)
    {
        N = myN;
        maxM = myMaxM;
        if(N==2)
        {
            species2S();
            reactions2S();
        }
        else if(N==3)
        {
            species3S();
            reactions3S();
        }
        else
        {
            cout << "ERROR: Model non defined." << endl;
        }
    }
    
    mat TransM(Par p)
    {
        mat A;
        if(N==2)
        {
            A = transM_2S(p.kON, p.kOFF, p.mu0, p.mu, p.d);
        }
        else if(N==3)
        {
            A = transM_3S(p.kON, p.kOFF, p.kONs, p.kOFFs, p.mu0, p.mu, p.muS, p.d);
        }
        return A;
    }
    
    void species2S()
    {
        S.set_size(3*(maxM+1),2);
        int i = 0;
        for(int p = 0; p <= 2; p++)
        {
            for(int m = 0; m <= maxM; m++)
            {
                S(i,0) = p;
                S(i,1) = m;
                i++;
            }
        }
    }
    
    void reactions2S()
    {
        R.set_size((5*3*(maxM+1))-(2*(3+(maxM+1))),3);
        Mat<int> nu(1,2);
        int ij = 0;
        int r;
        for(int i = 0; i < S.n_rows; i++)
        {
            for(int j = 0; j < S.n_rows; j++)
            {
                r = -1;
                nu = S.row(j) - S.row(i);
                if(nu(0,1)==0)
                {
                    if(nu(0,0)==0)
                    {
                        r = 0;
                    }
                    else if(nu(0,0)==1)
                    {
                        r = 1;
                    }
                    else if(nu(0,0)==-1)
                    {
                        r = 2;
                    }
                }
                else if(nu(0,0)==0)
                {
                    if(nu(0,1)==1)
                    {
                        r = 3;
                    }
                    else if(nu(0,1)==-1)
                    {
                        r = 4;
                    }
                }
                if(r >= 0)
                {
                    R(ij,0) = r;
                    R(ij,1) = j;
                    R(ij,2) = i;
                    ij++;
                }
            }
        }
    }

    mat transM_2S(double kON, double kOFF, double mu0, double mu, double d)
    {
        mat A(S.n_rows,S.n_rows,fill::zeros);
        int row,col;
        for(int i = 0; i < R.n_rows; i++)
        {
            row = R(i,1);
            col = R(i,2);
            if(R(i,0)==0)       // Diagonal, i.e. all negative
            {
                A(row,col) = -(kON*(2-S(col,0)))  // Promoter activation
                        -(kOFF*S(col,0))      // Promoter deactivation
                        -(mu0*(2-S(col,0)))   // mRNA synthesis from OFF promoters
                        -(mu*S(col,0))        // mRNA synthesis from ON promoters
                        -(d*S(col,1));        // mRNA degradation
            }
            else if(R(i,0)==1)   // Promoter activation
            {
                A(row,col) = kON*(2-S(col,0));
            }
            else if(R(i,0)==2)   // Promoter deactivation
            {
                A(row,col) = kOFF*S(col,0);
            }
            else if(R(i,0)==3)   // mRNA synthesis
            {
                A(row,col) = mu0*(2-S(col,0)) + mu*S(col,0);
            }
            else if(R(i,0)=4)   // mRNA degradation
            {
                A(row,col) = d*S(col,1);
            }
        }
        return A;
    }
    
    void species3S()
    {
        S.set_size(6*(maxM+1),3);
        int i = 0;
        for(int pS = 0; pS <= 2; pS++)
        {
            for(int p = 0; p <= (2-pS); p++)
            {
                for(int m = 0; m <= maxM; m++)
                {
                    S(i,0) = p;
                    S(i,1) = pS;
                    S(i,2) = m;
                    i++;
                }
            }
        }
    }
    
    void reactions3S()
    {
        R.set_size((6*(maxM+1))+(2*6*maxM)+(12*(maxM+1)),3);
        Mat<int> nu(1,3);
        int ij = 0;
        int r;
        for(int i = 0; i < S.n_rows; i++)
        {
            for(int j = 0; j < S.n_rows; j++)
            {
                r = -1;
                nu = S.row(j) - S.row(i);
                if(nu(0,0)==0 && nu(0,1)==0 && nu(0,2)==0)      // Diagonal
                {
                    r = 0;
                }
                else if(nu(0,0)==1 && nu(0,1)==0 && nu(0,2)==0) // kON
                {
                    r = 1;
                }
                else if(nu(0,0)==-1 && nu(0,1)==0 && nu(0,2)==0) // kOFF
                {
                    r = 2;
                }
                else if(nu(0,0)==-1 && nu(0,1)==1 && nu(0,2)==0) // kONs
                {
                    r = 3;
                }
                else if(nu(0,0)==1 && nu(0,1)==-1 && nu(0,2)==0) // kOFFs
                {
                    r = 4;
                }
                else if(nu(0,0)==0 && nu(0,1)==0 && nu(0,2)==1) // mus
                {
                    r = 5;
                }
                else if(nu(0,0)==0 && nu(0,1)==0 && nu(0,2)==-1) // d
                {
                    r = 6;
                }
                if(r >= 0)
                {
                    R(ij,0) = r;
                    R(ij,1) = j;
                    R(ij,2) = i;
                    ij++;
                }
            }
        }
    }
    
    mat transM_3S(double kON, double kOFF, double kONs, double kOFFs, double mu0, double mu, double muS, double d)
    {
        mat A(S.n_rows,S.n_rows,fill::zeros);
        int row,col;
        for(int i = 0; i < R.n_rows; i++)
        {
            row = R(i,1);
            col = R(i,2);
            if(R(i,0)==0)       // Diagonal, i.e. all negative
            {
                A(row,col) = -(kON*(2-S(col,0)-S(col,1))) // Promoter activation (OFF->ON)
                        -(kOFF*S(col,0))      // Promoter deactivation (ON->OFF)
                        -(kONs*S(col,0))      // Promoter super-activation (ON->ONs)
                        -(kOFF*S(col,1))      // Promoter super-deactivation (ONs->ON)
                        -(mu0*(2-S(col,0)-S(col,1))) // mRNA synthesis from OFF promoters
                        -(mu*S(col,0))        // mRNA synthesis from ON promoters
                        -(muS*S(col,1))       // mRNA synthesis from ONs promoters
                        -(d*S(col,2));        // mRNA degradation
            }
            else if(R(i,0)==1)   // Promoter activation (OFF->ON)
            {
                A(row,col) = kON*(2-S(col,0)-S(col,1));
            }
            else if(R(i,0)==2)   // Promoter deactivation (ON->OFF)
            {
                A(row,col) = kOFF*S(col,0);
            }
            else if(R(i,0)==3)   // Promoter super-activation (ON->ONs)
            {
                A(row,col) = kONs*S(col,0);
            }
            else if(R(i,0)==4)   // Promoter super-deactivation (ONs->ON)
            {
                A(row,col) = kOFFs*S(col,1);
            }
            else if(R(i,0)==5)   // mRNA synthesis
            {
                A(row,col) = mu0*(2-S(col,0)-S(col,1)) + mu*S(col,0) 
                        + muS*S(col,1);
            }
            else if(R(i,0)==6)   // mRNA degradation
            {
                A(row,col) = d*S(col,2);
            }
        }
        return A;
    }
};

#endif /* MODEL_H */

