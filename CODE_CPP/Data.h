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
 * DATA: Create data matrix 
 * 
 * Created by Mariana Gómez-Schiavon 
 * June 2016 
 * 
 * Data : Data structure.
 * 
 *  class myData
 *      Mat<int> data : Observed data matrix where Xm(i,j) is the number of 
 *          individuals in de population with exactly promoters in (i) state 
 *          and (j-1) mRNA molecules. 
 *          If N='2S', the possible promoters states are 1:{OFF,OFF}, 
 *          2:{OFF,ON}, and 3:{ON,ON}. 
 *          If N='3S', the possible promoters states are 1:{OFF,OFF},  
 *          2:{OFF,ON}, 3:{ON,ON}, 4:{OFF,ONs}, 5:{ON,ONs}, and 6:{ONs,ONs}.
 *      loadData(int N, int maxM, double a, char* myDataCode, int t) : Read 
 *          "myData_[myDataCode_t[t]_List.txt" file and load it in the data 
 *          matrix.
 *          int N : Model familiy, i.e. number of states (e.g. 2).
 *          int maxM : Maximum mRNA number to consider (e.g. 300).
 *          double a : If N='3S', threshold to define third TS state (e.g. 10).
 *          char* myDataCode : Code for specific data file (e.g. "Fos").
 *          int t : Time point to load (e.g. 5).
 * 
 */

#ifndef DATA_H
#define DATA_H

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <armadillo>

using namespace std;
using namespace arma;

class myData
{
public:
    Mat<int> data;
    
    myData(){};
    
    void loadData(int N, int maxM, double a, char* myDataCode, int t)
    {
        data.zeros(maxM+1,N+(N*(N-1)/2));
        
        char myFileName [255];
        strcpy (myFileName,"myData_");
        strcat (myFileName,myDataCode);
        strcat (myFileName,"_t%d_List.txt");
        sprintf(myFileName,myFileName,t);
        ifstream inputFile(myFileName);
        string line;
        while (getline(inputFile, line))
        {
            istringstream ss(line);
            double ts1, ts2;
            int m;
            ss >> ts1 >> ts2 >> m;
            if(ts1==0) // TS1 OFF
            {
                if(ts2==0) // TS2 OFF
                {
                    data(m,0)++;
                }
                else if(N==2 || ts2<=a) // TS2 ON
                {
                    data(m,1)++;
                }
                else // TS2 ONs
                {
                    data(m,3)++;
                }
            }
            else if(N==2 || ts1<=a) // TS1 ON
            {
                if(ts2==0) // TS2 OFF
                {
                    data(m,1)++;
                }
                else if(N==2 || ts2<=a) // TS2 ON
                {
                    data(m,2)++;
                }
                else // TS2 ONs
                {
                    data(m,4)++;
                }
            }
            else // TS1 ONs
            {
                if(ts2==0) // TS2 OFF
                {
                    data(m,3)++;
                }
                else if(N==2 || ts2<=a) // TS2 ON
                {
                    data(m,4)++;
                }
                else // TS2 ONs
                {
                    data(m,5)++;
                }
            }
        }        
    }
            
};

#endif /* DATA_H */

