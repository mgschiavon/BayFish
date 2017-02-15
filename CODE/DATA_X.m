% (C) Copyright 2017 Mariana GÃ³mez-Schiavon
%
%    This file is part of BayFish.
%
%    BayFish is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    BayFish is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with BayFish.  If not, see <http://www.gnu.org/licenses/>.
%
% BayFish pipeline
% DATA: Load data from excel file
%
% Created by Mariana Gómez-Schiavon
% May 2016
%
% DATA_X : Load data from an excel file and correctly format it.
%
%   X = DATA_X(FileName,SheetName,NaNtoZero,ExpInfo)
%   FileName : File name (e.g., 'DataMerged.xlsx').
%   SheetName : Sheet name (e.g., 'B','I', or 'E').
%   NaNtoZero : If 1, substitute NaN data points to zero.
%   ExpInfo : If 1, extract experiment information (fourth column).
%
%   See also DATA_Xm.m

function X = DATA_X(FileName,SheetName,NaNtoZero,ExpInfo)
    [X.data,y,x] = xlsread(FileName,SheetName);
	X.Name = x{1,1};
	X.Columns = x(2,1:3);
    if(ExpInfo==1)
        X.ExpList = x(3:length(x),4);
    end
    clear x y
    
    if(NaNtoZero==1)
        X.data(isnan(X.data)) = 0;
    end
