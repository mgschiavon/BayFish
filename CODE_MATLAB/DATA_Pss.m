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
% DATA: Calculate the stationary distribution vector given a 
%            transition matrix.
%
% Created by Mariana Gómez-Schiavon
% May 2016
%
% DATA_Pss : Calculate the stationary distribution vector given a 
%            transition matrix.
%
%   Pss = DATA_Pss(A)
%   A : Transition matrix describing the model.
%   Pss : Stationary distribution with states as described in A structure
%       (e.g. 'myTransM(2S,300).mat').
%
%   See also DATA_TransM.m
%   See also DATA_Pxt.m

function Pss = DATA_Pss(A)
    [V,lambda] = eigs(A,1,0);
    Pss = abs(V/sum(V));
end
