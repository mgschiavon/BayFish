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
% DATA: Calculate the probability distribution at time t given a 
%       transition matrix and some initial condition.
%
% Created by Mariana Gómez-Schiavon
% May 2016
%
% DATA_Pxt : Calculate the probability distribution at time t given a 
%            transition matrix and some initial condition.
%
%   Pt = DATA_Pxt(A,P0,t)
%   A : Transition matrix describing the model.
%   P0 : Initial probability distribution with states as described in A 
%        structure (e.g. 'myTransM(2S,300).mat').
%   t : Time (e.g. 5).
%   Pt : Probability distribution at time t with states as described in A 
%        structure (e.g. 'myTransM(2S,300).mat').
%
%   See also DATA_Pss.m
%   See also DATA_P2M.m

function Pt = DATA_Pxt(A,P0,t)
    Pt = expm(A*t)*P0;
