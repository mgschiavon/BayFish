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
% DATA: Generate a random sample given a probability distribution.
%
% Created by Mariana Gómez-Schiavon
% May 2016
%
% DATA_Xrs : Generate a random sample of size n from given a probability 
%            distribution matrix M.
%
%   Xrs = DATA_Xrs(M,n,s)
%   M : Probability distribution matrix where M(i,j) is the probability of 
%       having the promoters state (i) and exactly (j-1) mRNA molecules.
%   n : Sample size, i.e. number of observations to generate.
%   s : Random number generator seed (e.g. 11)
%
%   Xr : Random sample from the given probability distribution M.
%        * NOTE: Xr has the same dimensions that M.
%
%   See also DATA_X.m

function Xrs = DATA_Xrs(M,n,s)
    P = vec2mat(M,1);	% Probability distribution vector.
    % Cumulative distribution
    Pc = P;
    for i = 2:length(P)
        Pc(i) = Pc(i) + Pc(i-1);
    end
    Pc(length(Pc)) = 1; % Correct from small numeric errors
    % Sample:
    Xrs = zeros(length(P),1);
    rng(s,'twister');   % Define random number generator.
    for i = 1:n
        r = sum([Pc<rand()])+1;
        Xrs(r) = Xrs(r) + 1;
    end
    Xrs = vec2mat(Xrs,size(M,2));
