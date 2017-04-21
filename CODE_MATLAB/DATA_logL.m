% (C) Copyright 2017 Mariana Gómez-Schiavon
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
% DATA: Calculate the log-likelihood of observing the data x given the
%       probability distribution M.
%
% Created by Mariana G�mez-Schiavon
% May 2016
%
% DATA_logL : Calculate the log-likelihood of observing the data x given the
%             probability distribution M.
%
%   L = DATA_logL(M,x)
%   M : Probability matrix where M(i,j) is the probability of 
%       having the promoters state (i) and exactly (j-1) mRNA molecules.
%   x : Observed data matrix where x(i,j) is the number of 
%       individuals in de population with promoters in i-state and exactly 
%       (j-1) mRNA molecules.
%     * NOTE: M and x must have the same dimensions.
%
%   L : log-likelihood of x given M, normalized by number of observations.
%
%   See also DATA_P2M.m
%   See also DATA_Xm.m

function L = DATA_logL(M,x)
    % Imposing a minimal value for the probability values to avoid
    % numerical errors:
    M([M<(1e-100)]) = 1e-100;
    L = x.*log(M);
    L = sum(sum(L))/sum(sum(x));
