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
% Stochastic dynamics of single-neuron early transcriptional response
% DATA: Calculate the Bayesian Information Criterion (BIC).
%
% Created by Mariana Gómez-Schiavon
% June 2016
%
% DATA_BIC : Calculate the Bayesian Information Criterion (BIC) given the 
%            maximized log-likelihood (maxL), number of free parameters 
%            (k), and sample size (n).
%
%   bic = DATA_BIC(maxL,k,n)
%   maxL : Maximized log-likelihood value.
%   k    : Number of free parameters.
%   n    : Sample size.
%
%   bic : Bayesian Information Criterion.
%
%   See also DATA_logL.m
%   See also DATA_AIC.m
%   See also DATA_DIC.m

function bic = DATA_BIC(maxL,k,n)
    bic = (-2*maxL) + (k*log(n));
