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
% DATA: Calculate the Akaike Information Criterion (AIC).
%
% Created by Mariana Gómez-Schiavon
% June 2016
%
% DATA_AIC : Calculate the Akaike Information Criterion (AIC) given the 
%            maximized log-likelihood (maxL), and number of free 
%            parameters (k); if you want the correction by finite sample 
%            size (n), c = 1.
%
%   aic = DATA_AIC(maxL,k,n,c)
%   maxL : Maximized log-likelihood value.
%   k    : Number of free parameters.
%   n    : Sample size.
%   c    : If 1, apply correction by sample size.
%
%   aic : Akaike Information Criterion.
%
%   See also DATA_logL.m
%   See also DATA_BIC.m
%   See also DATA_DIC.m

function aic = DATA_AIC(maxL,k,n,c)
    aic = (-2*maxL) + (2*k);
    if(c==1)
        aic = aic + ((2*k*(k+1))/(n-k-1));
    end
