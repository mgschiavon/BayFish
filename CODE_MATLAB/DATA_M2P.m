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
% DATA: Translate a probability distribution vector to a matrix according 
%       to its model structure.
%
% Created by Mariana Gómez-Schiavon
% May 2016
%
% DATA_M2P : Translate a probability distribution vector to a matrix 
%       according to its model structure.
%
%   P = DATA_M2P(N,M)
%   M : Probability distribution matrix where M(i,j) is the probability of 
%       having the promoters state (i) and exactly (j-1) mRNA molecules.
%       If N='2S', the possible promoters states are 1:{OFF,OFF},
%       2:{OFF,ON}, and 3:{ON,ON}.
%       If N='3S', the possible promoters states are 1:{OFF,OFF}, 
%       2:{OFF,ON}, 3:{ON,ON}, 4:{OFF,ONs}, 5:{ON,ONs}, and 6:{ONs,ONs}.
%
%   P : Probability distribution vector.
%
%   See also DATA_P2M.m

function P = DATA_M2P(M)
    P = vec2mat(M,1);
