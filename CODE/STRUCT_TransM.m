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
% STRUCTURE: Call the appropiate function to create instructions for the 
%            transition matrix according to the chosen model.
%
% Created by Mariana Gómez-Schiavon
% May 2016
%
% STRUCT_TransM : Call the appropiate function to create instructions for  
%            the transition matrix according to the chosen model.
%
%   [] = STRUCT_TransM(N,maxM)
%   N : Model ('2S' or '3S').
%   maxM : Maximum mRNA number to consider.
%
%   See also STRUCT_TransM_2S.m
%   See also STRUCT_TransM_3S.m

function [] = STRUCT_TransM(N,maxM)
    if(strcmp(N,'2S'))
        STRUCT_TransM_2S(maxM);
    elseif(strcmp(N,'3S'))
        STRUCT_TransM_3S(maxM)
    else
        cat(2,'ERROR: Model ',N,' has not been defined.')
    end
end
