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
% DATA: Calculate the model's transition matrix.
%
% Created by Mariana Gómez-Schiavon
% May 2016
%
% DATA_TransM : Calculate the model's transition matrix given the chosen 
%               kinetic parameters.
%
%   A = DATA_TransM(N,maxM,Par)
%   N : Model ('2S' or '3S').
%   maxM : Maximum mRNA number to consider.
%   Par : Structure with the model parameters.
%      .kON : Promoter activation rate (OFF->ON)
%      .kOFF : Promoter deactivation rate (ON->OFF)
%      .kONs : Promoter "super" activation rate (ON->ONs)
%      .kOFFs : Promoter "super" deactivation rate (ONs->ON)
%      .mu0 : mRNA synthesis rate of promoter in OFF state
%      .mu : mRNA synthesis rate of promoter in ON state
%      .muS : mRNA synthesis rate of promoter in ONs state
%      .d : mRNA degradation rate
%
%   See also STRUCT_TransM.m
%   See also DATA_Pss.m

function A = DATA_TransM(N,maxM,Par)
    load(cat(2,'myTransM(',N,',',num2str(maxM),').mat'))
    % Transition matrix
    A = zeros(length(x),length(x));
    for myI = 1:length(myNames)
        myR = myNames{myI};
        for myJ = 1:length(myReaction.(myR))
            i = myReaction.(myR)(myJ,1);
            j = myReaction.(myR)(myJ,2);
            A(j,i) = myPropensity.(myR)(Par,x(i));
        end
    end
end
