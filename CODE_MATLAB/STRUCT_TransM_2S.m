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
% STRUCTURE: Create instructions for the model's transition matrix.
%
% Created by Mariana Gómez-Schiavon
% May 2016
%
% STRUCT_TransM_2S : Create instructions to construct the transition matrix 
%           under the two-state (2S) model.
%
%   [] = STRUCT_TransM_2S(maxM)
%   maxM : Maximum mRNA number to consider.
%
%   See also DATA_X.m
%   See also DATA_Xm.m

function [] = STRUCT_TransM_2S(maxM)
    %% Vector of states
    i = 1;
    for p = 0:2
        for m = 0:maxM
            x(i).species = [p m];
            i = i + 1;
        end
    end
    clear i p m

    %% Transition matrix cases
    k = zeros(5,1);	% Indexes
    for i = 1:length(x)
        for j = 1:length(x)
            nu = x(j).species - x(i).species;
			% For all i==j:
            if(nu==[0,0])
                k(1) = k(1)+1;
                myReaction.allNeg(k(1),:) = [i j];
			% If "promoter activation" occurs:
            elseif(nu==[1,0])
                k(2) = k(2)+1;
                myReaction.kON(k(2),:) = [i j];
			% If "promoter deactivation" occurs:
            elseif(nu==[-1,0])
                k(3) = k(3)+1;
                myReaction.kOFF(k(3),:) = [i j];
			% If "mRNA synthesis" occurs:
            elseif(nu==[0,1])
                k(4) = k(4)+1;
                myReaction.mus(k(4),:) = [i j];
			% If "mRNA degradation" occurs:
            elseif(nu==[0,-1])
                k(5) = k(5)+1;
                myReaction.d(k(5),:) = [i j];
            end
        end
    end
    clear i j k nu

    %% Functions
	%%% LEAVING [xJ] STATE %%%
    myPropensity.allNeg = @(Par,xJ) [...
        - (Par.kON*(2-xJ.species(1)))...	% Promoter deactivation
        - (Par.kOFF*xJ.species(1))...		% Promoter activation
        - (Par.mu0*(2-xJ.species(1)))...	% mRNA synthesis from OFF promoters
        - (Par.mu*xJ.species(1))...			% mRNA synthesis from ON promoters
        - (Par.d*xJ.species(2))];			% mRNA degradation
	%%% ENTERING [xJ] STATE %%%
    myPropensity.kON = @(Par,xJ) [Par.kON*(2-xJ.species(1))];	% Promoter activation
    myPropensity.kOFF = @(Par,xJ) [Par.kOFF*xJ.species(1)];		% Promoter deactivation
    myPropensity.mus = @(Par,xJ) [(Par.mu0*(2-xJ.species(1)))...% mRNA synthesis from OFF promoters
        + (Par.mu*xJ.species(1))];								% mRNA synthesis from ON promoters
    myPropensity.d = @(Par,xJ) [Par.d*xJ.species(2)];			% mRNA degradation

    % Exclude cases when mRNA molecule number is in the limit of the matrix.
    % If mRNA is already the maximum, no synthesis can occur:
    myPropensity.Exc.allNeg = @(Par,xJ) [...
        - (Par.mu0*(2-xJ.species(1)))...
        - (Par.mu*xJ.species(1))];

    %% Save
    myNames = {'allNeg','kON','kOFF','mus','d'};
    save(cat(2,'myTransM(2S,',num2str(maxM),').mat'),'x',...
        'myReaction','myPropensity','myNames')
