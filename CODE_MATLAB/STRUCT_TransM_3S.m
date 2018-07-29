% (C) Copyright 2018 Mariana GÃ³mez-Schiavon
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
% STRUCTURE: Create instructions for the model's transition matrix.
%
% Created by Mariana Gómez-Schiavon
% July 2018
%
% STRUCT_TransM_3S : Create instructions to construct the transition matrix 
%           under the three-state (3S) model.
%
%   [] = STRUCT_TransM_3S(maxM)
%   maxM : Maximum mRNA number to consider.
%
%   See also DATA_X.m
%   See also DATA_Xm.m

function [] = STRUCT_TransM_3S(maxM)
    %% Vector of states
    i = 1;
    for pS = 0:2
        for p = 0:(2-pS)
            for m = 0:maxM
                x(i).species = [p pS m];
                i = i + 1;
            end
        end
    end
    clear i p m

    %% Transition matrix cases
    k = zeros(7,1);	% Indexes
    for i = 1:length(x)
        for j = 1:length(x)
            nu = x(j).species - x(i).species;
			% For all i==j:
            if(nu==[0 0 0])
                k(1) = k(1)+1;
                myReaction.allNeg(k(1),:) = [i j];
			% If "promoter activation" occurs:
            elseif(nu==[1,0,0])
                k(2) = k(2)+1;
                myReaction.kON(k(2),:) = [i j];
			% If "promoter deactivation" occurs:
            elseif(nu==[-1,0,0])
                k(3) = k(3)+1;
                myReaction.kOFF(k(3),:) = [i j];
			% If "promoter super-activation" occurs:
            elseif(nu==[-1,1,0])
                k(4) = k(4)+1;
                myReaction.kONs(k(4),:) = [i j];
			% If "promoter super-deactivation" occurs:
            elseif(nu==[1,-1,0])
                k(5) = k(5)+1;
                myReaction.kOFFs(k(5),:) = [i j];
			% If "mRNA synthesis" occurs:
            elseif(nu==[0,0,1])
                k(6) = k(6)+1;
                myReaction.mus(k(6),:) = [i j];
			% If "mRNA degradation" occurs:
            elseif(nu==[0,0,-1])
                k(7) = k(7)+1;
                myReaction.d(k(7),:) = [i j];
            end
        end
    end
    clear i j k nu

    %% Functions
	%%% LEAVING [xJ] STATE %%%
    myPropensity.allNeg = @(Par,xJ) [...
        - (Par.kON*(2-xJ.species(1)-xJ.species(2)))...	% Promoter activation
        - (Par.kOFF*xJ.species(1))...                   % Promoter deactivation
        - (Par.kONs*(xJ.species(1)))...                 % Promoter super-activation
        - (Par.kOFFs*xJ.species(2))...                  % Promoter super-deactivation
        - (Par.mu0*(2-xJ.species(1)-xJ.species(2)))...	% mRNA synthesis from OFF promoters
        - (Par.mu*xJ.species(1))...                     % mRNA synthesis from ON promoters
        - (Par.muS*xJ.species(2))...                    % mRNA synthesis from ONs promoters
        - (Par.d*xJ.species(3))];                       % mRNA degradation
	%%% ENTERING [xJ] STATE %%%
    myPropensity.kON = @(Par,xJ) [Par.kON*(2-xJ.species(1)-xJ.species(2))];	% Promoter activation
    myPropensity.kOFF = @(Par,xJ) [Par.kOFF*xJ.species(1)];                 % Promoter deactivation
    myPropensity.kONs = @(Par,xJ) [Par.kONs*xJ.species(1)];                 % Promoter activation
    myPropensity.kOFFs = @(Par,xJ) [Par.kOFFs*xJ.species(2)];               % Promoter deactivation
    myPropensity.mus = @(Par,xJ) [(Par.mu0*(2-xJ.species(1)-xJ.species(2)))...% mRNA synthesis from OFF promoters
        + (Par.mu*xJ.species(1))+ (Par.muS*xJ.species(2))];                 % mRNA synthesis from ON & ONs promoters
    myPropensity.d = @(Par,xJ) [Par.d*xJ.species(3)];                       % mRNA degradation

    % Exclude cases when mRNA molecule number is in the limit of the matrix.
    % If mRNA is already the maximum, no synthesis can occur:
    myPropensity.Exc.allNeg = @(Par,xJ) [...
        - (Par.mu0*(2-xJ.species(1)-xJ.species(2)))...
        - (Par.mu*xJ.species(1))- (Par.mu*xJ.species(2))];

    %% Save
    myNames = {'allNeg','kON','kOFF','kONs','kOFFs','mus','d'};
    save(cat(2,'myTransM(3S,',num2str(maxM),').mat'),'x',...
        'myReaction','myPropensity','myNames')
