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
% SIMULATE: Calculate the log-likelihood of observing the data x  for a 
%           given model and set of kinetic parameters.
%
% Created by Mariana Gómez-Schiavon
% May 2016
%
% SIM_logL : Starting from a model structure, and its kinetic parameters, 
%            calculate the log-likelihood of observing the data x.
%
%   [M,L,eps] = SIM_logL(x,N,maxM,Par,ParS)
%   x : Observed data matrix where x(i,j) is the number of 
%       individuals in de population with promoters in i-state and exactly 
%       (j-1) mRNA molecules.
%   N : Model ('2S' or '3S').
%   maxM : Maximum mRNA number to consider.
%   Par : Structure with the kinetic parameters.
%      .kON : Promoter activation rate (OFF->ON)
%      .kOFF : Promoter deactivation rate (ON->OFF)
%      .kONs : Promoter "super" activation rate (ON->ONs)
%      .kOFFs : Promoter "super" deactivation rate (ONs->ON)
%      .mu0 : mRNA synthesis rate of promoter in OFF state
%      .mu : mRNA synthesis rate of promoter in ON state
%      .muS : mRNA synthesis rate of promoter in ONs state
%      .d : mRNA degradation rate
%   ParS : Structure with the kinetic parameters sensitive to stimulus as 
%      elements, and a 2-vector with the value of the given parameter 
%      in basal and stimulus conditions, [Basal,Stimulus] (e.g. 
%      ParS.kON = [0.0051,0.1373]).
%
%   M : Probability matrix matrix where M(i,j) is the probability of 
%       having the promoters state (i) and exactly (j-1) mRNA molecules.
%       If N='2S', the possible promoters states are 1:{OFF,OFF},
%       2:{OFF,ON}, and 3:{ON,ON}.
%       If N='3S', the possible promoters states are 1:{OFF,OFF}, 
%       2:{OFF,ON}, 3:{ON,ON}, 4:{OFF,ONs}, 5:{ON,ONs}, and 6:{ONs,ONs}.
%   L : log-likelihood of x given M, normalized by number of observations.
%   eps : FSP algorithm estimation error.
%
%   See also DATA_Xm.m
%   See also DATA_TransM.m
%   See also DATA_Pss.m
%   See also DATA_Pxt.m
%   See also DATA_P2M.m
%   See also DATA_logL.m

function [M,L,eps] = SIM_logL(x,N,maxM,Par,ParS)
    pS = fieldnames(ParS);  % Parameters sensitive to stimulus
    myT = fieldnames(x);    % Time points
    % Basal conditions:
        for p = 1:length(pS)
            Par.(pS{p}) = ParS.(pS{p})(1);
        end
        A = DATA_TransM(N,maxM,Par);
        % Probability distribution vector:
        P00 = DATA_Pss(A);
        % Probability distribution matrix:
        M.t00 = DATA_P2M(N,P00);
        % Calculate the log-likelihood:
        L.t00 = DATA_logL(M.t00,x.t00);
        % Calculate estimation error:
        eps.t00 = 1-sum(abs(P00));
    % Stimulus conditions:
        for p = 1:length(pS)
            Par.(pS{p}) = ParS.(pS{p})(2);
        end
        A = DATA_TransM(N,maxM,Par);
        for t = 2:length(myT)
            tt = regexp(myT{t},'\.*\d*','match');
            % Probability distribution vector:
            P = DATA_Pxt(A,P00,str2double(tt{1}));
            % Probability distribution matrix:
            M.(myT{t}) = DATA_P2M(N,P);
            % Calculate the log-likelihood:
            L.(myT{t}) = DATA_logL(M.(myT{t}),x.(myT{t}));
            % Calculate FSP algorithm estimation error:
            eps.(myT{t}) = 1-sum(P);
            if(eps.(myT{t})>0.001)
                cat(2,'CAUTION: Large FSP algorithm estimation error (',...
                    num2str(eps.(myT{t})),')')
            end
        end
    clear pS p myT t tt A P
