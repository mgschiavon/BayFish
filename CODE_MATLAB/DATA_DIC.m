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
% DATA: Calculate the Deviance Information Criterion (DIC).
%
% Created by Mariana Gómez-Schiavon
% June 2016
%
% DATA_DIC : Calculate the Deviance Information Criterion (DIC) given the 
%            log-likelihood chain resulting from a Metropolis Random
%            Walk run (mrw).
%
%   dic = DATA_DIC(mrw,Par,ParS,myGeneModel,buT)
%   mrw : Metropolis Random Walk result.
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
%   myGeneModel : Gene model to use (e.g. 'Fos(2S,300)')
%   buT : Threshold around the maximum likelihood to define the burn in 
%         period (e.g. 1.005).
%
%   dic : Deviance Information Criterion.
%
%   See also DATA_logL.m
%   See also DATA_BIC.m
%   See also DATA_AIC.m

function dic = DATA_DIC(mrw,Par,ParS,myGeneModel,buT)
    % Model (N), maximum mRNA number to consider (maxM), and data (x):
    load(cat(2,'myData_',myGeneModel,'.mat'),'N','maxM','x');

    % Define indexes:
    myI = [];
    I1 = fieldnames(mrw.U);
    for i = 1:length(I1)
        myI = [myI;{'U',I1{i},1}];
    end
    I1 = fieldnames(mrw.S);
    psN = ''; % ParS fieldnames to name results.
    for i = 1:length(I1)
        myI = [myI;{'S',I1{i},1};{'S',I1{i},2}];
        psN = cat(2,psN,I1{i});
    end
    clear i I1
    % Specify time points (myT) and sample size (myS):
    myT = fieldnames(x);
    myS = zeros(1,length(myT));
    for t = 1:length(myT)
        myS(t) = sum(sum(x.(myT{t})));
    end
    clear t
    
    % Log-Likelihood & mean parameters:
    L = sum(mrw.L,2);
    bu = find([L>=(max(L)*buT)],1);
    L = L(bu:length(L));
    mP = mean(mrw.P(bu:length(mrw.P),:));
    for i = 1:size(myI,1)
        if(strcmp(myI{i,1},'U'))
            Par.(myI{i,2}) = mP(i);
        else
            ParS.(myI{i,2})(myI{i,3}) = mP(i);
        end
    end
    [m,l,e] = SIM_logL_Opt(x,'2S',maxM,Par,ParS);
    lT = 0;
    for i = 1:length(myT)
        lT = myS(i)*l.(myT{i});
    end
    
    Dm = sum(-2*L)/length(L);
    Dt = (-2*lT);
    dic = Dm + (Dm - Dt);
