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
% SIMULATE: Generate synthetic data.
%
% Created by Mariana Gómez-Schiavon
% May 2016
%
% SIM_srs : Generate synthetic data using the chosen model and sample size.
%
%   [] = SIM_srs(N,maxM,srs,Par,ParS)
%   N : Model ('2S' or '3S').
%   maxM : Maximum mRNA number to consider.
%   srs : Structure with the Synthetic Random Sampling parameters:
%      .t : Time points (e.g. {'t00','t05','t15','t25'}).
%      .n : Sample size (e.g. [173,174,122,115]).
%      .s : Random number generator seed (e.g. 1).
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
%   OUTPUT FILE : 'myData_RS([N],[maxM],[fieldnames(ParS)])(n[srs.n(1)]s[srs.s]).mat'
%
%   See also DATA_Xm.m
%   See also SIM_logL.m
%   See also DATA_Xrs.m

function [] = SIM_srs(N,maxM,srs,Par,ParS)
    % Define probability distribution to sample from:
    for t = 1:length(srs.t)
        x.(srs.t{t}) = DATA_Xm(zeros(2,3),N,maxM,0)*0;
    end
    [M,l,e] = SIM_logL(x,N,maxM,Par,ParS);
    clear t l e

    % Generate synthetic data:
    for t = 1:length(srs.t)
        x.(srs.t{t}) = DATA_Xrs(M.(srs.t{t}),srs.n(t),(srs.s*10)+t);
    end
    clear t
    
    % Save:
    ps = fieldnames(ParS);
    psN = '';
    for i = 1:length(ps)
        psN = cat(2,psN,ps{i});
    end
    clear ps i
    save(cat(2,'myData_RS(',N,',',num2str(maxM),',',psN,')',...
        '(n',num2str(srs.n(1)),'s',num2str(srs.s),').mat'),...
        'N','maxM','srs','Par','ParS','x')