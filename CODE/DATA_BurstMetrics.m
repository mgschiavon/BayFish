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
% DATA: Find burst metrics in MRW replicas.
%
% Created by Mariana Gómez-Schiavon
% November 2016
%
% DATA_BurstMetrics : Find burst metrics in MRW replicas.
%
%   [BM, PM] = DATA_BurstMetrics(myMRW,myS,buT)
%   myMRW : Results files common name (e.g. 'MRW_Fos(2S,300)(kON)')
%   myS : Replica numbers (e.g. [1:3])
%   buT : Threshold around the maximum likelihood to define the burn-out 
%         period (e.g. 1.005)
%   
%   BM : Burst metrics structure.
%       .fB   : Fraction of ON promoters
%       .tON  : OFF state duration
%       .tOFF : Burst (ON state) duration
%       .amp  : Burst amplitude
%   PM : Parameter metrics structure.
%       .mean : Vector of mean values
%       .std  : Vector of standard deviation values
%
%   See also SIM_mrw.m
%   See also DATA_UniquePar.m

function [BM, PM] = DATA_BurstMetrics(myMRW,myS,buT)
    % Load results:
    P = [];
    for s = myS
        load(cat(2,myMRW,'(s',num2str(s),').mat'),'mrw','Par')
        bu = find([sum(mrw.L,2)>=(max(sum(mrw.L,2))*buT)],1);
        P = [P;mrw.P(bu:mrw.I,:)];
    end
    clear s bu

    %% Parameter metrics:
    PM.mean = mean(P);
    PM.std = std(P);
    
    %% Define indexes:
    p.U = Par;
    p.S = Par;
    I1 = fieldnames(mrw.U);
    j = 1;
    for i = 1:length(I1)
        p.U.(I1{i}) = P(:,j);
        p.S.(I1{i}) = P(:,j);
        j = j + 1;
    end
    I1 = fieldnames(mrw.S);
    for i = 1:length(I1)
        p.U.(I1{i}) = P(:,j);
        p.S.(I1{i}) = P(:,j+1);
        j = j + 2;
    end
    clear i j I1
    
    %% Burst metrics:
    c = {'U','S'};
    for i = 1:2
        BM.(c{i}).fB = [mean(p.(c{i}).kON./(p.(c{i}).kON+p.(c{i}).kOFF)) ...
            std(p.(c{i}).kON./(p.(c{i}).kON+p.(c{i}).kOFF))];
        BM.(c{i}).tON = [mean(1./p.(c{i}).kOFF) ...
            std(1./p.(c{i}).kOFF)];
        BM.(c{i}).tOFF = [mean(1./p.(c{i}).kON) ...
            std(1./p.(c{i}).kON)];
        BM.(c{i}).amp = [mean((p.(c{i}).mu./p.(c{i}).kOFF)...
            .*(1-exp(-p.(c{i}).d./p.(c{i}).kOFF))) ...
            std((p.(c{i}).mu./p.(c{i}).kOFF)...
            .*(1-exp(-p.(c{i}).d./p.(c{i}).kOFF)))];
    end
    clear c i
    