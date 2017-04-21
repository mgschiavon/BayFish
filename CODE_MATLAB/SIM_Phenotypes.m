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
% SIMULATE: Calculate phenotype (i.e. probability distributions) for the given 
%       model for each unique parameter set in the MRW runs.
%
% Created by Mariana Gómez-Schiavon
% May 2016
%
% SIM_Phenotypes : Calculate phenotype (i.e. probability distributions) 
%                  for the given model for each unique parameter set in 
%                  the MRW runs.
%
%   [] = SIM_Phenotypes(myMRW,myS,buT,pT)
%   myMRW : Results files common name (e.g. 'MRW_Fos(2S,300)(kON)')
%   myS : Replica numbers (e.g. [1:3])
%   buT : Threshold around the maximum likelihood to define the burn in 
%         period (e.g. 1.005)
%   pT : Parameter resolution to define unique parameter sets (e.g. 0.01)
%   
%   OUTPUT : File [myMRW]_Ms.mat.
%       uP : Unique parameter sets.
%       iP : Occurence of each unique parameter set.
%       sL : SIM_logL_Opt results per unique parameter set.
%
%   See also SIM_mrw.m
%   See also DATA_UniquePar.m
%   See also SIM_logL_Opt.m

function [] = SIM_Phenotypes(myMRW,myS,buT,pT)
    % Load observed data:
    myX = regexp(strrep(myMRW,'_',' '),'\w*\(\dS,\d*\)','match');
    load(cat(2,'myData_',myX{1},'.mat'));
    clear myX
    % Define unique parameter sets:
    [uP,iP] = DATA_UniquePar(myMRW,myS,buT,pT);
    % Calculate phenotype (probability distributions) for each unique set:
    %%% Define indexes:
    load(cat(2,myMRW,'(s',num2str(myS(1)),').mat'));
    myI = [];
    I1 = fieldnames(mrw.U);
    for i = 1:length(I1)
        myI = [myI;{'U',I1{i},1}];
    end
    I1 = fieldnames(mrw.S);
    for i = 1:length(I1)
        myI = [myI;{'S',I1{i},1};{'S',I1{i},2}];
    end
    clear i I1 mrw
    %%% Iterate over unique sets:
    for i = 1:length(uP)
        % Kinetic parameters:
        for ii = 1:size(myI,1)
            if(strcmp(myI{ii,1},'U'))
                Par.(myI{ii,2}) = uP(i,ii);
            else
                ParS.(myI{ii,2})(myI{ii,3}) = uP(i,ii);
            end
        end
        [sL(i).M,sL(i).L,sL(i).eps] = SIM_logL_Opt(x,N,maxM,Par,ParS);
        % Save progress:
        if(mod(i,round(length(uP)/50))==0)
            cat(2,num2str(round(100*i/length(uP))),'%')
            save(cat(2,'TEMP_',myMRW,'_Ms.mat'),'i','uP','iP','sL')
        end
    end
    clear i ans
    %%% Save results:
    save(cat(2,myMRW,'_Ms.mat'),'uP','iP','sL')
    delete(cat(2,'TEMP_',myMRW,'_Ms.mat'))
