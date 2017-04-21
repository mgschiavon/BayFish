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
% DATA: Find unique parameter sets in MRW replicas.
%
% Created by Mariana Gómez-Schiavon
% May 2016
%
% DATA_UniquePar : Find unique parameter sets in MRW replicas.
%
%   [uP,iP] = DATA_UniquePar(myMRW,myS,buT,pT)
%   myMRW : Results files common name (e.g. 'MRW_Fos(2S,300)(kON)')
%   myS : Replica numbers (e.g. [1:3])
%   buT : Threshold around the maximum likelihood to define the burn in 
%         period (e.g. 1.005)
%   pT : Parameter resolution to define unique parameter sets (e.g. 0.01)
%   
%   uP : Unique parameter sets.
%   iP : Occurence of each unique parameter set.
%
%   See also SIM_mrw.m

function [uP,iP] = DATA_UniquePar(myMRW,myS,buT,pT)
    % Load results:
    P = [];
    for s = myS
        load(cat(2,myMRW,'(s',num2str(s),').mat'))
        bu = find([sum(mrw.L,2)>=(max(sum(mrw.L,2))*buT)],1);
        P = [P;mrw.P(bu:mrw.I,:)];
    end
    clear s bu

    %% Unique parameters
    mP = mean(P,1);
    myE = zeros(1,length(mP))+1;
    r = floor((pT*mP)./(10.^myE));
    while(sum([r==0]))
        myE([r==0]) = myE([r==0]) - 1;
        r = floor((pT*mP)./(10.^myE));
    end
    clear r
    
    for i = 1:size(P,2)
        P(:,i) = round(P(:,i)/(10^myE(i)))*(10^myE(i));
    end
    clear i
    
    [uP,i,iP] = unique(P,'rows');
    iP = hist(iP,[1:length(uP)]);
    clear i
    
