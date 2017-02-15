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
% Stochastic dynamics of single-neuron early transcriptional response
% FIGURE: Evaluate & plot MRW for synthetic random samples.
%
% Created by Mariana Gómez-Schiavon
% August 2016
%
% FIG_EvalRS : Evaluate & plot MRW for synthetic random samples efficiency 
%              recapitulating the "true" model parameters.
%
%   [] = FIG_EvalRS(myModel,myRS,OptP,myR,myS,buT,Hb,SaveFigs)
%   myModel : Code for MRW model to use as reference (e.g. 
%             'Fos(2S,300)(kON)').
%   myRS : Code for Random Sample MRW to use (e.g.
%          'RS(2S,300,kON)(n1000s#)(kON)').
%   OptP : Set of "real" parameters to evaluate MRW performance (i.e. 
%          [0.0459,0.026,5.65,0.0049,0.111]).
%   myR : Replica numbers for random samples (e.g. [1:3])
%   myS : Replica numbers for MRW runs (e.g. [1:3])
%   buT : Threshold around the maximum likelihood to define the burn-out 
%         period (e.g. 1.005)
%   Hb : Number of bins for 2D-histogram plot (e.g. 25)
%   SaveFigs : If 1, save figures as MRW_[myModel]_[myRS]_EvalRS.png.
%
%   See also FIG_CompRS.m
%   See also FIG_CompGenes.m
%   See also FIG_CompModels.m

function [] = FIG_EvalRS(myModel,myRS,OptP,myR,myS,buT,Hb,SaveFigs)
    myC = [0.6 0.8 0;
        1 0.6 0;
        0 0.8 0.8;
        0.8 0.6 1;
        1 0 0.4;
        0.08 0.17 0.55];
    
    %% Load main model:
    for s = myS
        load(cat(2,'MRW_',myModel,'(s',num2str(s),').mat'),'mrw')
        bu = find([sum(mrw.L,2)>=(max(sum(mrw.L,2))*buT)],1);
        Data(1,s).P = mrw.P(bu:mrw.I,:);
        clear bu
    end
    % Define parameter labels:
    myPn = [];
    I1 = fieldnames(mrw.U);
    for i = 1:length(I1)
        myPn = [myPn;{I1{i}}];
    end
    I1 = fieldnames(mrw.S);
    for i = 1:length(I1)
        myPn = [myPn;{cat(2,I1{i},'^{(B)}')};{cat(2,I1{i},'^{(S)}')}];
    end
    myPn = strrep(myPn,'ON','_{ON}');
    myPn = strrep(myPn,'OFF','_{OFF}');
    myPn = strrep(myPn,'mu','\mu');
    myPn = strrep(myPn,'0','_0');
    clear i I1 mrw
    
    %% Load random samples:
    for r = myR
        for s = myS
            load(cat(2,'MRW_',strrep(myRS,'#',num2str(r)),'(s',...
                num2str(s),')','.mat'),'mrw')
            bu = find([sum(mrw.L,2)>=(max(sum(mrw.L,2))*buT)],1);
            Data(r+1,s).P = mrw.P(bu:mrw.I,:);
        end
    end
    clear r s mrw bu
    
    %% Parameters
    ssP = zeros(size(Data(1).P,2),length(Data));
    for l = 1:size(Data(1).P,2)
        p = [];
        for s = myS
            p = [p;Data(1,s).P(:,l)];
        end
        ssP(l,1) = sum(((p-OptP(l))/OptP(l)).^2)/length(p);
        
        for r = myR
            p = [];
            for s = myS
                p = [p;Data(r+1,s).P(:,l)];
            end
            ssP(l,r+1) = sum(((p-OptP(l))/OptP(l)).^2)/length(p);
        end
    end
    clear l s
    
    %% Figure:
    fig = figure();
    fig.Units = 'inches';
    fig.PaperPosition = [1 1 6 4];
    fig.Position = fig.PaperPosition;
    % Title:
    annotation('textbox',[0.15 0.96 0.8 0.04],...
        'String',{cat(2,'MRW :',myModel,' vs ',myRS)},...
        'LineStyle','none','HorizontalAlignment','center',...
        'FontWeight','bold','FontSize',12);

    hold on;
    plot(ssP(:,1),'MarkerSize',8,'Marker','o','linestyle',':')
    plot(ssP(:,myR+1),'MarkerSize',15,'Marker','.','linestyle',':')
        xlim([0.5 5.5])
        ylabel('\langle(\rho-\rho_*)^2/\rho_*\rangle')
        legend({'Fos','RS_1','RS_2','RS_3'})
        set(gca,'YGrid','on','XTick',[1:length(myPn)],...
            'XTickLabel',myPn)

    if(SaveFigs==1)
        print(fig,cat(2,'MRW_',myModel,'_',myRS,'_EvalRS.png'),'-dpng','-r300')
    end
