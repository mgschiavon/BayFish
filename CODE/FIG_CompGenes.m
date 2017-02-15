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
% FIGURE: Plot and compare MRW for different genes or constructs.
%
% Created by Mariana Gómez-Schiavon
% August 2016
%
% FIG_CompGenes : Plot and compare MRW for different genes or constructs with 
%            the same underlying model.
%
%   [] = FIG_CompGenes(myModel,myGenes,myS,buT,Hb,SaveFigs)
%   myModel : Code for sensitive parameters to stimulus (e.g. '(kON)')
%   myGenes : List of MRW-gene codes to compare (e.g. {'p300e2(2S,300)',
%             'p300DYe2(2S,300)','p300NC(2S,300)'})
%   myS : Replica numbers (e.g. [1:3])
%   buT : Threshold around the maximum likelihood to define the burn-out 
%         period (e.g. 1.005)
%   Hb : Number of bins for 2D-histogram plot (e.g. 25)
%   SaveFigs : If 1, save figures as MRW[myN][myModel]_CompGenes.png,
%              where myN is extracted from the gene lists and corresponds
%              to the number of states considered (i.e. '(2S)' or '(3S)').
%
%   See also SIM_mrwS.m
%   See also FIG_CompModels.m
%   See also FIG_CompRS.m

function [] = FIG_CompGenes(myModel,myGenes,myS,buT,Hb,SaveFigs)
    myC = [0.6 0.8 0;
        1 0.6 0;
        0 0.8 0.8;
        0.08 0.17 0.55;
        0.8 0.6 1;
        1 0 0.4];
    
    myN = strsplit(myGenes{1},{'(',','});
    myN = cat(2,'(',myN{2},')');

    fig = figure();
    fig.Units = 'inches';
    fig.PaperPosition = [1 1 16 4];
    fig.Position = fig.PaperPosition;
    % Title:
    annotation('textbox',[0.15 0.96 0.8 0.04],...
        'String',{cat(2,'MRW',myN,myModel)},...
        'LineStyle','none','HorizontalAlignment','center',...
        'FontWeight','bold','FontSize',12);

    for mrwI = 1:length(myGenes)
        %% Load results:
        for s = myS
            load(cat(2,'MRW_',myGenes{mrwI},myModel,'(s',num2str(s),').mat'))
            Data(s).bu = find([sum(mrw.L,2)>=(max(sum(mrw.L,2))*buT)],1);
            Data(s).L = mrw.L;
            Data(s).P = mrw.P;
        end
        clear s ans 

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

        %% Parameters
        myBs = [0.02 0.008 0.2 0.002 0.02]/10;
        for l = 1:size(Data(1).P,2)
            subplot(1,size(Data(1).P,2),l)
            hold on;
            for s = myS
                i = Data(s).bu:length(Data(s).P);
                [h hI] = hist(Data(s).P(i,l),...
                    [min(Data(s).P(i,l)):myBs(l):max(Data(s).P(i,l))]);
                h = h/length(i);
                plot(hI,h,'LineWidth',2,...
                    'Color',myC(mod(mrwI,length(myC))+1,:),...
                    'DisplayName',myGenes{mrwI})
            end
            xlabel(myPn{l})
            ylabel('Occurence')
        end
        clear l s i h hI
    end

    if(SaveFigs==1)
        print(fig,cat(2,'MRW',myN,myModel,'_CompGenes.png'),'-dpng','-r300')
    end
