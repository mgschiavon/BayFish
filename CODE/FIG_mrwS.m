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
% FIGURE: Plot and compare MRW replicas.
%
% Created by Mariana G�mez-Schiavon
% May 2016
%
% FIG_mrwS : Plot and compare MRW replicas.
%
%   [] = FIG_mrwS(myMRW,myS,buT,Hb,SaveFigs)
%   myMRW : Results files common name (e.g. 'MRW_Fos(2S,300)(kON)')
%   myS : Replica numbers (e.g. [1:3])
%   buT : Threshold around the maximum likelihood to define the burn-out 
%         period (e.g. 1.005)
%   Hb : Number of bins for 2D-histogram plot (e.g. 25)
%   SaveFigs : If 1, save figures as [myMRW]_LogL.png, and [myMRW]_Par.png.
%
%   See also SIM_mrw.m

function [] = FIG_mrwS(myMRW,myS,myT,buT,Hb,SaveFigs)
    % Load results:
    for s = myS
        load(cat(2,myMRW,'(s',num2str(s),').mat'))
        [min(mrw.e) max(mrw.e)]
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

    %% LogL
    fig = figure();
    fig.Units = 'inches';
    fig.PaperPosition = [1 1 16 4];
    fig.Position = fig.PaperPosition;
    myC = colormap('parula');
    % Title:
    annotation('textbox',[0.15 0.96 0.8 0.04],...
        'String',{strrep(myMRW,'_','\_')},...
        'LineStyle','none','HorizontalAlignment','center',...
        'FontWeight','bold','FontSize',12); 
    for l = 1:size(Data(1).L,2)
        subplot(1,size(Data(1).L,2),l)
        hold on;
        for s = myS
            i = Data(s).bu:length(Data(s).L);
            [h hI] = hist(Data(s).L(i,l),Hb);
            h = h/length(i);
            plot(hI,h,'LineWidth',2,'Color',myC(s*floor(64/length(myS)),:))
        end
        xlabel(cat(2,'logL_{',myT{l},'}'))
        ylabel('Occurence')
    end
    clear l s i h hI
    if(SaveFigs==1)
        print(fig,cat(2,myMRW,'_LogL.png'),'-dpng','-r300')
    end

    %% Parameters
    fig = figure();
    fig.Units = 'inches';
    fig.PaperPosition = [1 1 16 4];
    fig.Position = fig.PaperPosition;
    % Title:
    annotation('textbox',[0.15 0.96 0.8 0.04],...
        'String',{strrep(myMRW,'_','\_')},...
        'LineStyle','none','HorizontalAlignment','center',...
        'FontWeight','bold','FontSize',12);
    
    
    myBs = [0.02 0.02 1 0.02 0.2]/20;
    
    
    for l = 1:size(Data(1).P,2)
        subplot(1,size(Data(1).P,2),l)
        hold on;
        for s = myS
            i = Data(s).bu:length(Data(s).P);
            [h hI] = hist(Data(s).P(i,l),...
                [min(Data(s).P(i,l)):myBs(l):max(Data(s).P(i,l))]);
            h = h/length(i);
            plot(hI,h,'LineWidth',2,'Color',myC(s*floor(64/length(myS)),:))
        end
        xlabel(myPn{l})
        ylabel('Occurence')
    end
    clear l s i h hI
    if(SaveFigs==1)
        print(fig,cat(2,myMRW,'_Par.png'),'-dpng','-r300')
    end
