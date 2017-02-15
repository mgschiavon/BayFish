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
% FIGURE: Plot phenotype distribution for a given MRW simulation.
%
% Created by Mariana Gómez-Schiavon
% May 2016
%
% FIG_Phenotype : Plot phenotype distribution for a given MRW simulation.
%
%   [] = FIG_Phenotype(myMRW,Hb,Hc,Hm,SaveFigs)
%   myMRW : Results file name (e.g. 'MRW_Fos(2S,300)(kON)(s1)')
%   Hb : Number of bins on the histogram for the probability values (e.g. 
%        100).
%   Hc : Maximum value for the color bar (e.g. 0.25)
%   Hm : Upper limits for the histogram for the probability and mRNA 
%        values (e.g. [0.02,200])
%   SaveFigs : If 1, save figure as [myMRW]_Phenotypes.png.
%
%   See also SIM_Phenotypes.m

function [] = FIG_Phenotype(myMRW,Hb,Hc,Hm,SaveFigs)
    % Load phenotypes & unique parameter sets:
    load(cat(2,myMRW,'_Ms.mat'))
    myT = fieldnames(sL(1).M);
    myPS = size(sL(1).M.t00,1);
    
    %%
    fig = figure();
    fig.Units = 'inches';
    fig.PaperPosition = [0 0 18 10];
    fig.Position = fig.PaperPosition;
    fig.Colormap = colormap('parula');
    fig.Colormap = fig.Colormap(64:-1:1,:);
    fig.Colormap(1,:) = [1 1 1];
    
    for t = 1:length(myT)
        annotation('textbox',[0.13+((t-1)*0.77/length(myT)) 0.92...
            (0.77/length(myT)) 0.03],...
            'String',{myT{t}},...
            'LineStyle','none','HorizontalAlignment','center',...
            'FontWeight','bold','FontSize',12);
        for ps = 1:myPS
            H = zeros(Hb,Hm(2));
            for i = 1:length(iP)
                p = sL(i).M.(myT{t})(ps,:);
                h = histcounts2(p,[1:length(p)]-1,[0:(Hm(1)/Hb):Hm(1)],[0:Hm(2)]);
                H = H + (iP(i)*h);
            end
            H = H/sum(iP);
            subplot(myPS,length(myT),(t+((ps-1)*length(myT))))
                imagesc([0:Hm(2)],[0:(Hm(1)/Hb):Hm(1)],H,[0 Hc])
                    set(gca,'YDir','normal')
                    xlim([0 Hm(2)])
                    ylim([0 Hm(1)])
                    if(ps==myPS)
                        xlabel('mRNA')
                    end
                    if(t==1)
                        ylabel(cat(2,'P(promoter state [',num2str(ps),...
                            '], # mRNAs)'))
                    end
        end
    end
    clear t ps H i p h
    
    % Title:
    annotation('textbox',[0.13 0.95 0.77 0.03],...
        'String',{strrep(myMRW,'_','\_')},...
        'LineStyle','none','HorizontalAlignment','center',...
        'FontWeight','bold','FontSize',12);
        
    if(SaveFigs==1)
        print(fig,cat(2,myMRW,'_Phenotypes.png'),'-dpng','-r300')
    end
