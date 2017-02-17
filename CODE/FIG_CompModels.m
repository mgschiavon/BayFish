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
% FIGURE: Plot and compare MRW log-likelihood for different models.
%
% Created by Mariana Gómez-Schiavon
% August 2016
%
% FIG_CompModels : Plot and compare MRW log-likelihood for different models.
%
%   [] = FIG_CompModels(myGene,myM,myS,buT,Hb,xT,SaveFigs)
%   myGene : MRW-gene code (e.g. 'Fos(2S,300)')
%   myM : Models (e.g.
%   {'kON','kOFF','mu','kONkOFF','kONmu','kOFFmu','kONkOFFmu'}).
%   myS : Replica numbers (e.g. [1:3])
%   buT : Threshold around the maximum likelihood to define the burn in 
%         period (e.g. 1.005)
%   Hb : Number of bins for 2D-histogram plot (e.g. 25)
%   xT : Flag to plot individual time measurements (1) or total
%        log-likelihood (0).
%   SaveFigs : If 1, save figures as MRW[myGene]_CompModels_xT.png, if 
%              xT is 1, MRW[myGene]_CompModels.png, otherwise.
%
%   See also SIM_mrwS.m
%   See also FIG_CompGenes.m
%   See also FIG_CompRS.m

function [] = FIG_CompModels(myGene,myM,myS,buT,Hb,xT,SaveFigs)
    myC = [0.6 0.8 0;
        1 0.6 0;
        0 0.8 0.8;
        0.08 0.17 0.55;
        0.8 0.6 1;
        1 0 0.4];
    
    % Define time labels
    load(cat(2,'myData_',myGene,'.mat'),'x')
    myT = fieldnames(x);
    mySx = zeros(1,length(myT));
    for t = 1:length(myT)
        mySx(t) = sum(sum(x.(myT{t})));
    end
    clear t x

    fig = figure();
    fig.Units = 'inches';
    if(xT==1)
        fig.PaperPosition = [1 1 16 4];
    else
        fig.PaperPosition = [1 1 6 4];
    end
    fig.Position = fig.PaperPosition;
    % Title:
    annotation('textbox',[0.15 0.96 0.8 0.04],...
        'String',{cat(2,'MRW : ',myGene)},...
        'LineStyle','none','HorizontalAlignment','center',...
        'FontWeight','bold','FontSize',12);

    for m = 1:length(myM)
        %% Load results:
        for s = myS
            load(cat(2,'MRW_',myGene,'(',myM{m},')(s',num2str(s),').mat'));
            Data(s).bu = find([sum(mrw.L,2)>=(max(sum(mrw.L,2))*buT)],1);
            Data(s).L = mrw.L;
            Data(s).P = mrw.P;
        end
        clear s ans

        %% Log-likelihood
        if(xT==1)
            for l = 1:size(Data(1).L,2)
                subplot(1,size(Data(1).L,2),l)
                hold on;
                for s = myS
                    i = Data(s).bu:length(Data(s).L);
                    [h hI] = hist(Data(s).L(i,l)/mySx(l),Hb);
                    h = h/length(i);
                    plot(hI,h,'LineWidth',2,...
                        'Color',myC(mod(m,length(myC))+1,:),...
                        'DisplayName',myM{m})
                end
                xlabel(cat(2,'logL_{',myT{l},'}'))
                ylabel('Occurence')
            end
            clear l s i h hI
        else
            hold on;
            for s = myS
                i = Data(s).bu:length(Data(s).L);
                [h hI] = hist(Data(s).L(i),Hb);
                h = h/length(i);
                plot(hI,h,'LineWidth',2,...
                    'Color',myC(mod(m,length(myC))+1,:),...
                    'DisplayName',myM{m})
            end
            xlabel(cat(2,'\langlelogL\rangle'))
            ylabel('Occurence')
            clear l s i h hI
        end
    end
    legend('show')

    if(SaveFigs==1)
        if(xT==1)
            print(fig,cat(2,'MRW',myGene,'_CompModels_xT.png'),'-dpng','-r300')
        else
            print(fig,cat(2,'MRW',myGene,'_CompModels.png'),'-dpng','-r300')
        end
    end
