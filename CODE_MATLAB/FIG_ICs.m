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
% FIGURE: Calculate information criteria and compare them.
%
% Created by Mariana Gómez-Schiavon
% June 2016
%
% FIG_ICs : Plot and compare BIC, AIC, and DIC metrics.
%
%   [] = FIG_ICs(myGeneModel,myM,buT,sT,SaveFig)
%   myGeneModel : Gene model to use (e.g. 'Fos(2S,300)')
%   myM : Models (e.g. {'kON','kOFF','mu','kONkOFF','kONmu','kOFFmu','kONkOFFmu'})
%   buT : Threshold around the maximum likelihood to define the burn in 
%         period (e.g. 1.005)
%   sT : Number of replicas (e.g. 3)
%   SaveFig : If 1, save figure as MRW_[myGeneModel]_ICs.png.
%
%   See also DATA_BIC.m
%   See also DATA_AIC.m
%   See also DATA_DIC.m

function [] = FIG_ICs(myGeneModel,myM,buT,sT,SaveFig)
    %% Load data:
    load(cat(2,'myData_',myGeneModel,'.mat'),'x')
    myT = fieldnames(x);
    n = 0;
    for t = 1:length(myT)
        n = n + sum(sum(x.(myT{t})));
    end
    clear t

    %% Calculate metrics:
    BIC = zeros(length(myM),sT);
    AIC = zeros(length(myM),sT);
    DIC = zeros(length(myM),sT);

    fig = figure();
        fig.Units = 'inches';
        fig.PaperPosition = [2 1 6 3];
        fig.Position = fig.PaperPosition;
        left_color = [0 0 0];
        right_color = [.5 .5 .5];
        set(fig,'defaultAxesColorOrder',[left_color; right_color]);

    for m = 1:length(myM)
        for s = 1:sT
            load(cat(2,'MRW_',myGeneModel,'(',myM{m},')(s',num2str(s),').mat'));
    
            L = sum(mrw.L,2);
            bu = find([L>=(max(L)*buT)],1);
            maxL = max(L(bu:length(L)));
            clear L bu
    
            k = length(fieldnames(mrw.U)) + (2*length(fieldnames(mrw.S)));
            BIC(m,s) = DATA_BIC(maxL,k,n);
            AIC(m,s) = DATA_AIC(maxL,k,n,0);
            DIC(m,s) = DATA_DIC(mrw,Par,ParS,myGeneModel,buT);
        end

        hold on;
        yyaxis left
            plot(zeros(1,sT)+m,BIC(m,:),'DisplayName','BIC',...
                'Marker','o','LineWidth',1,'LineStyle','none','Color',[0.6 0.8 0])
            plot(zeros(1,sT)+m,AIC(m,:),'DisplayName','AIC',...
                'Marker','s','LineWidth',1,'LineStyle','none','Color',[0.8 0 0.6])
            ylabel('BIC,AIC');
        yyaxis right
            plot(zeros(1,sT)+m,DIC(m,:),'DisplayName','DIC',...
                'Marker','v','LineWidth',1,'LineStyle','none','Color',[0 0.6 0.8])
            ylabel('DIC');
            if(m==1)
                legend('show')
            end
    end
    clear m s k 
    a = myM; a = strrep(a,'ON','_{ON}'); a = strrep(a,'OFF','_{OFF}');
    a = strrep(a,'mu','\mu'); a = strrep(a,'}k','},k');
    a = strrep(a,'}\','},\');
    set(gca,'XTick',[1:length(myM)],'XTickLabel',a)
    xlim([0.8 length(myM)+0.2])
    set(gca,'YGrid','on')
    
    if(SaveFig==1)
        print(fig,cat(2,'MRW_',myGeneModel,'_ICs.png'),'-dpng','-r300')
    end
