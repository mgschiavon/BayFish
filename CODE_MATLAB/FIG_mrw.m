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
% FIGURE: Plot Metropolis-Hastings results.
%
% Created by Mariana Gómez-Schiavon
% May 2016
%
% FIG_mrw : Plot Metropolis-Hastings results.
%
%   [] = FIG_mrw(myMRW,buT,Hc,Hb,cj,SaveFigs)
%   myMRW : Results file name (e.g. 'MRW_Fos(2S,300)(kON)(s1)')
%   buT : Threshold around the maximum likelihood to define the burn in 
%         period (e.g. 1.005)
%   Hc : Maximum population fraction for color bar in 2D-histogram plot 
%        (e.g. 1/50).
%   Hb : Number of bins for 2D-histogram plot (e.g. 25)
%   cj : Correlation jump, i.e. sample time traces every cj iterations 
%        (e.g. 1)
%   SaveFigs : If 1, save figures as [myMRW]_LogL.png, [myMRW]_Par.png, 
%              and [myMRW]_Hist.png.
%
%   See also SIM_mrw.m

function [] = FIG_mrw(myMRW,buT,Hc,Hb,cj,SaveFigs)
    % Load results:
    load(myMRW)
    if(isempty(regexp(myMRW,'TEMP_','match'))==0)
        mrw.I = j0 - 1;
        mrw.P = mrw.P([1:mrw.I],:);
        mrw.L = mrw.L([1:mrw.I],:);
        mrw.e = mrw.e([1:mrw.I],:);
        clear r j0
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
    clear i I1
    
    % Determine burn-out period:
    L = sum(mrw.L,2);
    bu = find([L>=(max(L)*buT)],1);
    
    %% Log-likelihood & burn-in
    fig = figure();
    fig.Units = 'inches';
    fig.PaperPosition = [1 1 5 4];
    fig.Position = fig.PaperPosition;
    % Title:
    annotation('textbox',[0.15 0.94 0.8 0.05],...
        'String',{strrep(myMRW,'_','\_')},...
        'LineStyle','none','HorizontalAlignment','center',...
        'FontWeight','bold','FontSize',12);
    % Burn-out period:
    axes('Position',[0.15 0.575 0.8 0.325])
    hold on;
        plot(L([1:min(round(bu*1.2),length(L))]),'DisplayName','\langleL\rangle',...
            'LineWidth',2,'Color',[0.5 0.5 0.5])
        plot([bu bu],[0 min(L)],'DisplayName','burn-out',...
            'LineWidth',2,'LineStyle',':','Color',[1 0 0])
        xlim([0 min(round(bu*1.2),length(L))])
        ylim([min(L) max(L)*0.95])
        xlabel('Iteration')
    % Log-Likelihood:
    axes('Position',[0.15 0.125 0.8 0.325])
    hold on;
        plot(L,'DisplayName','\langleL\rangle',...
            'LineWidth',2,'Color',[0.5 0.5 0.5])
        plot([bu bu],[0 min(L)],'DisplayName','burn-out',...
            'LineWidth',2,'LineStyle',':','Color',[1 0 0])
        xlim([0 mrw.I])
        ylim([min(L) max(L)*0.95])
        xlabel('Iteration')
        legend('show','Location','southeast')
    if(SaveFigs==1)
        print(fig,cat(2,myMRW,'_LogL.png'),'-dpng','-r300')
    end
    
    %% Parameters:
    fig = figure();
    fig.Units = 'inches';
    fig.PaperPosition = [1 0 6 12];
    fig.Position = fig.PaperPosition;
    myC = colormap('parula');
    % Title:
    annotation('textbox',[0.15 0.94 0.8 0.05],...
        'String',{strrep(myMRW,'_','\_')},...
        'LineStyle','none','HorizontalAlignment','center',...
        'FontWeight','bold','FontSize',12);
    % Parameters:
    pl = length(myPn);
    for p = 1:length(myPn)
        axes('Position',[0.125,0.05+((pl-p)*((0.8/pl)+0.025)),0.85,(0.8/pl)])
        plot([bu:cj:mrw.I],mrw.P([bu:cj:mrw.I],p),...
            'LineWidth',1,'Color',myC(floor(64/length(myPn))*p,:))
        ylabel(myPn{p})
        set(gca,'XTick',[10000:10000:mrw.I],'XTickLabel','')
    end
    set(gca,'XTick',[10000:10000:mrw.I],'XTickLabel',[10000:10000:mrw.I])
    xlabel('Iterations')
    clear pl p myC
    if(SaveFigs==1)
        print(fig,cat(2,myMRW,'_Par.png'),'-dpng','-r300')
    end
        
    %% Histograms
    fig = figure();
    fig.Units = 'inches';
    fig.PaperPosition = [2 0 12 12];
    fig.Position = fig.PaperPosition;
    fig.Colormap = colormap('parula');
    fig.Colormap = fig.Colormap(64:-1:1,:);
    
    myPn = [{'logL'};myPn];
    pl = length(myPn);
    x = [L([bu:cj:mrw.I]),mrw.P([bu:cj:mrw.I],:)];
    mX = [min(x); max(x)];
    
    % Title:
    annotation('textbox',[0.15 0.97 0.8 0.03],...
        'String',{strrep(myMRW,'_','\_')},...
        'LineStyle','none','HorizontalAlignment','center',...
        'FontWeight','bold','FontSize',12);

    for i = 1:pl
        for j = 1:pl
            k = 0.9/pl;
            axes('Position',[0.075+((i-1)*k) 0.075+((j-1)*k) k k])
            h = histogram2(x(:,i),x(:,j),...
                [mX(1,i):((mX(2,i)-mX(1,i))/Hb):mX(2,i)],...
                [mX(1,j):((mX(2,j)-mX(1,j))/Hb):mX(2,j)],...
                'FaceColor','flat','DisplayStyle','tile');
            set(gca,'CLim',[0 length(x)*Hc],...
                'XGrid','off','YGrid','off','Box','on')
            xlim([mX(1,i) mX(2,i)])
            ylim([mX(1,j) mX(2,j)])
            if(j==1)
                xlabel(myPn{i})
            else
                set(gca,'XTick',[])
            end
            if(i==1)
                ylabel(myPn{j})
            else
                set(gca,'YTick',[])
            end

            if(j==i)
                axes('Position',[0.075+((i-1)*k) 0.075+((j-1)*k) k k])
                [h,hI] = hist(x(:,i),[mX(1,i):((mX(2,i)-mX(1,i))/Hb):mX(2,i)]);
                bar(hI,h/length(x),'FaceColor',[0.078 0.169 0.549],...
                    'EdgeColor',[0.6 0.6 0.6])
                xlim([mX(1,i) mX(2,i)])
                ylim([0 0.2])
                set(gca,'XTick',[],'YTick',[])
            end
        end
    end

    axes('Position',[0.075+((i-1)*k)+0.015 0.075+((j-1)*k)+(k-0.025) (k-0.025) 0.005])
        scatter([1:64],zeros(1,64),30,fig.Colormap,'filled')
            xlim([1 64])
            xlabel('%','Position',[0 -2.25 0])
            set(gca,'XTick',[16:16:64],...
                'XTickLabel',100*[1:4]*Hc/4)
    clear i j k h hI
    if(SaveFigs==1)
        print(fig,cat(2,myMRW,'_Hist.png'),'-dpng','-r300')
    end
