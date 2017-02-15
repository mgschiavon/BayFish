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
% FIGURE: Plot estimated distribution vs data.
%
% Created by Mariana Gómez-Schiavon
% May 2016
%
% FIG_MData :  Plot estimated distribution vs data.
%
%   [] = FIG_MData(M,x,AxLim)
%   M : Probability distribution.
%   x : Data occurence (same structure as M).
%   AxLim : Upper limit for the axes ineach plot, [x-axis,yLeft-axis,
%           yRight-axis] (e.g. [300,30,0.05]).
%
%   See also FIG_Data.m

function [] = FIG_MData(M,x,AxLim)
    myT = fieldnames(x);
    myPS = size(x.t00,1);
    myMs = [1:size(x.t00,2)]-1;
    
    fig = figure();
    fig.Units = 'inches';
    fig.PaperPosition = [0 0 18 10];
    fig.Position = fig.PaperPosition;
    fig.Colormap = colormap('parula');
    fig.Colormap = fig.Colormap(64:-1:1,:);
    fig.Colormap(1,:) = [1 1 1];
    left_color = [0 0 0];
    right_color = [0.8 0 0.4];
    set(fig,'defaultAxesColorOrder',[left_color; right_color]);
    
    for t = 1:length(myT)
        annotation('textbox',[0.13+((t-1)*0.77/length(myT)) 0.92...
            (0.77/length(myT)) 0.03],...
            'String',{cat(2,myT{t},' | n = ',...
            num2str(sum(sum(x.(myT{t})))))},...
            'LineStyle','none','HorizontalAlignment','center',...
            'FontWeight','bold','FontSize',12);
        for ps = 1:myPS
            a = x.(myT{t})(ps,:);
            b = zeros(sum(a),2);
            i = 0;
            for m = 1:length(a)
                for j = 1:a(m)
                    i = i + 1;
                    b(i,:) = [m-1,j];
                end
            end
            clear m i j
            subplot(myPS,length(myT),(t+((ps-1)*length(myT))))
            yyaxis left
                scatter(b(:,1),b(:,2),...
                    'MarkerFaceColor',[0 0.45 0.74],...
                    'MarkerEdgeColor',[0.08 0.17 0.55])
                    xlim([0 AxLim(1)])
                    ylim([0 AxLim(2)])
                    if(t==1)
                        ylabel('Occurence')
                    end
            yyaxis right
                plot(myMs,M.((myT{t}))(ps,:),...
                    'LineWidth',2,'Color',[0.8 0 0.4])
                    xlim([0 AxLim(1)])
                    ylim([0 AxLim(3)])
                    if(ps==myPS)
                        xlabel('mRNA')
                    end
                    if(t==length(myT))
                        ylabel(cat(2,'P([\rho[',num2str(ps),...
                            '], # mRNAs])'))
                    end
                    set(gca,'FontSize',14)
        end
    end
    clear t ps
