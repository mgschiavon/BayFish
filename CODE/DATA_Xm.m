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
% DATA: Create data matrix
%
% Created by Mariana Gómez-Schiavon
% May 2016
%
% DATA_Xm : Sort the data in a matrix form where each row represents a 
%           promoter state (assuming either a 2-states or a 3-states model) 
%           and each row the number of free mRNA molecules in the cell.
%
%   Xm = DATA_Xm(X,N,maxM,a)
%   X : Data list (['TS1','TS2','mRNA']).
%   N : Model ('2S' or '3S').
%   maxM : Maximum mRNA number to consider.
%   a : If N='3S', threshold to define a TS super active, i.e. ONs (e.g. 10).
%
%   Xm : Observed data matrix where Xm(i,j) is the number of 
%        individuals in de population with exactly promoters in (i) state
%        and (j-1) mRNA molecules.
%        If N='2S', the possible promoters states are 1:{OFF,OFF},
%        2:{OFF,ON}, and 3:{ON,ON}.
%        If N='3S', the possible promoters states are 1:{OFF,OFF}, 
%        2:{OFF,ON}, 3:{ON,ON}, 4:{OFF,ONs}, 5:{ON,ONs}, and 6:{ONs,ONs}.
%
%   See also DATA_X.m

function Xm = DATA_Xm(X,N,maxM,a)
    if(strcmp(N,'2S'))
        X = [([X(:,1)>0]+[X(:,2)>0]) X(:,3)];
        Xm = zeros(3,maxM+1);
        for i = 1:length(X)
            Xm(X(i,1)+1,X(i,2)+1) = Xm(X(i,1)+1,X(i,2)+1) + 1;
        end
        clear i
    elseif(strcmp(N,'3S'))
        X = [([X(:,1)>0 & X(:,1)<=a]+[X(:,2)>0 & X(:,2)<=a]),...
            ([X(:,1)>a]+[X(:,2)>a]),X(:,3)];
        Xm = zeros(6,maxM+1);
        for i = 1:length(X)
            if(X(i,1)==0 & X(i,2)==0)
                Xm(1,X(i,3)+1) = Xm(1,X(i,3)+1) + 1;
            elseif(X(i,1)==1 & X(i,2)==0)
                Xm(2,X(i,3)+1) = Xm(2,X(i,3)+1) + 1;
            elseif(X(i,1)==2 & X(i,2)==0)
                Xm(3,X(i,3)+1) = Xm(3,X(i,3)+1) + 1;
            elseif(X(i,1)==0 & X(i,2)==1)
                Xm(4,X(i,3)+1) = Xm(4,X(i,3)+1) + 1;
            elseif(X(i,1)==1 & X(i,2)==1)
                Xm(5,X(i,3)+1) = Xm(5,X(i,3)+1) + 1;
            elseif(X(i,1)==0 & X(i,2)==2)
                Xm(6,X(i,3)+1) = Xm(6,X(i,3)+1) + 1;
            end
        end
        clear i
    else
        cat(2,'ERROR: Model ',N,' has not been defined.')
    end
end
