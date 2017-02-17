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
% SIMULATE: Run Metropolis-Hastings & obtain "best" kinetic parameters.
%
% Created by Mariana Gómez-Schiavon
% May 2016
%
% SIM_mrw : Run the Metropolis Random Walk algorithm over a particular 
%           model and experimental data.
%
%   [] = SIM_mrw(myGeneModel,Par,mrw,cont)
%   myGeneModel : Gene model to simulate (e.g. 'Fos(2S,300)')
%   Par : Structure with the kinetic parameters.
%      .kON : Promoter activation rate (OFF->ON)
%      .kOFF : Promoter deactivation rate (ON->OFF)
%      .kONs : Promoter "super" activation rate (ON->ONs)
%      .kOFFs : Promoter "super" deactivation rate (ONs->ON)
%      .mu0 : mRNA synthesis rate of promoter in OFF state
%      .mu : mRNA synthesis rate of promoter in ON state
%      .muS : mRNA synthesis rate of promoter in ONs state
%      .d : mRNA degradation rate
%   mrw : Structure with the Metropolis Random Walk algorithm parameters:
%      .I : Number of iterations (e.g. 100000)
%      .s : Random number generator seed (e.g. 1)
%      .([U/S]).(par) : [sigma,min,max] (e.g. [1e-4,1e-3,10])
%         [U/S] : S if the parameter changes between 'basal' and
%                 'stimulus'; U, otherwise.
%         sigma : Covariance of (par) to calculate parameter
%                 perturbations.
%         min,max : Initial value range of (par).
%      .P : OUTPUT - Parameters per iteration
%      .L : OUTPUT - Log-likelihood per iteration
%      .e : OUTPUT - Maximum error observed in the distribution estimation
%   cont : If 1, continue stopped simulation saved on file 
%            'TEMP_MRW_[myGene]([N],[maxM],[fieldnames(ParS)])(s[mrw.s]).mat'
%
%   OUTPUT FILE : 'MRW_[myGene]([N],[maxM],[fieldnames(ParS)])(s[mrw.s]).mat'
%
%   See also SIM_logL_Opt.m
%   See also SIM_logL_SLOW.m

function [] = SIM_mrw(myGeneModel,Par,mrw,cont)
    % Model (N), maximum mRNA number to consider (maxM), and data (x):
    load(cat(2,'myData_',myGeneModel,'.mat'),'N','maxM','x');

    % Define indexes:
    myI = [];
    I1 = fieldnames(mrw.U);
    for i = 1:length(I1)
        myI = [myI;{'U',I1{i},1}];
    end
    I1 = fieldnames(mrw.S);
    psN = ''; % ParS fieldnames to name results.
    for i = 1:length(I1)
        myI = [myI;{'S',I1{i},1};{'S',I1{i},2}];
        psN = cat(2,psN,I1{i});
    end
    clear i I1
    % Specify time points (myT) and sample size (myS):
    myT = fieldnames(x);
    myS = zeros(1,length(myT));
    for t = 1:length(myT)
        myS(t) = sum(sum(x.(myT{t})));
    end
    clear t
    
    if(cont==1)
        load(cat(2,'TEMP_MRW_',myGeneModel,'(',psN,')',...
            '(s',num2str(mrw.s),')(EfosB).mat'));
    else
        % (1) Define random number generator:
        rng(mrw.s,'twister');
        r.p0 = rand(size(myI,1),1);
        for i = 1:size(myI,1)
            Z(i) = mrw.(myI{i,1}).(myI{i,2})(myI{i,3},1);
        end
        r.Pe = mvnrnd(zeros(mrw.I,length(Z)),Z);
        r.tL = rand(mrw.I,1);
        clear i Z

        % (2) Initialize system:
        mrw.P = zeros(mrw.I,size(myI,1));
        mrw.L = zeros(mrw.I,length(myT));
        mrw.e = zeros(mrw.I,1);
        % Initial parameter set:
        for i = 1:size(myI,1)
            if(strcmp(myI{i,1},'U'))
                Par.(myI{i,2}) = mrw.U.(myI{i,2})(2) + ...
                    (r.p0(i)*(mrw.U.(myI{i,2})(3)-mrw.U.(myI{i,2})(2)));
                mrw.P(1,i) = Par.(myI{i,2});
            else
                ParS.(myI{i,2})(myI{i,3}) = mrw.S.(myI{i,2})(myI{i,3},2)...
                    + (r.p0(i)...
                    *(mrw.S.(myI{i,2})(myI{i,3},3)-mrw.S.(myI{i,2})(myI{i,3},2)));
                mrw.P(1,i) = ParS.(myI{i,2})(myI{i,3});
            end
        end
        [m,l,e] = SIM_logL_Opt(x,N,maxM,Par,ParS);
        for t = 1:length(myT)
            mrw.L(1,t) = myS(t)*l.(myT{t});
            mrw.e(1) = max(mrw.e(1),abs(e.(myT{t})));
        end
        clear i t m l e
        j0 = 2;
    end

    % (3) Iterate:
    for j = j0:mrw.I
        mrw.P(j,:) = mrw.P(j-1,:);
        mrw.L(j,:) = mrw.L(j-1,:);
        mrw.e(j) = mrw.e(j-1,:);
        % Alternative parameter set:
        myP = mrw.P(j,:) + r.Pe(j,:);
        if(min(myP)>=1e-8)
            for i = 1:size(myI,1)
                if(strcmp(myI{i,1},'U'))
                    Par.(myI{i,2}) = myP(i);
                else
                    ParS.(myI{i,2})(myI{i,3}) = myP(i);
                end
            end
            % Generate proposal:
            [m,l,e] = SIM_logL_Opt(x,N,maxM,Par,ParS);
            myL = zeros(1,length(myT));
            myE = 0;
            for i = 1:length(myT)
                myL(i) = myS(i)*l.(myT{i});
                myE = max(myE,abs(e.(myT{i})));
            end
            clear i t m l e
            % If proposal is accepted, update system:
            if(r.tL(j) <= min(1,exp(sum(myL)-sum(mrw.L(j,:)))))
                mrw.P(j,:) = myP;
                mrw.L(j,:) = myL;
                mrw.e(j) = myE;
            end
        end
        % Save progress:
        if(mod(j,100)==0)
            j0 = j + 1
            save(cat(2,'TEMP_MRW_',myGeneModel,'(',psN,')',...
            '(s',num2str(mrw.s),').mat'),...
                'mrw','j0','r','Par','ParS');
        end
    end
    clear j myP myL myE

    % (4) Save:
    save(cat(2,'MRW_',myGeneModel,'(',psN,')',...
            '(s',num2str(mrw.s),').mat'),...
        'mrw','Par','ParS');
    delete(cat(2,'TEMP_MRW_',myGeneModel,'(',psN,')',...
            '(s',num2str(mrw.s),').mat'));
