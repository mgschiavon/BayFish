# BayFish - MATLAB version
Bayesian inference of transcription dynamics from population snapshots of smFISH

BayFish is a computational pipeline to infer kinetic parameters of gene expression from sparse single-molecule RNA fluorescence *in situ* hybridization (smFISH) data at multiple time points after induction. Given an underlying model of gene expression, BayFish uses a Markov Chanin Monte Carlo method to estimate the posterior probability of the model parameters and quantify the parameter uncertainty given the observed smFISH data.

NOTE: To use the C++ version, please refer to CODE_CPP/README.md

## Instructions
BayFish pipeline general description. A name or flag must be assigned to the used data set (*myGene*, e.g. `Npas4`).

1. Define data.
2. Define mathematical model.
3. Calculate probability distributions.
4. Compare model & data (i.e. calculate/define posteriors).
5. Run Metropolis-Hastings & obtain "best" kinetic parameters.
6. Explore sampling effect (i.e. synthetic data).

### (1) Define data:

(1.1) Import data using DATA_X.m. Experimental data should be in an excel file (e.g. `DataMerged_Npas4.xlsx`) where each time point measure is a sheet (e.g. `t15`), and the first line has the "name" of the experiment (e.g. `KCl (5 min) + CM (10 min)`), the second line the name of the columns (`TS1`, `TS2`, `mRNA`, and, optionally,  `Exp.` with the experiment reference or other notes), and then the list of individual cell measurements (e.g. `[6.55,4.76,46,'II']`).

(1.2) Convert data into a matrix form.

(1.3) Save file.

```matlab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Exp = cat(2,'DataMerged_',myGene,'.xlsx');  % Experimental data file
NtZ = 0;                                    % If 1, substitute NaN to zero.
ExI = 0;                                    % If 1, extract experiment information.
N = '2S';                                   % Model ('2S' or '3S').
maxM = 300;                                 % Maximum mRNA number to consider.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% (1.1) Import data:
    [i,t] = xlsfinfo(Exp);
    for i = 1:length(t)
        X = DATA_X(Exp,t{i},NtZ,ExI);
        
        % (1.2) Convert to matrix:
        x.(t{i}) = DATA_Xm(X.data,N,maxM);
    end
    clear i t X

% (1.3) Save file:
    save(cat(2,'myData_',myGene,'(',N,',',num2str(maxM),').mat'),...
        'x','N','maxM')
    clear x
```

### (2) Define mathematical model:

(2.1) If needed, create instructions for the model's transition matrix given the chosen limit for mRNA molecule number.

(2.2) Calculate the transition matrixes given the model and its kinetic parameters.

```matlab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model (N) and maximum mRNA number to consider (maxM):
load(cat(2,'myData_',myGeneModel,'.mat'),'N','maxM');
% Kinetic parameters:
Par.kON = NaN;      % Promoter activation rate (OFF->ON)
Par.kOFF = 0.0616;	% Promoter deactivation rate (ON->OFF)
Par.kONs = NaN;     % Promoter "super" activation rate (ON->ONs)
Par.kOFFs = NaN;	% Promoter "super" deactivation rate (ONs->ON)
Par.mu0 = 0.0266;	% mRNA synthesis rate of promoter in OFF state
Par.mu = 5.0851;	% mRNA synthesis rate of promoter in ON state
Par.muS = NaN;   	% mRNA synthesis rate of promoter in ONs state
Par.d = 0.0462;     % mRNA degradation rate
% Kinetic parameters sensitive to stimulus:
ParS.kON = [0.0051,0.1373]; % [Basal,Stimulus]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% (2.1) Model structure:
    STRUCT_TransM(N,maxM);
    
% (2.2) Calculate transition matrixes:
    pS = fieldnames(ParS);
    % Basal conditions:
    for i1 = 1:length(pS)
        Par.(pS{i1}) = ParS.(pS{i1})(1);
    end
    A.b = DATA_TransM(N,maxM,Par);
    % Stimulus conditions:
    for i1 = 1:length(pS)
        Par.(pS{i1}) = ParS.(pS{i1})(2);
    end
    A.s = DATA_TransM(N,maxM,Par);
    clear pS i1
```

### (3) Calculate probability distributions

(3.1) Calculate the probability distributions for the chosen time points the transition matrixes.

(3.2) Transform the probability distributions vectors to a matrix form where M(i,j) is the probability of having the promoters state (i) and exactly (j-1) mRNA molecules.

```matlab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model (N):
load(cat(2,'myData_',myGeneModel,'.mat'),'N');
% Time points:
load(cat(2,'myData_',myGeneModel,'.mat'),'x');
myT = fieldnames(x);
clear x
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% (3.1) Calculate probability distributions vectors:
    P.t00 = DATA_Pss(A.b);
    for t = 2:length(myT)
        tt = regexp(myT{t},'\.*\d*','match');
        P.(myT{t}) = DATA_Pxt(A.s,P.t00,str2double(tt{1}));
        % Calculate FSP algorithm estimation error:
        eps.(myT{t}) = 1-sum(P.(myT{t}));
        if(eps.(myT{t})>0.001)
            cat(2,'CAUTION: Large FSP algorithm estimation error (',...
                num2str(eps.(myT{t})),')')
        end
    end
    clear t tt
    
% (3.2) Transform probability distributions vectors to matrixes:
    for t = 1:length(myT)
        M.(myT{t}) = DATA_P2M(N,P.(myT{t}));
    end
    clear t
```

### (4) Compare model & data:

(4.1) Calculate the log-likelihood of observing the data x given the probability distribution M.

```matlab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data:
load(cat(2,'myData_',myGeneModel,'.mat'),'x');
% Time points:
myT = fieldnames(x);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% (4.1) Calculate the log-likelihood:
    for t = 1:length(myT)
        L.(myT{t}) = DATA_logL(M.(myT{t}),x.(myT{t}));
    end
    clear t
```

### (5) Run Metropolis-Hastings & obtain "best" kinetic parameters:

(5.1) Define random numbers.

(5.2) Initialize system.

(5.3) Iterate.

(5.4) Save.

```matlab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model (N), maximum mRNA number to consider (maxM), and data (x):
load(cat(2,'myData_',myGeneModel,'.mat'),'N','maxM','x');
% Kinetic parameters:
Par.kON = NaN;          % Promoter activation rate (OFF->ON)
Par.kOFF = NaN;         % Promoter deactivation rate (ON->OFF)
Par.mu0 = NaN;          % mRNA synthesis rate of promoter in OFF state
Par.mu = NaN;           % mRNA synthesis rate of promoter in ON state
Par.kONs = NaN;         % Promoter "super" activation rate (ON->ONs)
Par.kOFFs = NaN;        % Promoter "super" deactivation rate (ONs->ON)
Par.muS = NaN;          % mRNA synthesis rate of promoter in ONs state
Par.d = 0.0559;         % mRNA degradation rate
% Metropolis Random Walk algorithm parameters:
mrw.I = 100000;         % Number of iterations
mrw.s = 1;              % Random number generator seed
%%% mrw.([U/S]).(par) = [sigma,min,max] %%%
mrw.U.kOFF = [1e-5,	1e-5, 0.1];
mrw.U.mu0  = [1e-5,	1e-5, 0.1];
mrw.U.mu   = [1e-3,	1e-3, 10];
mrw.S.kON  = [1e-5,	1e-6, 0.01;...
              1e-5,	1e-4, 1];
%%% OUTPUT %%%
mrw.P = [];             % Parameters per iteration
mrw.L = [];             % Log-likelihood per iteration
mrw.e = [];             % Maximum error in the distribution estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% (5.1) Define random number generator:
rng(mrw.s,'twister');
r.p0 = rand(5,1);
for i = 1:size(myI,1)
    Z(i) = mrw.(myI{i,1}).(myI{i,2})(myI{i,3},1);
end
r.Pe = mvnrnd(zeros(mrw.I,length(Z)),Z);
r.tL = rand(mrw.I,1);
clear i Z
    
% (5.2) Initialize system:
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

% (5.3) Iterate:
for j = 2:mrw.I
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
        [m,l,e] = SIM_logL(x,N,maxM,Par,ParS);
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
        save(cat(2,'TEMP_MRW_',myGeneModel,'(',psN,')',...
            '(s',num2str(mrw.s),').mat'),...
            'mrw','j','r','Par','ParS');
    end
end
clear j myP myL myE

% (5.4) Save:
save(cat(2,'MRW_',myGeneModel,'(',psN,')',...
    '(s',num2str(mrw.s),').mat'),...
    'mrw','Par','ParS');
delete(cat(2,'TEMP_MRW_',myGeneModel,'(',psN,')',...
            '(s',num2str(mrw.s),').mat'));
```

### (6) Explore sampling effect (i.e. synthetic data):
(6.1) Define probability distribution to sample from.
(6.2) Generate synthetic data.
(6.3) Run Metropolis Random Walk algorithm.

```matlab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Synthetic random sampling parameters:
N = '2S';                          % Model ('2S' or '3S').
maxM = 300;                        % Maximum mRNA number to consider.
srs.t = {'t00','t05','t15','t25'}; % Time points.
srs.n = [173 174 122 115];         % Sample size; if 0, use x sample size.
srs.s = 1;                         % Random number generator seed.
% Model kinetic parameters:
Par.kON = NaN;      % Promoter activation rate (OFF->ON)
Par.kOFF = 0.0459;	% Promoter deactivation rate (ON->OFF)
Par.kONs = NaN;     % Promoter "super" activation rate (ON->ONs)
Par.kOFFs = NaN;	% Promoter "super" deactivation rate (ONs->ON)
Par.mu0 = 0.0260;	% mRNA synthesis rate of promoter in OFF state
Par.mu = 5.65;	    % mRNA synthesis rate of promoter in ON state
Par.muS = NaN;   	% mRNA synthesis rate of promoter in ONs state
Par.d = 0.0559;     % mRNA degradation rate
% Kinetic parameters sensitive to stimulus:
ParS.kON = [0.00492,0.1110]; % [Basal,Stimulus]

% Metropolis Random Walk algorithm parameters:
mrw.I = 100000;     % Number of iterations
mrw.s = 1;          % Random number generator seed
%%% mrw.([U/S]).(par) = [sigma,min,max] %%%
mrw.U.kOFF = [1e-5,	1e-5, 0.1];
mrw.U.mu0  = [1e-5,	1e-5, 0.1];
mrw.U.mu   = [1e-3,	1e-3, 10];
mrw.S.kON  = [1e-5,	1e-6, 0.01;...
              1e-5,	1e-4, 1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% (6.1) Define probability distribution to sample from:
clear x
for t = 1:length(srs.t)
    x.(srs.t{t}) = DATA_Xm(zeros(2,3),N,maxM,0)*0;
end
[M,l,e] = SIM_logL_Opt(x,N,maxM,Par,ParS);
clear t l e

% (6.2) Generate synthetic data:
for t = 1:length(srs.t)
    x.(srs.t{t}) = DATA_Xrs(M.(srs.t{t}),srs.n(t),(srs.s*10)+t);
end
clear t
% Save:
ps = fieldnames(ParS);
psN = '';
for i = 1:length(ps)
    psN = cat(2,psN,ps{i});
end
clear ps i
save(cat(2,'myData_RS(',N,',',num2str(maxM),',',psN,')',...
    '(n',num2str(srs.n(1)),'s',num2str(srs.s),').mat'),...
    'N','maxM','srs','Par','ParS','x')

% (6.3) Run Metropolis Random Walk algorithm (multiple chains):
myGeneModel = cat(2,'RS(',N,',',num2str(maxM),',',psN,')',...
    '(n',num2str(srs.n(1)),'s',num2str(srs.s),')');
SIM_mrw(myGeneModel,Par,mrw,0)
```

## Referencing

If you use this code or the data associated with it please cite:

Gómez-Schiavon *et al.* (2017); https://doi.org/10.1101/109603.

## Copyright

(C) Copyright 2017 Mariana Gómez-Schiavon

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
