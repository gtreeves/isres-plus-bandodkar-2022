% script_isres
clear
clc
close all 



%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
% CONTROLS
%
modelname           = 'schmierer';

% evo options
evo.G               = 50;        % maximum number of generations
evo.lambda          = 150;          % population size (number of offspring) (100 to 400) 
evo.mue             = 20;          % ~round((lambda)/7)
evo.nIslands        = 1;           % number of islands
evo.migGen          = -1;           % ~10 times migration is executed

evo.alphaa          = 0.2;         % Smoothing factor
evo.gammaa          = 0.85;        % Recombination parameter
evo.pf              = .45;         % pressure on fitness in [0 0.5] try 0.45
evo.varphi          = 1;           % expected rate of convergence (usually 1)
evo.tmax            = 72*3600;     % hrs * seconds/hr
evo.mm              = 'min';       % Minimize ('min') or Maximize('max')

% plus options
plus.yeslin         = true;
plus.yesnew         = true;
plus.nlin           = 2;        % 1-2       % 2  
plus.nnewt          = 1;        % 1-2       % 1
plus.nlinPar        = 1;        % 1-2       % 1
plus.nnewtPar       = 1;        % 1-2       % 1
plus.startlin       = 1;
plus.endlin         = evo.G;
plus.startnew       = 1;
plus.endnew         = evo.G;
plus.useFullNewtonStep      = true;     % applies only to newton step.
plus.sortPrevParamsByError  = false;     % applies to both linstep and newton step
plus.betalin = 1;

% model options
model.modelname     =  modelname;

% other options
options.reshuffGen  = 0; 
options.liveUpdates = true;
options.restart     = false;
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


% Create options super-structure
options.evo         = evo;
options.plus        = plus;
%options.model       = model;


%
% Run isres
%
fhandle  = @error_calc_RECI;

lu       = [[0.0001 0.00001 0.0001 0.0001 0.00001 0.0001 0.00001 0.00001 0.0001 0.0001];[10 10 10 10 10 100 10 10 10 1000]];
lu       = log10(lu);

[xb,BestMin,Gm,Stats,opts] = isres_plus(fhandle,lu,options);
xb = 10.^xb;
BestMin;
Gm;
Stats;
opts;

filename = ['./Results_isres-plus_RECI',datestr(now,'yyyy-mm-dd-HH-MM-SS'),'-',char(randi([97,106],1,2)),num2str(randi(10))];
save([filename,'.mat'])   




