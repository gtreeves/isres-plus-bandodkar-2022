% script_isres_plus
clear
clc
close all 

addpath(genpath('./'))

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
% CONTROLS
%
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% model options
modelname     = 'dlcact';

% % evo options
evo.G               = 100;        % maximum number of generations
evo.lambda          = 150;          % population size (number of offspring) (100 to 400) 
evo.mue             = 20;          % ~round((lambda)/7)

% plus options
plus.yeslin         = true;
plus.yesnew         = true;
plus.nlin           = 2;        % 1-2      
plus.nnewt          = 1;        % 1-2      
% plus.startlin       = 1;
% plus.endlin         = 100;
% plus.startnew       = 1;
% plus.endnew         = 100;

% other options
options.liveUpdates = true;
model.modelname     = modelname;

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


%
% Lower and upper bounds specifications for the models. Add any parameters
% that do not need to be optimized to the addParams variable. Also specify
% the function handle for each model.
%
if strcmp(modelname,'dlcact')
    delta       = 1e-100;                   % neg fbk params 1
    K           = 0.05;                     % neg fbk params 2
    nParams     = 16;                       % 16 free parameters
    lu          = [-4*ones(1,nParams); 4*ones(1,nParams)];
    lu(1,14)    = -1;
    lu(2,14)    = 0;                        % Limits on phi
    addParams   = [log10(delta), log10(K)];
    fhandle     = @calc_error_dc;
elseif strcmp(modelname,'gap')
    lu(1,:)     = [repmat(10,1,4),zeros(1,4),repmat(0.01,1,4),repmat(-1,1,32)];
    lu(2,:)     = [repmat(100,1,4),repmat(0.5,1,4),repmat(0.1,1,4),ones(1,32)];
    lu(:,1:4)   = log10(lu(:,1:4));      % shift bounds to ensure that all values are of order 0
    lu(:,9:12)  = log10(lu(:,9:12));
    addParams   = [];
    fhandle     = @calc_error_gap;
elseif strcmp(modelname,'smad')
    lu          = [[0.0001 0.00001 0.0001 0.0001 0.00001 0.0001 0.00001 0.00001 0.0001 0.0001 10 10 10 1 0.001];[10 10 10 10 10 100 10 10 10 1000 1000 1000 1000 1000 1000]];
    lu          = log10(lu); 
    addParams   = [];
    fhandle     = @calc_error_smad;
else
    error('Please specify correct modelname')
end


% Add additional parameters to model structure
model.addParams     = addParams;


% Create options super-structure
options.evo         = evo;
options.plus        = plus;
options.model       = model;



%
% Run isres plus
%
[xb,BestMin,Gm,Stats,options] = isres_plus(fhandle,lu,options);



% Save Results
if ~exist('Results','dir')
    mkdir('Results')
end
filename = ['Results/Results_isres-plus_',modelname,'_',datestr(now,'yyyy-mm-dd-HH-MM-SS')];
save([filename,'.mat'])   




