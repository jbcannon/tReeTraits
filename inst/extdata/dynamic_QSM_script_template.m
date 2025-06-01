%------------------------------------------------------------------------->
% USE TREE QSM TO EXPLORE PARAMETERS FOR MODELING BRANCHES
%
% Jeffery B. Cannon
% The Jones Center at Ichuaway
% jeffery.cannon@jonesctr.org
% 31 July 2022
%
% Purpose: The goal of this script is to run TreeQSM (InverseTampere) on 
% individual tree scans (cleaned using pre-processing script). The script
% runs TreeQSM on a range of different parameters, then selects the best
% model by looking at deviance between individual points and fitted
% cylinders. Because we are focusing on bole and primary branches, the 
% fit is evaluated only on these segments.
%------------------------------------------------------------------------->
%------------------------------------------------------------------------->
% SET INPUTS
%------------------------------------------------------------------------->
if ~exist('results/', 'dir')
    mkdir('results/');
end
addpath(genpath('%%TREE_QSM_DIRECTORY%%'));% add path to include all functions
tree_filename = '%%TREE_MAT_PATH%%'; %no .mat extension
tree_name = '%%TREE_ID%%';  % For naming output files

local= parcluster('local');
rmdir(local.JobStorageLocation, 's');

% Set the optimization method based on distances between trunk and 
% primary branches. See select_optimum.m for other settings
opt_method = 'trunk+1branch_mean_dis';

create_input;
inputs.PatchDiam1 = [%%INPUTS_PATCHDIAM1%%];
% These two are the most important to vary per TreeQSM Manual
inputs.PatchDiam2Min = [%%INPUTS_PATCHDIAM2MIN%%];
inputs.PatchDiam2Max = [%%INPUTS_PATCHDIAM2MAX%%];
inputs.lcyl = [%%INPUTS_LCYL%%];

%------------------------------------------------------------------------->
% RUN MODELS AND OPTIMIZATION
%------------------------------------------------------------------------->
% Set remaining inputs automatically (don't need to be adjusted)
% BallRad1 and BallRad2 should be slightly larger than corresponding patch
% diameter per TreeQSM Manual
inputs.BallRad1 = inputs.PatchDiam1 + 0.01;
inputs.BallRad2 = inputs.PatchDiam2Max + 0.01;
inputs.Tria=0;
% Preprocessing removed ground vegetation so OnlyTree=1
inputs.OnlyTree = 1;
inputs.name=tree_name

% Load tree
load([tree_filename,'.mat'])

% Run tons of models
output='tmp';
make_models_parallel(tree_filename, output, %%N_CORES%%, inputs);

% Select the best model
load('results//tmp.mat')
select_optimum(QSMs, opt_method, output)
opt = load(['results//OptimalQSMs_', output, '.mat'])

% Run best model
opt.OptInputs.plot = 0; opt.OptInputs.Tria=0
QSM = treeqsm(tree, opt.OptInputs)

%------------------------------------------------------------------------->
% OUTPUTS
%------------------------------------------------------------------------->
%save results
fn1 = ['%%OUT_DIR%%', tree_name, '_qsm.mat'];
save(fn1,'QSM'); %save outputs
fn2 = ['%%OUT_DIR%%', tree_name, '_qsm_optInputs.xml'];
writestruct(opt.OptInputs, fn2); %save optimum inputs
fn3 = ['%%OUT_DIR%%', tree_name, '_qsm_cylFits.csv'];
writecell(opt.OptDist, fn3); %save 3 best cylinder fits

opt.OptInputs %print optimum inputs
opt.OptDist{1}(1)