%% Kappa Coefficient of Agreement

clear;

sampler = 0; % Choose 0=WinBUGS, 1=JAGS

%% Data
% Choose dataset
dataset = 1;

switch dataset
    case 1, y = [14 4 5 210]; % Influenza
    case 2, y = [20 7 103 417]; % Hearing Loss
    case 3, y = [0 0 13 157]; % Rare Disease
end;

% Constants
n = sum(y);

%% Sampling
% MCMC Parameters
nchains = 2; % How Many Chains?
nburnin = 5e3; % How Many Burn-in Samples?
nsamples = 5e3;  %How Many Recorded Samples?
nthin = 1; % How Often is a Sample Recorded?
doparallel = 0; % Parallel Option

% Assign Matlab Variables to the Observed Nodes
datastruct = struct('y',y,'n',n);

% Initialize Unobserved Variables
for i=1:nchains
    S.alpha = 0.5;
    S.beta = 0.5;
    S.gamma = 0.5;
    init0(i) = S;
end

if ~sampler
    % Use WinBUGS to Sample
    tic
[samples, stats] = matbugs(datastruct, ...
    fullfile(pwd, 'Kappa.txt'), ...
    'init', init0, ...
    'nChains', nchains, ...
    'view', 1, 'nburnin', nburnin, 'nsamples', nsamples, ...
    'thin', nthin, 'DICstatus', 0, 'refreshrate',100, ...
    'monitorParams', {'kappa','xi','psi','alpha','beta','gamma','pi'}, ...
    'Bugdir', 'C:/Program Files/WinBUGS14');
else
    % Use JAGS to Sample
    tic
    fprintf( 'Running JAGS ...\n' );
    [samples, stats] = matjags( ...
        datastruct, ...
        fullfile(pwd, 'Kappa.txt'), ...
        init0, ...
        'doparallel' , doparallel, ...
        'nchains', nchains,...
        'nburnin', nburnin,...
        'nsamples', nsamples, ...
        'thin', nthin, ...
        'monitorparams', {'kappa','xi','psi','alpha','beta','gamma','pi'}, ...
        'savejagsoutput' , 1 , ...
        'verbosity' , 1 , ...
        'cleanup' , 0 , ...
        'workingdir' , 'tmpjags' );
    toc
end;