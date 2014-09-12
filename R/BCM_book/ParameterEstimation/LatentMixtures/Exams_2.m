%% Exam Scores With Individual Differences

clear;

sampler = 0; % Choose 0=WinBUGS, 1=JAGS

%% Data
k = [21 17 21 18 22 31 31 34 34 35 35 36 39 36 35];n=40;

% Constants
p=length(k); % Number of people


%% Sampling
% MCMC Parameters
nchains = 2; % How Many Chains?
nburnin = 1e3; % How Many Burn-in Samples?
nsamples = 1e3;  %How Many Recorded Samples?
nthin = 1; % How Often is a Sample Recorded?
doparallel = 0; % Parallel Option

% Assign Matlab Variables to the Observed Nodes
datastruct = struct('p',p,'k',k,'n',n);

% Initial Values to Supply to WinBugs
for i=1:nchains
    S.z = round(rand(1,p)); % Initial Random Group Assignments
    S.mu = 0.75;
    S.lambda = 10;
    init0(i) = S;
end

if ~sampler
    % Use WinBUGS to Sample
    tic
    [samples, stats] = matbugs(datastruct, ...
        fullfile(pwd, 'Exams_2.txt'), ...
        'init', init0, ...
        'nChains', nchains, ...
        'view', 1, 'nburnin', nburnin, 'nsamples', nsamples, ...
        'thin', nthin, 'DICstatus', 0, 'refreshrate',100, ...
        'monitorParams', {'predphi','theta','z','mu','sigma'}, ...
        'Bugdir', 'C:/Program Files/WinBUGS14');
else
    % Use JAGS to Sample
    tic
    fprintf( 'Running JAGS with chains serially...\n' );
    [samples, stats] = matjags( ...
        datastruct, ...
        fullfile(pwd, 'Exams_2J.txt'), ...
        init0, ...
        'doparallel' , doparallel, ...
        'nchains', nchains,...
        'nburnin', nburnin,...
        'nsamples', nsamples, ...
        'thin', nthin, ...
        'monitorparams', {'predphi','theta','z','mu','sigma'}, ...
        'savejagsoutput' , 1 , ...
        'verbosity' , 1 , ...
        'cleanup' , 0 , ...
        'workingdir' , 'tmpjags' );
    toc
end;