%% Two Country Quiz

clear;

sampler = 0; % Choose 0=WinBUGS, 1=JAGS

%% Data
dataset = 1; %Choose Dataset
switch dataset
    case 1, k = [1 0 0 1 1 0 0 1;
            1 0 0 1 1 0 0 1;
            0 1 1 0 0 1 0 0;
            0 1 1 0 0 1 1 0;
            1 0 0 1 1 0 0 1;
            0 0 0 1 1 0 0 1;
            0 1 0 0 0 1 1 0;
            0 1 1 1 0 1 1 0];
    case 2,k = [1 0 0 1 1 0 0 1;
            1 0 0 1 1 0 0 1;
            0 1 1 0 0 1 0 0;
            0 1 1 0 0 1 1 0;
            1 0 0 1 1 0 0 1;
            0 0 0 1 1 0 0 1;
            0 1 0 0 0 1 1 0;
            0 1 1 1 0 1 1 0;
            1 0 0 1 nan nan nan nan;
            0 nan nan nan nan nan nan nan;
            nan nan nan nan nan nan nan nan];
    case 3,k = [1 0 0 1 1 0 0 1;
            1 0 0 1 1 0 0 1;
            1 0 0 1 1 0 0 1;
            1 0 0 1 1 0 0 1;
            1 0 0 1 1 0 0 1;
            1 0 0 1 1 0 0 1;
            1 0 0 1 1 0 0 1;
            1 0 0 1 1 0 0 1;
            1 0 0 1 1 0 0 1;
            1 0 0 1 1 0 0 1;
            1 0 0 1 1 0 0 1;
            1 0 0 1 1 0 0 1;
            0 1 1 0 0 1 0 0;
            0 1 1 0 0 1 1 0;
            1 0 0 1 1 0 0 1;
            0 0 0 1 1 0 0 1;
            0 1 0 0 0 1 1 0;
            0 1 1 1 0 1 1 0;
            1 0 0 1 nan nan nan nan;
            0 nan nan nan nan nan nan nan;
            nan nan nan nan nan nan nan nan];
end;

% Constants
[nx,nz] = size(k); % Number of people and questions

%% Sampling
% MCMC Parameters
nchains = 2; % How Many Chains?
nburnin = 1e3; % How Many Burn-in Samples?
nsamples = 2e3;  %How Many Recorded Samples?
nthin = 1; % How Often is a Sample Recorded?
doparallel = 0; % Parallel Option

% Assign Matlab Variables to the Observed Nodes
datastruct = struct('nx',nx,'nz',nz,'k',k);

% Initial Values to Supply to WinBugs
for i=1:nchains
    S.x = round(rand(1,nx));
    S.z = round(rand(1,nz));
    S.alpha = 1/2;
    S.beta = 1/2;
    init0(i) = S;
end

% Use WinBUGS to Sample
switch dataset
    case 1,
        if ~sampler
            % Use WinBUGS to Sample
            tic
            [samples, stats] = matbugs(datastruct, ...
                fullfile(pwd, 'TwoCountryQuiz.txt'), ...
                'init', init0, ...
                'nChains', nchains, ...
                'view', 1, 'nburnin', nburnin, 'nsamples', nsamples, ...
                'thin', nthin, 'DICstatus', 0, 'refreshrate',100, ...
                'monitorParams', {'x','z','alpha','beta'}, ...
                'Bugdir', 'C:/Program Files/WinBUGS14');
            toc
        else
            % Use JAGS to Sample
            tic
            fprintf( 'Running JAGS ...\n' );
            [samples, stats] = matjags( ...
                datastruct, ...
                fullfile(pwd, 'TwoCountryQuiz.txt'), ...
                init0, ...
                'doparallel' , doparallel, ...
                'nchains', nchains,...
                'nburnin', nburnin,...
                'nsamples', nsamples, ...
                'thin', nthin, ...
                'monitorparams', {'x','z','alpha','beta'}, ...
                'savejagsoutput' , 1 , ...
                'verbosity' , 1 , ...
                'cleanup' , 0 , ...
                'workingdir' , 'tmpjags' );
            toc
        end;
    case {2,3},
        if ~sampler
            % Use WinBUGS to Sample
            tic
            [samples, stats] = matbugs(datastruct, ...
                fullfile(pwd, 'TwoCountryQuiz.txt'), ...
                'init', init0, ...
                'nChains', nchains, ...
                'view', 1, 'nburnin', nburnin, 'nsamples', nsamples, ...
                'thin', nthin, 'DICstatus', 0, 'refreshrate',100, ...
                'monitorParams', {'x','z','alpha','beta','k'}, ...
                'Bugdir', 'C:/Program Files/WinBUGS14');
            toc
        else
            % Use JAGS to Sample
            tic
            fprintf( 'Running JAGS ...\n' );
            [samples, stats] = matjags( ...
                datastruct, ...
                fullfile(pwd, 'TwoCountryQuiz.txt'), ...
                init0, ...
                'doparallel' , doparallel, ...
                'nchains', nchains,...
                'nburnin', nburnin,...
                'nsamples', nsamples, ...
                'thin', nthin, ...
                'monitorparams', {'x','z','alpha','beta','k'}, ...
                'savejagsoutput' , 1 , ...
                'verbosity' , 1 , ...
                'cleanup' , 0 , ...
                'workingdir' , 'tmpjags' );
            toc
        end;
end;