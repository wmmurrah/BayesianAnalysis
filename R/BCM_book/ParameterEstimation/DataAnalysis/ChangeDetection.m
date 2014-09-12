%% Change Detection

clear;

sampler = 0; % Choose 0=WinBUGS, 1=JAGS

%% Data
load Changepointdata;c = data;
n = length(c);
t = 1:n;
        
%% Sampling
% MCMC Parameters
nchains = 2; % How Many Chains?
nburnin = 1e3; % How Many Burn-in Samples?
nsamples = 1e3;  %How Many Recorded Samples?
nthin = 1; % How Often is a Sample Recorded?
doparallel = 0; % Parallel Option

% Assign Matlab Variables to the Observed Nodes
datastruct = struct('c',c,'t',t,'n',n);

% Initialize Unobserved Variables
for i=1:nchains
    S.mu = ones(1,2); % An Intial Value for both means
    S.lambda = 1; % An Intial Value for the common precision
    S.tau= n/2; % An initial value for the change-point
    init0(i) = S;
end

if ~sampler
    % Use WinBUGS to Sample
    tic
[samples, stats] = matbugs(datastruct, ...
    fullfile(pwd, 'ChangeDetection.txt'), ...
    'init', init0, ...
    'nChains', nchains, ...
    'view', 1, 'nburnin', nburnin, 'nsamples', nsamples, ...
    'thin', nthin, 'DICstatus', 0, 'refreshrate',100, ...
    'monitorParams', {'mu','sigma','tau'}, ...
    'Bugdir', 'C:/Program Files/WinBUGS14');
else
    % Use JAGS to Sample
    tic
    fprintf( 'Running JAGS ...\n' );
    [samples, stats] = matjags( ...
        datastruct, ...
        fullfile(pwd, 'ChangeDetection.txt'), ...
        init0, ...
        'doparallel' , doparallel, ...
        'nchains', nchains,...
        'nburnin', nburnin,...
        'nsamples', nsamples, ...
        'thin', nthin, ...
        'monitorparams', {'mu','sigma','tau'}, ...
        'savejagsoutput' , 1 , ...
        'verbosity' , 1 , ...
        'cleanup' , 0 , ...
        'workingdir' , 'tmpjags' );
    toc
end;

%% Analysis
figure(1);clf;hold on;
ph = plot(1:n,c,'k-');
set(ph,'color',[.5 .5 .5]);
mt = stats.mean.tau;
ph = plot([1 mt],mean(samples.mu(1,:,1))*[1 1],'k-');
set(ph,'linewidth',2);
ph = plot([mt n],mean(samples.mu(1,:,2))*[1 1],'k-');
set(ph,'linewidth',2);
set(gca,'box','on','xlim',[0 n+1],'fontsize',14);
xlabel('Time','fontsize',16);
ylabel('Count','fontsize',16);



