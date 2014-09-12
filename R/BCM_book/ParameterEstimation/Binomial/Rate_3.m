%% Inferring A Common Rate

clear;

sampler = 0; % Choose 0=WinBUGS, 1=JAGS

%% Data
k1 = 5;
n1 = 10;
k2 = 7;
n2 = 10;

%% Sampling
% MCMC Parameters
nchains = 2; % How Many Chains?
nburnin = 0; % How Many Burn-in Samples?
nsamples = 1e4;  %How Many Recorded Samples?
nthin = 1; % How Often is a Sample Recorded?
doparallel = 0; % Parallel Option

% Assign Matlab Variables to the Observed Nodes
datastruct = struct('k1',k1,'n1',n1,'k2',k2,'n2',n2);

% Initialize Unobserved Variables
for i=1:nchains
    S.theta = 0.5; % An Intial Value for the Success Rate
    init0(i) = S;
end

if ~sampler
    % Use WinBUGS to Sample
tic
% Use WinBUGS to Sample
[samples, stats] = matbugs(datastruct, ...
    fullfile(pwd, 'Rate_3.txt'), ...
    'init', init0, ...
    'nChains', nchains, ...
    'view', 1, 'nburnin', nburnin, 'nsamples', nsamples, ...
    'thin', nthin, 'DICstatus', 0, 'refreshrate',100, ...
    'monitorParams', {'theta'}, ...
    'Bugdir', 'C:/Program Files/WinBUGS14');
    toc
else
    % Use JAGS to Sample
    tic
    fprintf( 'Running JAGS ...\n' );
    [samples, stats] = matjags( ...
        datastruct, ...
        fullfile(pwd, 'Rate_3.txt'), ...
        init0, ...
        'doparallel' , doparallel, ...
        'nchains', nchains,...
        'nburnin', nburnin,...
        'nsamples', nsamples, ...
        'thin', nthin, ...
        'monitorparams', {'theta'}, ...
        'savejagsoutput' , 1 , ...
        'verbosity' , 1 , ...
        'cleanup' , 0 , ...
        'workingdir' , 'tmpjags' );
    toc
end;

%% Analysis
figure(3);clf;hold on;
eps = .01; binsc = eps/2:eps:1-eps/2; binse = 0:eps:1;
count = histc(reshape(samples.theta,1,[]),binse);
count = count(1:end-1);
count = count/sum(count)/eps;
ph = plot(binsc,count,'k-');
set(gca,'box','on','fontsize',14,'xtick',[0:.2:1],'ytick',[1:ceil(max(get(gca,'ylim')))]);
xlabel('Rate','fontsize',16);
ylabel('Posterior Density','fontsize',16);
