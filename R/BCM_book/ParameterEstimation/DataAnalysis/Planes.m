%% Planes
   
clear;

sampler = 0; % Choose 0=WinBUGS, 1=JAGS

%% Data
x = 10;
k = 4;
n = 5;
tmax = 50;
        
%% Sampling
% MCMC Parameters
nchains = 2; % How Many Chains?
nburnin = 1e3; % How Many Burn-in Samples?
nsamples = 2e3;  %How Many Recorded Samples?
nthin = 1; % How Often is a Sample Recorded?
doparallel = 0; % Parallel Option

% Assign Matlab Variables to the Observed Nodes
datastruct = struct('k',k,'x',x,'n',n,'tmax',tmax);

% Initialize Unobserved Variables
for i=1:nchains
    S.t = tmax;
    init0(i) = S;
end

if ~sampler
    % Use WinBUGS to Sample
    tic
[samples, stats] = matbugs(datastruct, ...
    fullfile(pwd, 'Planes.txt'), ...
    'init', init0, ...
    'nChains', nchains, ...
    'view', 1, 'nburnin', nburnin, 'nsamples', nsamples, ...
    'thin', nthin, 'DICstatus', 0, 'refreshrate',100, ...
    'monitorParams', {'t'}, ...
    'Bugdir', 'C:/Program Files/WinBUGS14');
else
    % Use JAGS to Sample
    tic
    fprintf( 'Running JAGS ...\n' );
    [samples, stats] = matjags( ...
        datastruct, ...
         fullfile(pwd, 'PlanesJ.txt'), ...
        init0, ...
        'doparallel' , doparallel, ...
        'nchains', nchains,...
        'nburnin', nburnin,...
        'nsamples', nsamples, ...
        'thin', nthin, ...
        'monitorparams', {'t'}, ...
        'savejagsoutput' , 1 , ...
        'verbosity' , 1 , ...
        'cleanup' , 0 , ...
        'workingdir' , 'tmpjags' );
    toc
end;

%% Analysis
figure(1);clf;hold on;
set(gcf,'units','norm','pos',[.2 .2 .6 .4],'paperpositionmode','auto');
bins = x+n-k:tmax;
count = hist(reshape(samples.t,1,[]),bins);
count = count/sum(count);
ph = bar(bins,count);
set(ph,'facecolor','k');
set(gca,'box','on','xlim',[x+n-k-1 tmax+1],'xtick',[0 min(samples.t(:)) tmax],'fontsize',14);
xlabel('Number of Planes','fontsize',16);
ylabel('Posterior Mass','fontsize',16);



