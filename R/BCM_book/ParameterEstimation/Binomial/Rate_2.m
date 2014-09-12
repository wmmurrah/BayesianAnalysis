%% Difference Between Two Rates

clear;

sampler = 0; % Choose 0=WinBUGS, 1=JAGS

%% Data (Observed Variables)
k1 = 0;
n1 = 5;
k2 = 5;
n2 = 10;

%% Sampling
% MCMC Parameters
nchains = 3; % How Many Chains?
nburnin = 0; % How Many Burn-in Samples?
nsamples = 5e4;  %How Many Recorded Samples?
nthin = 1; % How Often is a Sample Recorded?
doparallel = 0; % Parallel Option

% Assign Matlab Variables to the Observed Nodes
datastruct = struct('k1',k1,'n1',n1,'k2',k2,'n2',n2);

%Initialize Unobserved Variables
for i=1:nchains
    S.theta1 = 0.5; % An Intial Value for the Success Rate
    S.theta2 = 0.5; % An Intial Value for the Success Rate
    init0(i) = S;
end

if ~sampler
    % Use WinBUGS to Sample
tic
[samples, stats] = matbugs(datastruct, ...
    fullfile(pwd, 'Rate_2.txt'), ...
    'init', init0, ...
    'nChains', nchains, ...
    'view', 1, 'nburnin', nburnin, 'nsamples', nsamples, ...
    'thin', nthin, 'DICstatus', 0, 'refreshrate',100, ...
    'monitorParams', {'theta1','theta2','delta'}, ...
    'Bugdir', 'C:/Program Files/WinBUGS14');
    toc
else
    % Use JAGS to Sample
    tic
    fprintf( 'Running JAGS ...\n' );
    [samples, stats] = matjags( ...
        datastruct, ...
        fullfile(pwd, 'Rate_2.txt'), ...
        init0, ...
        'doparallel' , doparallel, ...
        'nchains', nchains,...
        'nburnin', nburnin,...
        'nsamples', nsamples, ...
        'thin', nthin, ...
        'monitorparams', {'theta1','theta2','delta'}, ...
        'savejagsoutput' , 1 , ...
        'verbosity' , 1 , ...
        'cleanup' , 0 , ...
        'workingdir' , 'tmpjags' );
    toc
end;

%% Analysis
% Plot Posterior
figure(2);clf;hold on;
eps = .025; binsc = -1+eps/2:eps:1-eps/2; binse = -1:eps:1;
count = histc(reshape(samples.delta,1,[]),binse);
count = count(1:end-1);
count = count/sum(count)/eps;
ph = plot(binsc,count,'k-');
set(gca,'box','on','fontsize',14,'xtick',[-1:.2:1],'ytick',[1:ceil(max(get(gca,'ylim')))]);
xlabel('Difference in Rates','fontsize',16);
ylabel('Posterior Density','fontsize',16);

% Summaries of Posterior
% MEAN
disp(sprintf('Mean is %1.2f',stats.mean.delta));
% MODE
[~,ind] = max(count);
disp(sprintf('Mode is %1.2f',binsc(ind)));
% MEDIAN
disp(sprintf('Median is %1.2f',median(reshape(samples.delta,1,[]))));
% CREDIBLE INTERVAL
cred = 0.95;
b1 = (1-cred)/2;b2=1-b1;
val = sort(reshape(samples.delta,1,[]));
disp(sprintf('%d percent credible interval is [%1.2f, %1.2f]',cred*100,val(round(b1*nsamples*nchains)),val(round(b2*nsamples*nchains))));



