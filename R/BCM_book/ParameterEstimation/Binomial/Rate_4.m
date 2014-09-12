%% Prior and Posterior Predictive

clear;

sampler = 0; % Choose 0=WinBUGS, 1=JAGS

%% Data
dataset = 1;
switch dataset
    case 1, k = 1; n = 15; % Toy data
    case 2, k = 24; n = 121; % Trouw nursing Home Data
end;

%% Sampling
% MCMC Parameters
nchains = 2; % How Many Chains?
nburnin = 1e2; % How Many Burn-in Samples?
nsamples = 1e4;  %How Many Recorded Samples?
nthin = 1; % How Often is a Sample Recorded?
doparallel = 0; % Parallel Option


% Assign Matlab Variables to the Observed Nodes
datastruct = struct('k',k,'n',n);

% Initialize Unobserved Variables
for i=1:nchains
    S.theta = 0.5; % Intial Value
    init0(i) = S;
end

if ~sampler
    % Use WinBUGS to Sample
    tic
    [samples, stats] = matbugs(datastruct, ...
        fullfile(pwd, 'Rate_4.txt'), ...
        'init', init0, ...
        'nChains', nchains, ...
        'view', 1, 'nburnin', nburnin, 'nsamples', nsamples, ...
        'thin', nthin, 'DICstatus', 0, 'refreshrate',100, ...
        'monitorParams', {'theta','thetaprior','postpredk','priorpredk'}, ...
        'Bugdir', 'C:/Program Files/WinBUGS14');
    toc
else
    % Use JAGS to Sample
    tic
    fprintf( 'Running JAGS ...\n' );
    [samples, stats] = matjags( ...
        datastruct, ...
        fullfile(pwd, 'Rate_4.txt'), ...
        init0, ...
        'doparallel' , doparallel, ...
        'nchains', nchains,...
        'nburnin', nburnin,...
        'nsamples', nsamples, ...
        'thin', nthin, ...
        'monitorparams', {'theta','thetaprior','postpredk','priorpredk'}, ...
        'savejagsoutput' , 1 , ...
        'verbosity' , 1 , ...
        'cleanup' , 0 , ...
        'workingdir' , 'tmpjags' );
    toc
end;

%% Analysis
figure(4);clf;

% Parameter Space
subplot(211);hold on;
eps=.015;binse=[0:eps:1];binsc=[eps/2:eps:1-eps/2];
count=histc(reshape(samples.thetaprior,1,[]),binse);
count=count(1:end-1);
count=count/sum(count)/eps;
ph=plot(binsc,count,'k--');
count=histc(reshape(samples.theta,1,[]),binse);
count=count(1:end-1);
count=count/sum(count)/eps;
ph=plot(binsc,count,'k-');
set(gca,'xlim',[0 1],'box','on','fontsize',14,'xtick',[0:.2:1]);
legend('Prior','Posterior');
set(gca,'box','on','fontsize',14);
xlabel('Rate','fontsize',16);
ylabel('Density','fontsize',16);

% Data Space
subplot(212);hold on;
kbins=[0:n];
count1=histc(samples.priorpredk,kbins);
count1=count1/sum(count1);
count2=histc(samples.postpredk,kbins);
count2=count2/sum(count2);
switch dataset
    case 1,
        ph=bar(kbins,count1,.6);set(ph,'facecolor','none','linewidth',1.5,'linestyle','--')
        ph=bar(kbins,count2,.8);set(ph,'facecolor','none')
        set(gca,'xlim',[-1 n+1],'box','on','fontsize',14,'xtick',[0:n]);
    otherwise,
          ph=bar(kbins,count1,.6);set(ph,'facecolor','none','linewidth',0.5,'linestyle',':')
        ph=bar(kbins,count2);set(ph,'facecolor','none');
        set(gca,'xlim',[-1 n+1],'box','on','fontsize',14,'xtick',[0 n],'ticklength',[0 0]);
end;
[lh oh] = legend('Prior','Posterior');
set(lh,'box','off');
xlabel('Success Count','fontsize',16);
ylabel('Mass','fontsize',16);