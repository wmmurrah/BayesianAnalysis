%% Prior and Posterior Predictive, Second Example

clear;

sampler = 0; % Choose 0=WinBUGS, 1=JAGS

%% Data
k1 = 0;
n1 = 10;
k2 = 10;
n2 = 10;

%% Sampling
% MCMC Parameters
nchains = 2; % How Many Chains?
nburnin = 1e2; % How Many Burn-in Samples?
nsamples = 1e4;  %How Many Recorded Samples?
nthin = 1; % How Often is a Sample Recorded?
doparallel = 0; % Parallel Option

% Assign Matlab Variables to the Observed Nodes
datastruct = struct('k1',k1,'k2',k2,'n1',n1,'n2',n2);

% Initialize Unobserved Variables
for i=1:nchains
    S.theta = 0.5; % Intial Value
    init0(i) = S;
end

if ~sampler
    % Use WinBUGS to Sample
    tic
    [samples, stats] = matbugs(datastruct, ...
        fullfile(pwd, 'Rate_5.txt'), ...
        'init', init0, ...
        'nChains', nchains, ...
        'view', 1, 'nburnin', nburnin, 'nsamples', nsamples, ...
        'thin', nthin, 'DICstatus', 0, 'refreshrate',100, ...
        'monitorParams', {'theta','postpredk1','postpredk2'}, ...
        'Bugdir', 'C:/Program Files/WinBUGS14');
    toc
else
    % Use JAGS to Sample
    tic
    fprintf( 'Running JAGS ...\n' );
    [samples, stats] = matjags( ...
        datastruct, ...
        fullfile(pwd, 'Rate_5.txt'), ...
        init0, ...
        'doparallel' , doparallel, ...
        'nchains', nchains,...
        'nburnin', nburnin,...
        'nsamples', nsamples, ...
        'thin', nthin, ...
        'monitorparams', {'theta','postpredk1','postpredk2'}, ...
        'savejagsoutput' , 1 , ...
        'verbosity' , 1 , ...
        'cleanup' , 0 , ...
        'workingdir' , 'tmpjags' );
    toc
end;

%% Analysis
figure(5);clf;
set(gcf,'units','norm','pos',[.2 .2 .6 .4],'paperpositionmode','auto');

% Parameter Space
subplot(121);hold on;
axis square;
eps=.015;bins=[0:eps:1];binsc=[eps/2:eps:1-eps/2];
count=histc(reshape(samples.theta,1,[]),bins);
count=count(1:end-1);
count=count/sum(count)/eps;
ph=plot(binsc,count,'k-');
set(gca,'xlim',[0 1],'box','on','fontsize',14,'xtick',[0:.2:1],'ytick',[1:ceil(max(get(gca,'ylim')))]);
set(gca,'box','on','fontsize',14);
xlabel('Rate','fontsize',16);
ylabel('Density','fontsize',16);

% Data Space
subplot(122);hold on;
axis equal;
axis([-1 n1+1 -1 n2+1]);
sc=70;
for i=0:n1
    for j=0:n2
        match=length(find(samples.postpredk1==i&samples.postpredk2==j))/nsamples/nchains;
        if match>0
            ph=plot(i,j,'ks');
        set(ph,'markersize',sc*sqrt(match));
        end;
    end;
end;
ph=plot(k1,k2,'kx');
set(ph,'markersize',16,'linewidth',4);
set(gca,'box','on','fontsize',14,'xtick',[0:n1],'ytick',[0:n2]);
xlabel('Success Count 1','fontsize',16);
ylabel('Success Count 2','fontsize',16);