% One Sample Comparison of Means
clear;

sampler = 0; % Choose 0=WinBUGS, 1=JAGS
run_model = 1; % set 0 to load samples, or 1 to run WinBUGS

%% Data
Winter = [-0.05 0.41 0.17 -0.13 0.00 -0.05 0.00 0.17 0.29 0.04 0.21 0.08 0.37 0.17 0.08 -0.04 -0.04 0.04 -0.13 -0.12 0.04 0.21 0.17 0.17 0.17 0.33 0.04 0.04 0.04 0.00 0.21 0.13 0.25 -0.05 0.29 0.42 -0.05 0.12 0.04 0.25 0.12];
Summer = [0.00 0.38 -0.12 0.12 0.25 0.12 0.13 0.37 0.00 0.50 0.00 0.00 -0.13 -0.37 -0.25 -0.12 0.50 0.25 0.13 0.25 0.25 0.38 0.25 0.12 0.00 0.00 0.00 0.00 0.25 0.13 -0.25 -0.38 -0.13 -0.25 0.00 0.00 -0.12 0.25 0.00 0.50 0.00];

% Allowed because it is a within-subjects design
x = Winter - Summer;
x = x/std(x); % Standardize

% Constants
ndata = length(Winter);

%% Sampling
% MCMC Parameters
nchains = 4; % How Many Chains?
nburnin = 1e3; % How Many Burn-in Samples?
nsamples = 1e5;  %How Many Recorded Samples?
nthin = 1; % How Often is a Sample Recorded?
doparallel = 0; % Parallel Option

% Assign Matlab Variables to the Observed WinBUGS Nodes
datastruct = struct('x',x,'ndata',ndata);

% Initialize Unobserved Variables
for i=1:nchains
    S.delta = randn;
    S.sigmatmp = rand;
    init0(i) = S;
end

if ~run_model
    load OneSample samples stats
else
    if ~sampler
        % Use WinBUGS to Sample
        [samples, stats] = matbugs(datastruct, ...
            fullfile(pwd, 'OneSample.txt'), ...
            'init', init0, ...
            'nChains', nchains, ...
            'view', 0, 'nburnin', nburnin, 'nsamples', nsamples, ...
            'thin', nthin, 'DICstatus', 0, 'refreshrate',100, ...
            'monitorParams', {'delta','deltaprior'}, ...
            'Bugdir', 'C:/Program Files/WinBUGS14');
    else
        % Use JAGS to Sample
        tic
        fprintf( 'Running JAGS ...\n' );
        [samples, stats] = matjags( ...
            datastruct, ...
            fullfile(pwd, 'OneSample.txt'), ...
            init0, ...
            'doparallel' , doparallel, ...
            'nchains', nchains,...
            'nburnin', nburnin,...
            'nsamples', nsamples, ...
            'thin', nthin, ...
            'monitorparams',  {'delta','deltaprior'}, ...
            'savejagsoutput' , 1 , ...
            'verbosity' , 1 , ...
            'cleanup' , 0 , ...
            'workingdir' , 'tmpjags' );
        toc
    end;
    save OneSample samples stats
end;

%% Analysis
figure(1);clf;hold on;
set(gcf,'units','norm','pos',[.2 .2 .4 .3],'paperpositionmode','auto');
set(gca,'fontsize',14);
crit = 0;
eps =.01;binsc = [-3+eps/2:eps:3-eps/2];binse = [-3:eps:3];

ph = plot([-100 -99],[-100 -100],'k--');
set(ph,'linewidth',1.5);
ph = plot([-100 -99],[-100 -100],'k-');
[lh oh] = legend('Posterior','Prior','location','northwest');
set(ph,'linewidth',1.5);
set(lh,'box','off');

% Posterior on Delta
if exist('ksdensity')
    [f,xi] = ksdensity(reshape(samples.delta,1,[]),'kernel','normal','support',[binse(1) binse(end)]);
    ph = plot(xi,f,'k--');
    set(ph,'linewidth',1.5);
    [valk indk] = min(abs(xi-crit));
    ph = plot(crit,f(indk),'ko');
    set(ph,'markerfacecolor','w','markersize',8,'linewidth',1.5);
else
    count = histc(reshape(samples.delta,1,[]),binse);
    count = count(1:end-1);
    count = count/sum(count)/eps;
    ph = plot(binsc,count,'k--');
    set(ph,'linewidth',1.5);
    [val ind] = min(abs(binsc-crit));
    ph = plot(crit,count(ind),'ko');
    set(ph,'markerfacecolor','w','markersize',8,'linewidth',1.5);
end;

% Prior on Delta
if exist('ksdensity')
    tmp = reshape(samples.deltaprior,1,[]);
    tmp = tmp(find(tmp>binse(1)&tmp<binse(end)));
    [f2,xi] = ksdensity(tmp,'kernel','normal','support',[binse(1) binse(end)]);
    ph = plot(xi,f2,'k-');
    set(ph,'linewidth',1.5);
    [valk2 indk2] = min(abs(xi-crit));
    ph = plot(crit,f2(indk2),'ko');
    set(ph,'markerfacecolor','k','markersize',8,'linewidth',1.5);
else
    count2 = histc(reshape(samples.deltaprior,1,[]),binse);
    count2 = count2(1:end-1);
    count2 = count2/sum(count2)/eps;
    ph = plot(binsc,count2,'k-');
    set(ph,'linewidth',1.5);
    [val2 ind2] = min(abs(binsc-crit));
    ph = plot(crit,count2(ind2),'ko');
    set(ph,'markerfacecolor','k','markersize',8,'linewidth',1.5);
end;
set(gca,'xlim',[binse(1) binse(end)]);
set(gca,'box','on','fontsize',14,'xtick',[binse(1):binse(end)],'ticklength',[0 0],'ytick',[]);
xlabel('Delta','fontsize',16);
ylabel('Density','fontsize',16);

if exist('ksdensity')
    v1 = f(indk); v2 = f2(indk2);
else
    v1 = count(ind); v2 = count2(ind2);
end;
    bf = [v1/v2 log(v1)-log(v2)]

sc = .9; sc2 = .975;
arrow([crit max(v1,v2)*sc2],[crit min(v1,v2)*1/sc],10,30);
arrow([crit min(v1,v2)*1/sc],[crit max(v1,v2)*sc2],10,30);
str = sprintf(' $%d\\times$',round(max([bf(1) 1/bf(1)])));
th = text(crit,mean([v1 v2]),str);
set(th,'vert','mid','hor','left','fontsize',14,'interp','latex');

