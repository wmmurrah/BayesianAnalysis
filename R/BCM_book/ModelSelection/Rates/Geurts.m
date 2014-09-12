% Geurts
clear;

sampler = 0; % Choose 0=WinBUGS, 1=JAGS
run_model = 1; % set 0 to load samples, or 1 to run WinBUGS


%% Data
% Control
num_errors = [15,10, 61,11, 60, 44, 63, 70, 57,11, 67, 21, 89,12, 63, 11, 96,10, 37,19, 44,18, 78, 27, 60,14];
nc = [89,74,128,87,128,121,128,128,128,78,128,106,128,83,128,100,128,73,128,86,128,86,128,100,128,79];
kc = nc-num_errors;
nsc = length(kc);
% ADHD
num_errors = [88, 50, 58,17, 40, 18,21, 50, 21, 69,19, 29,11, 76, 46, 36, 37, 72,27, 92,13, 39, 53, 31, 49, ...
    57,17,10,12,21, 39, 43, 49,17, 39,13, 68, 24, 21,27, 48, 54, 41, 75, 38, 76,21, 41, 61,24, 28,21];
na =[128,128,128,86,128,117,89,128,110,128,93,107,87,128,128,113,128,128,98,128,93,116,128,116,128, ...
    128,93,86,86,96,128,128,128,86,128,78,128,111,100,95,128,128,128,128,128,128,98,127,128,93,110,96];
ka = na-num_errors;
nsa = length(ka);

%% Sampling
% MCMC Parameters
nchains = 3; % How Many Chains?
nburnin = 1e4; % How Many Burn-in Samples?
nsamples = 1e5;  %How Many Recorded Samples?
nthin = 1; % How Often is a Sample Recorded?
doparallel = 0; % Parallel Option

% Assign Matlab Variables to the Observed WinBUGS Nodes
datastruct = struct('kc',kc,'ka',ka,'nc',nc,'na',na,'nsc',nsc,'nsa',nsa);

% Initialize Unobserved Variables
for i=1:nchains
    S.mu = randn;
    S.sigma = 1;
    S.delta = randn;
    init0(i) = S;
end

if ~run_model
    load Geurts samples stats
else
    if ~ sampler
        % Use WinBUGS to Sample
        [samples, stats] = matbugs(datastruct, ...
            fullfile(pwd, 'Geurts.txt'), ...
            'init', init0, ...
            'nChains', nchains, ...
            'view', 1, 'nburnin', nburnin, 'nsamples', nsamples, ...
            'thin', nthin, 'DICstatus', 0, 'refreshrate',100, ...
            'monitorParams', {'delta','deltaprior'}, ...
            'Bugdir', 'C:/Program Files/WinBUGS14');
    else
        % Use JAGS to Sample
        tic
        fprintf( 'Running JAGS ...\n' );
        [samples, stats] = matjags( ...
            datastruct, ...
            fullfile(pwd, 'Geurts.txt'), ...
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
    save Geurts samples stats
end;

%% Analysis
crit = 0;
eps =.05;mineffect=-3;maxeffect = 3;
binsc = [mineffect+eps/2:eps:maxeffect-eps/2];binse = [mineffect:eps:maxeffect];

figure(1);clf;hold on;
ff=1.3;set(gcf,'units','norm','pos',[.2 .2 .4*ff .4*ff],'paperpositionmode','auto','color','w');

ph = plot([-100 -99],[-100 -100],'k--');
set(ph,'linewidth',1.5);
ph = plot([-100 -99],[-100 -100],'k-');
[lh oh] = legend('Posterior','Prior','location','northeast');
set(ph,'linewidth',1.5);
set(lh,'box','off');

% Posterior on Delta
if exist('aksdensity')
    [f,xi] = ksdensity(reshape(samples.delta,1,[]),'kernel','normal','npoints',5e3);
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
if exist('aksdensity')
    tmp = reshape(samples.deltaprior,1,[]);
    tmp = tmp(find(tmp>binse(1)&tmp<binse(end)));
    [f2,xi] = ksdensity(tmp,'kernel','normal','support',[binse(1) binse(end)],'npoints',5e3);
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
set(gca,'xlim',[binse(1) binse(end)],'ytick',[]);
set(gca,'box','on','fontsize',14,'xtick',[binse(1):binse(end)],'ticklength',[0 0]);
xlabel('Delta','fontsize',16);
ylabel('Density','fontsize',16);


if exist('aksdensity')
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

