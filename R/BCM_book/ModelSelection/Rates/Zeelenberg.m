% Zeelenberg
clear;

sampler = 0; % Choose 0=WinBUGS, 1=JAGS
run_model = 1; % set 0 to load samples, or 1 to run WinBUGS


%%  Data
% Study Both
sb = [15,11,15,14,15,18,16,16,18,16,15,13,18,12,11,13,17,18,16,11,17,18,12,18,18,14,21,18,17,10, ...
    11,12,16,18,17,15,19,12,21,15,16,20,15,19,16,16,14,18,16,19,17,11,19,18,16,16,11,19,18,12, ...
    15,18,20, 8,12,19,16,16,16,12,18,17,11,20]; 
nb = 21;
% Study Neither
sn =[15,12,14,15,13,14,10,17,13,16,16,10,15,15,10,14,17,18,19,12,19,18,10,18,16,13,15,20,13,15, ...
    13,14,19,19,19,18,13,12,19,16,14,17,15,16,15,16,13,15,14,19,12,11,17,13,18,13,13,19,18,13, ...
    13,16,18,14,14,17,12,12,16,14,16,18,13,13]; 
nn = 21;

% Constants
ns = length(sb);

%% Sampling
% MCMC Parameters
nchains = 3; % How Many Chains?
nburnin = 1e4; % How Many Burn-in Samples?
nsamples = 2e5;  %How Many Recorded Samples?
nthin = 1; % How Often is a Sample Recorded?
doparallel = 0; % Parallel Option

% Assign Matlab Variables to the Observed WinBUGS Nodes
datastruct = struct('sb',sb,'sn', sn, 'nb', nb, 'nn',nn, 'ns',ns);

% Initialize Unobserved Variables
for i=1:nchains
    S.mu = rand;
    S.sigma = 1;
    S.delta = 1;
    S.sigmaalpha = 1;
    init0(i) = S;
end
         
         
if ~run_model
    load Zeelenberg samples stats
else
    if ~ sampler
        % Use WinBUGS to Sample
        [samples, stats] = matbugs(datastruct, ...
            fullfile(pwd, 'Zeelenberg.txt'), ...
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
            fullfile(pwd, 'Zeelenberg.txt'), ...
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
    save Zeelenberg samples stats
end;

%% Analysis
crit = 0;
eps =.05;maxeffect = 4;
binsc = [eps/2:eps:maxeffect];binse = [0:eps:maxeffect];

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