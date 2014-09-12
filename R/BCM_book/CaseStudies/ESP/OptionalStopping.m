% Correlation Analysis of Bem Optional Stopping

clear;

sampler = 0; % Choose 0=WinBUGS, 1=JAGS
run_model = 1; % set 0 to load samples, or 1 to run WinBUGS

%% Data
% Sample size n and effect size e in the Bem experiments
x(:,1) = [100, 150, 97, 99, 100, 150, 200, 100, 50];
x(:,2) = [0.25, 0.20, 0.25, 0.20, 0.22, 0.15, 0.09, 0.19, 0.42];
[n,~] = size(x);

%% Sampling
% MCMC Parameters
nchains = 2; % How Many Chains?
nburnin = 1e3; % How Many Burn-in Samples?
nsamples = 5e4;  %How Many Recorded Samples?
nthin = 1; % How Often is a Sample Recorded?
doparallel = 0; % Parallel Option

% Assign Matlab Variables to the Observed Nodes
datastruct = struct('x',x,'n',n);

% Initial Values to Supply to WinBugs
for i=1:nchains
    S.r = 0;
    init0(i) = S;
end

if ~run_model
    load OptionalStopping samples stats
else
    if ~sampler
    % Use WinBUGS to Sample
    tic
[samples, stats] = matbugs(datastruct, ...
    fullfile(pwd, 'Correlation_1.txt'), ...
    'init', init0, ...
    'nChains', nchains, ...
    'view', 1, 'nburnin', nburnin, 'nsamples', nsamples, ...
    'thin', nthin, 'DICstatus', 0, 'refreshrate',100, ...
    'monitorParams', {'r'}, ...
    'Bugdir', 'C:/Program Files/WinBUGS14');
else
    % Use JAGS to Sample
    tic
    fprintf( 'Running JAGS ...\n' );
    [samples, stats] = matjags( ...
        datastruct, ...
        fullfile(pwd, 'Correlation_1.txt'), ...
        init0, ...
        'doparallel' , doparallel, ...
        'nchains', nchains,...
        'nburnin', nburnin,...
        'nsamples', nsamples, ...
        'thin', nthin, ...
        'monitorparams', {'r'}, ...
        'savejagsoutput' , 1 , ...
        'verbosity' , 1 , ...
        'cleanup' , 0 , ...
        'workingdir' , 'tmpjags' );
    toc
end;
    save OptionalStopping samples stats
end;

%% Analysis
crit= 0;eps = .015; binsc=[-1+eps/2:eps:1-eps/2];binse=[-1:eps:1];
figure(1);clf;hold on;
set(gcf,'units','norm','pos',[.2 .2 .75/1.2 .6/1.2],'paperpositionmode','auto','color','w');
subplot(121);hold on;
         ph = plot(x(:,1),x(:,2),'ko');
         set(ph,'markersize',6,'markerfacecolor',.5*ones(1,3),'linewidth',1,'markeredgecolor','k');
axis([0 210 0 .5]);
set(gca,'box','on','fontsize',14,'xtick',[0:50:200],'ytick',[0:.1:.5]);
xlabel('Number of Subjects','fontsize',16);
ylabel('Effect Size','fontsize',16);

subplot(122);hold on;
count = histc(reshape(samples.r,1,[]),binse);
count = count(1:end-1);
count = count/sum(count)/eps;
ph = plot(binsc,count,'k-');
set(ph,'linewidth',1);
if exist('aksdensity')
   [f,xi] = ksdensity(reshape(samples.r,1,[]),'support','unbounded','kernel','normal');
ph = plot(xi,f,'k-');
set(ph,'linewidth',1.5);
[val ind] = min(abs(xi-crit));
ph = plot(crit,f(ind),'ko');
set(ph,'markerfacecolor','w','markersize',8,'linewidth',1.5);
else
    [val ind] = min(abs(binsc-crit));
ph = plot(crit,count(ind),'ko');
set(ph,'markerfacecolor','w','markersize',8,'linewidth',1.5);
end;
ph = plot(crit,1/2,'ko');
set(ph,'markerfacecolor','k','markersize',8,'linewidth',1.5);
tmp = corrcoef(x(:,1),x(:,2));
ph = plot(ones(1,2)*tmp(1,2),[0 ceil(max(count))],'k--');
ph = plot([-1 1],[1 1]/2,'k:');
set(ph,'linewidth',1.5);
axis([-1 1 0 ceil(max(count))]);
set(gca,'box','on','fontsize',14,'xtick',[-1:.5:1],'ytick',[],'ticklength',[0 0]);
xlabel('Correlation','fontsize',16);
ylabel('Density','fontsize',16);

if exist('aksdensity')
    v1 = 1/2; v2 = f(ind);
else
    v1 = 1/2; v2 = count(ind);
end;
bf = [v1/v2 log(v1)-log(v2)];
disp(bf);

sc = .3; sc2 = .95;
arrow([crit max(v1,v2)*sc2],[crit min(v1,v2)*1/sc],10,30);
arrow([crit min(v1,v2)*1/sc],[crit max(v1,v2)*sc2],10,30);
str = sprintf(' $%d\\times$',round(max([bf(1) 1/bf(1)])));
th = text(crit,mean([v1 v2]),str);
set(th,'vert','mid','hor','left','fontsize',14,'interp','latex');
