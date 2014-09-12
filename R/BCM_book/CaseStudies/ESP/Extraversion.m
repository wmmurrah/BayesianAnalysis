% Correlation Coefficient With Measurement Error for Performance vs Extraversion

clear;

sampler = 0; % Choose 0=WinBUGS, 1=JAGS
run_model = 1; % set 0 to load samples, or 1 to run WinBUGS


%% Data
k= [36 32 36 36 28 40 40 24 36 36 28 40 28 36 20 24 24 16 20 32 40 32 36 24 28 44 40 36 40 32 32 ...
    40 28 20 24 32 24 24 20 28 24 28 28 32 20 44 16 36 32 28 24 32 40 28 32 32 28 24 28 40 28 20 ...
    20 20 24 24 36 28 20 20 40 32 20 36 28 28 24 20 28 32 48 24 32 32 40 40 40 36 36 32 20 28 40 ...
    32 20 20 16 16 28 40];
x = [50, 80, 79, 56, 50, 80, 53, 84, 74, 67, 50, 45, 62, 65, 71, 71, 68, 63, 67, 58, 72, 73, 63, 54, ...
    63, 70, 81, 71, 66, 74, 70, 84, 66, 73, 78, 64, 54, 74, 62, 71, 70, 79, 66, 64, 62, 63, 60, 56, 72, ...
    72, 79, 67, 46, 67, 77, 55, 63, 44, 84, 65, 41, 62, 64, 51, 46, 53, 26, 67, 73, 39, 62, 59, 75, 65, ...
    60, 69, 63, 69, 55, 63, 86, 70, 67, 54, 80, 71, 71, 55, 57, 41, 56, 78, 58, 76, 54, 50, 61, 60, 32, 67];
 nsubjs = length(k);
ntrials = 60;
sigmax = 3;

%% Sampling
% MCMC Parameters
nchains = 2; % How Many Chains?
nburnin = 5e3; % How Many Burn-in Samples?
nsamples = 1e4;  %How Many Recorded Samples?
nthin = 1; % How Often is a Sample Recorded?
doparallel = 0; % Parallel Option

% Assign Matlab Variables to the Observed Nodes
datastruct = struct('k',k,'x',x,'lambdax',1/sigmax^2,'ntrials',ntrials,'nsubjs',nsubjs);

% Initial Values to Supply to WinBugs
for i=1:nchains
    S.r = 0;
    S.mu = zeros(1,2);
    S.lambda = ones(1,2);
    init0(i) = S;
end

if ~run_model
    load(['Extraversion_' int2str(sigmax)],'samples','stats');
else
    if ~sampler
    % Use WinBUGS to Sample
    tic
[samples, stats] = matbugs(datastruct, ...
    fullfile(pwd, 'Extraversion.txt'), ...
    'init', init0, ...
    'nChains', nchains, ...
    'view', 1, 'nburnin', nburnin, 'nsamples', nsamples, ...
    'thin', nthin, 'DICstatus', 0, 'refreshrate',100, ...
    'monitorParams', {'r','mu','sigma','theta'}, ...
    'Bugdir', 'C:/Program Files/WinBUGS14');
else
    % Use JAGS to Sample
    tic
    fprintf( 'Running JAGS ...\n' );
    [samples, stats] = matjags( ...
        datastruct, ...
        fullfile(pwd, 'Extraversion.txt'), ...
        init0, ...
        'doparallel' , doparallel, ...
        'nchains', nchains,...
        'nburnin', nburnin,...
        'nsamples', nsamples, ...
        'thin', nthin, ...
        'monitorparams', {'r','mu','sigma','theta'}, ...
        'savejagsoutput' , 1 , ...
        'verbosity' , 1 , ...
        'cleanup' , 0 , ...
        'workingdir' , 'tmpjags' );
    toc
end;
    save(['Extraversion_' int2str(sigmax)],'samples','stats')
end;

%% Analysis
crit = 0;
nplot = 10;eps = .03; binsc=[-1+eps/2:eps:1-eps/2];binse=[-1:eps:1];
figure(3);clf;hold on;
set(gcf,'units','norm','pos',[.2 .2 .75/1.2 .6/1.2],'paperpositionmode','auto','color','w');
subplot(121);hold on;
tmp1 = randperm(nsamples); tmp1 = tmp1(1:nplot);
tmp2 = ceil(rand(nplot,1)*nchains);
for i = 1:nsubjs
   for j = 1:nplot
    ph = plot([stats.mean.theta(i,2) samples.theta(tmp2(j),tmp1(j),i,2)],[stats.mean.theta(i,1) samples.theta(tmp2(j),tmp1(j),i,1)],'k-');
    set(ph,'markersize',5,'color',.7*ones(1,3));
end;
end;
for i = 1:nsubjs
         ph = plot(stats.mean.theta(i,2),stats.mean.theta(i,1),'ko');
         set(ph,'markersize',4,'markerfacecolor','w','linewidth',1.5);
end;
axis([0 100 0 1]);
set(gca,'box','on','fontsize',14,'xtick',[0:10:100],'ytick',[0:.2:1]);
xlabel('Extraversion','fontsize',16);
ylabel('Accuracy Session 1','fontsize',16);
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
tmp = corrcoef(x,k);
ph = plot(ones(1,2)*tmp(1,2),[0 ceil(max(count))],'k--');
ph = plot([-1 1],[1 1]/2,'k:');
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

sc = .95; sc2 = .95;
arrow([crit max(v1,v2)*sc2],[crit min(v1,v2)*1/sc],10,30);
arrow([crit min(v1,v2)*1/sc],[crit max(v1,v2)*sc2],10,30);
str = sprintf(' $%d\\times$',round(max([bf(1) 1/bf(1)])));
th = text(crit,mean([v1 v2]),str);
set(th,'vert','mid','hor','left','fontsize',14,'interp','latex');
