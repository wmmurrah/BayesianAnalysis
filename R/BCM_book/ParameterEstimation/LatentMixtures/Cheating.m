%% Cheater

clear;

sampler = 0; % Choose 0=WinBUGS, 1=JAGS

%% Data
load Cheatingdata d truez; % Alzheimer's Cheating Data

return

% Constants
k = sum(d,2); % Total correct per participant
p = length(k); % Number of participants
n = 40; % Total trials

%% Sampling
% MCMC Parameters
nchains = 2; % How Many Chains?
nburnin = 1e3; % How Many Burn-in Samples?
nsamples = 2e3;  %How Many Recorded Samples?
nthin = 1; % How Often is a Sample Recorded?
doparallel = 0; % Parallel Option

% Initial Values to Supply to WinBugs
for i=1:nchains
    S.z = round(rand(1,p));
    S.phi = 0.5;
    S.mubon = 0.5;
    S.mudiff = 0.1;
    S.lambdabon = 10;
    S.lambdache = 10;
    init0(i) = S;
end

% Assign Matlab Variables to the Observed Nodes
datastruct = struct('k',k,'p',p,'n',n,'truth',truez);

% Sampling
if ~sampler
    % Use WinBUGS to Sample
    tic
[samples, stats, structarray] = matbugs(datastruct, ...
    fullfile(pwd, 'Cheating.txt'), ...
    'init', init0, ...
    'nChains', nchains, ...
    'view',1, 'nburnin', nburnin, 'nsamples', nsamples, ...
    'thin', nthin, 'DICstatus', 0, 'refreshrate',100, ...
    'monitorParams', {'z','phi','mubon','muche','lambdabon','lambdache','pc'}, ...
    'Bugdir', 'C:/Program Files/WinBUGS14');
toc
else
    % Use JAGS to Sample
    tic
    fprintf( 'Running JAGS ...\n' );
    [samples, stats] = matjags( ...
        datastruct, ...
        fullfile(pwd, 'Cheating.txt'), ...
        init0, ...
        'doparallel' , doparallel, ...
        'nchains', nchains,...
        'nburnin', nburnin,...
        'nsamples', nsamples, ...
        'thin', nthin, ...
        'monitorparams',  {'z','phi','mubon','muche','lambdabon','lambdache','pc'},...
        'savejagsoutput' , 1 , ...
        'verbosity' , 1 , ...
        'cleanup' , 0 , ...
        'workingdir' , 'tmpjags' );
    toc
end;

%% Analysis
% Proportion Correct Analysis
figure(1);clf;
set(gcf,'units','norm','paperpositionmode','auto','pos',[.2 .2 .65 .8]);

subplot(211);hold on; % Distribution of Total Scores
bins = 0:n;
x = hist(k(find(truez==0)),bins);
y = hist(k(find(truez==1)),bins);
ph = bar(bins,[x;y]');
set(ph(1),'facecolor','k','barwidth',1);
set(ph(2),'facecolor','w','barwidth',1);
set(gca,'ylim',[0 1.05*max([x y])]);
set(gca,'fontsize',15,'xlim',[0 n+1],'box','on','xtick',[0:5:n],'ytick',[0:5:20],'ticklength',[0 0]);
[lh oh] = legend({'Bona fide','Cheater'},'location','northwest');
set(lh,'box','off');
ylabel('Number of People','fontsize',18);

subplot(212);hold on; % Accuracy as a Function of Cutoff
pc=zeros(size(bins));
for i = 1:length(bins)
    t = zeros(p,1);
    t(find(k >= bins(i))) = 1;
    pc(i) = sum(t == truez);
end;
pc = pc/p;
[val ind]=max(pc);
ind=21;
bins2= 0:p ;
sc=.01;
count = hist(reshape(samples.pc,1,[]),bins2);
oount = count/sum(count);
for i = 1:length(bins2)
    if count(i) > 0
        ph=plot([0 count(i)*sc],ones(1,2)*bins2(i)/p,'k-');
        set(ph,'linewidth',2,'color',.7*ones(1,3));
    end;
end;
ph = plot(bins,pc,'k-');
set(ph,'linewidth',2);
axis([0 n+1 .4 1]);
set(gca,'fontsize',15,'box','on','xtick',[0:5:n],'ytick',[.4:.1:1],'ticklength',[0 0]);
set(gca,'ygrid','on','yaxisloc','left');
xlabel('Number of Items Recalled Correctly','fontsize',18);
ylabel('Proportion Correct','fontsize',18);
   
% Classification Analysis
figure(2);clf;hold on;
ph = plot(k,stats.mean.z,'kx');
set(ph,'markersize',6,'linewidth',2);
axis([0 41 0 1]); axis square;
set(gca,'fontsize',13,'box','on','xtick',[0:5:n],'ytick',[0:.1:1],'ticklength',[0 0]);
xlabel('Number of Items Recalled Correctly','fontsize',16);
ylabel('Cheater Classification','fontsize',16);
ph = plot([0 30],[.2 .2],'k--');
ph = plot([30 30],[0 .2],'k--');
ph = plot([0 35],[.5 .5],'k--');
ph = plot([35 35],[0 .5],'k--');

