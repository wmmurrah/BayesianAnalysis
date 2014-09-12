% Multinomial Processing Trees

clear;

sampler = 1; % Choose 0=WinBUGS, 1=JAGS

% Expected values used for nuisance parameter u
expectedu = [0.3207  0.4261 0.6582];

%% Riefer et al (2002) Data
nsubjs = 21;
ntrials = 3;
npairs = 20;
trial = zeros(nsubjs,4,ntrials);
trial(:,:,1) = [2,4,4,10;2,1,3,14;2,2,5,11;6,0,4,10;1,0,4,15;1,0,2,17;1,2,4,13;...
    4,1,6,9;5,1,4,10;1,0,9,10;5,0,3,12;0,1,6,13;1,5,7,7;1,1,4,14;2,2,3,13;...
    2,1,5,12;2,0,6,12;1,0,5,14;2,1,8,9;3,0,2,15;1,2,3,14];
trial(:,:,2) = [7,5,3,5;5,2,3,10;6,2,7,5;9,4,2,5;2,2,7,9;1,3,3,13;5,0,5,10;...
    7,3,4,6;7,3,6,4;4,1,10,5;9,1,2,8;3,1,6,10;3,5,9,3;2,0,6,12;8,0,3,9;...
    3,2,7,8;7,1,5,7;2,1,6,11;5,3,5,7;5,0,6,9;6,2,2,10];
trial(:,:,3) = [14,3,1,2;12,3,1,4;18,0,1,1;15,3,0,2;7,1,10,2;3,6,11,0;8,4,3,5;...
    17,1,1,1;13,4,3,0;11,6,1,2;16,1,2,1;10,1,3,6;7,13,0,0;8,4,3,5;16,1,1,2;...
    5,4,7,4;15,0,5,0;6,3,6,5;17,2,0,1;17,1,0,2;8,3,6,3];

%% WinBUGS
% WinBUGS Sampling Parameters
nchains=2; % number of chains
nburnin=1e3; % number of burn-in samples
nsamples=1e4; % number of samples

for triali = 1:ntrials
    % Data to Supply to WinBugs
datastruct = struct('k',sum(trial(:,:,triali),1),'n',nsubjs*npairs,'u',expectedu(triali));

% Initial Values to Supply to WinBugs
for i=1:nchains
	S.c = rand;
	S.r = rand;
	init0(i) = S;
end

% Sampling
if ~sampler
    % Use WinBUGS to Sample
    tic
[samples{triali}, stats{triali}] = matbugs(datastruct, ...
	fullfile(pwd, 'MPT_1_answer.txt'), ...
	'init', init0, ...
	'nChains', nchains, ...
	'view', 0, 'nburnin', nburnin, 'nsamples', nsamples, ...
	'thin', 1, 'DICstatus', 0, 'refreshrate',100, ...
	'monitorParams', {'c','r'}, ...
	'Bugdir', 'C:/Program Files/WinBUGS14');
toc
else
    % Use JAGS to Sample
    tic
    doparallel = 0;
    fprintf( 'Running JAGS with chains serially...\n' );
    [samples{triali}, stats{triali}] = matjags( ...
        datastruct, ...
        fullfile(pwd, 'MPT_1_answer.txt'), ...
        init0, ...
        'doparallel' , doparallel, ...
        'nchains', nchains,...
        'nburnin', nburnin,...
        'nsamples', nsamples, ...
        'thin', 1, ...
        'monitorparams',    {'c','r'},...
        'savejagsoutput' , 1 , ...
        'verbosity' , 1 , ...
        'cleanup' , 0 , ...
        'workingdir' , 'tmpjags' );
    toc
end;
end;


%% Analysis
% Draw The Posteriors
figure(1);clf;
set(gcf,'units','norm','pos',[.2 .2 .7 .45],'paperpositionmode','auto');
eps=.01;binsc=[eps/2:eps:1-eps/2];binse=[0:eps:1];
ls=char('-','--',':');
col=char('k','k','k');
% c
subplot(131);hold on;
for i=1:3
	count=histc(reshape(samples{i}.c,1,[]),binse);
count=count(1:end-1);
count=count/sum(count)/eps;
ph=plot(binsc,count,'k-');
	set(ph,'linewidth',1.5,'linestyle',deblank(ls(i,:)),'color',col(i,:));
end;
[lh oh] = legend('Trial 1','Trial 2','Trial 6','location','northwest');
set(lh,'box','off','fontsize',16);
xlim([0 1]);
set(gca,'fontsize',14,'xtick',[0:.2:1],...
	'ytick',[],'box','on','ticklength',[0 0]);
xlabel('c','fontsize',16);
ylabel('Probability Density','fontsize',16);
% r
subplot(132);hold on;
for i=1:3
	count=histc(reshape(samples{i}.r,1,[]),binse);
count=count(1:end-1);
count=count/sum(count)/eps;
ph=plot(binsc,count,'k-');
	set(ph,'linewidth',1.5,'linestyle',deblank(ls(i,:)),'color',col(i,:));
end;
xlim([0 1]);
set(gca,'fontsize',14,'xtick',[0:.2:1],...
	'ytick',[],'box','on','ticklength',[0 0]);
xlabel('r','fontsize',16);

