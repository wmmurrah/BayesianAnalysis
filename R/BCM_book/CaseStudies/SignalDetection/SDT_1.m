% Signal Detection Theory

clear;

sampler = 0; % Choose 0=WinBUGS, 1=JAGS

%% Data
dataset = 1; % Choose data
switch dataset
	case 1, % Demo
		data=[70 50 30 50;
			7 5 3 5;
			10 0 0 10];
	case 2; % Lehrner Et Al (1995) data
		data=[148 29 32 151;
			150 40 30 140;
			150 51 40 139];
end; % switch data set

% Constants
[k,~] = size(data);
% Hit, False-Alarm, Miss, and Correct-Rejection Counts
h = data(:,1);
f = data(:,2);
m = data(:,3);
c = data(:,4);
s = h+m;
n = f+c;

%% Sampling
% MCMC Parameters
nchains = 2; % How Many Chains?
nburnin = 1e3; % How Many Burn-in Samples?
nsamples = 1e4;  %How Many Recorded Samples?
nthin = 1; % How Often is a Sample Recorded?
doparallel = 0; % Parallel Option

% Assign Matlab Variables to the Observed Nodes
datastruct = struct('h',h,'f',f,'s',s,'n',n,'k',k);

% Initial Values to Supply to WinBugs
for i=1:nchains
	S.d = zeros(k,1);
	S.c = zeros(k,1);
	init0(i) = S;
end

% Sampling
if ~sampler
    % Use WinBUGS to Sample
    tic
[samples, stats] = matbugs(datastruct, ...
	fullfile(pwd, 'SDT_1.txt'), ...
	'init', init0, ...
	'nChains', nchains, ...
	'view', 1, 'nburnin', nburnin, 'nsamples', nsamples, ...
	'thin', n, 'DICstatus', 0, 'refreshrate',100, ...
	'monitorParams', {'d','c','thetah','thetaf'}, ...
	'Bugdir', 'C:/Program Files/WinBUGS14');
toc
else
    % Use JAGS to Sample
    tic
    doparallel = 1;
    fprintf( 'Running JAGS ...\n' );
    [samples, stats] = matjags( ...
        datastruct, ...
        fullfile(pwd, 'SDT_1.txt'), ...
        init0, ...
        'doparallel' , doparallel, ...
        'nchains', nchains,...
        'nburnin', nburnin,...
        'nsamples', nsamples, ...
        'thin', 1, ...
        'monitorparams',   {'d','c','thetah','thetaf'},...
        'savejagsoutput' , 1 , ...
        'verbosity' , 1 , ...
        'cleanup' , 0 , ...
        'workingdir' , 'tmpjags' );
    toc
end;

%% Analysis
% Draw The Posteriors
figure(1);clf;
set(gcf,'units','norm','pos',[.2 .2 .6 .6],'paperpositionmode','auto');
deps = 0.05;
ceps = 0.05;
reps = 0.01;
dbinse=[-2:deps:6];
cbinse=[-2:ceps:2];
ratebinse=[0:reps:1];
dbinsc=[-2+deps/2:deps:6-deps/2];
cbinsc=[-2+ceps/2:ceps:2-ceps/2];
ratebinsc=[reps/2:reps:1-reps/2];
fs = 16;
ls=char('-','--',':');
col=char('k','k','k');
lw=[1.5 1.5 2.5];
% Discriminability
subplot(221);hold on;
for i=1:3
	count=histc(reshape(samples.d(:,:,i),1,[]),dbinse);
    count = count(1:end-1);
	ph=plot(dbinsc,count/sum(count),'k-');
	set(ph,'linewidth',lw(i),'linestyle',deblank(ls(i,:)),'color',col(i,:));
end;
xlim([dbinse(1) dbinse(end)]);
set(gca,'fontsize',fs,'xtick',[dbinse(1):dbinse(end)],...
	'ytick',[],'box','on','ticklength',[0 0]);
xlabel('Discriminability','fontsize',fs);
ylabel('Probability Density','fontsize',fs);
% Bias
subplot(222);hold on;
for i=1:3
	count=histc(reshape(samples.c(:,:,i),1,[]),cbinse);
	    count = count(1:end-1);
ph=plot(cbinsc,count/sum(count),'k-');
	set(ph,'linewidth',lw(i),'linestyle',deblank(ls(i,:)),'color',col(i,:));
end;
xlim([cbinse(1) cbinse(end)]);
set(gca,'fontsize',fs,'xtick',[cbinse(1):cbinse(end)],...
	'ytick',[],'box','on','ticklength',[0 0]);
xlabel('Bias','fontsize',fs);
ylabel('Probability Density','fontsize',fs);
% Hit Rate
subplot(223);hold on;
for i=1:3
	count=histc(reshape(samples.thetah(:,:,i),1,[]),ratebinse);
	    count = count(1:end-1);
ph=plot(ratebinsc,count/sum(count),'k-');
	set(ph,'linewidth',lw(i),'linestyle',deblank(ls(i,:)),...
		'color',col(i,:));
end;
xlim([ratebinse(1) ratebinse(end)]);
set(gca,'fontsize',fs,'xtick',[0:.2:1],'ytick',[],...
	'box','on','ticklength',[0 0]);
xlabel('Hit Rate','fontsize',fs);
ylabel('Probability Density','fontsize',fs);
% Currently Need To Set This Legend Manually
switch dataset
	case 1, lh=legend('h=70/100, f=50/100',...
			'h=7/10,     f=5/10','h=10/10,   f=0/10','location','northwest');
	case 2, lh=legend('Controls','Group I',...
			'Group II','location','northwest');
end;
set(lh,'fontsize',fs,'box','off');
% False-Alarm Rate
subplot(224);hold on;
for i=1:3
	count=histc(reshape(samples.thetaf(:,:,i),1,[]),ratebinse);
	    count = count(1:end-1);
ph=plot(ratebinsc,count/sum(count),'k-');
	set(ph,'linewidth',lw(i),'linestyle',deblank(ls(i,:)),'color',col(i,:));
end;
xlim([ratebinse(1) ratebinse(end)]);
set(gca,'fontsize',12,'xtick',[0:.2:1],'ytick',[],...
	'box','on','ticklength',[0 0]);
xlabel('False Alarm Rate','fontsize',fs);
ylabel('Probability Density','fontsize',fs);

print -depsc ../../../Content/CaseStudies/SignalDetection/Figures/SDT_2.eps