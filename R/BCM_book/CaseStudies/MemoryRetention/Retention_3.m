% Retention With Structured Individual Differences

clear;

sampler = 0; % Choose 0=WinBUGS, 1=JAGS

%% Data
t = [1 2 4 7 12 21 35 59 99 200];
slist = 1:4;
ns = length(slist);
k = [18    18    16    13     9     6     4     4     4 nan;
	17    13     9     6     4     4     4     4     4 nan;
	14    10     6     4     4     4     4     4     4 nan;
	nan   nan   nan   nan   nan   nan   nan    nan   nan nan];
nt = length(t);
n = 18;

%% Sampling
% MCMC Parameters
nchains = 2; % How Many Chains?
nburnin = 1e3; % How Many Burn-in Samples?
nsamples = 1e4;  %How Many Recorded Samples?
nthin = 1; % How Often is a Sample Recorded?
doparallel = 0; % Parallel Option

% Assign Matlab Variables to the Observed Nodes
datastruct = struct('k',k,'n',n,'t',t,'ns',ns,'nt',nt);

% Initial Values to Supply to WinBUGS
for i=1:nchains
	S.alphamu = 1/2;
	S.alphalambda = 1;
	S.betamu = .1;
	S.betalambda = 1;
	init0(i) = S;
end

% Sampling
if ~sampler
    % Use WinBUGS to Sample
[samples, stats] = matbugs(datastruct, ...
	fullfile(pwd, 'Retention_3.txt'), ...
	'init', init0, ...
	'nChains', nchains, ...
	'view', 1, 'nburnin', nburnin, 'nsamples', nsamples, ...
	'thin', nthin, 'DICstatus', 0, 'refreshrate',100, ...
	'monitorParams', {'alpha','beta','predk'}, ...
	'Bugdir', 'C:/Program Files/WinBUGS14');
toc
else
    % Use JAGS to Sample
    tic
    fprintf( 'Running JAGS ...\n' );
    [samples, stats] = matjags( ...
        datastruct, ...
        fullfile(pwd, 'Retention_3J.txt'), ...
        init0, ...
        'doparallel' , doparallel, ...
        'nchains', nchains,...
        'nburnin', nburnin,...
        'nsamples', nsamples, ...
        'thin', nthin, ...
        'monitorparams',   {'alpha','beta','predk'},...
        'savejagsoutput' , 1 , ...
        'verbosity' , 1 , ...
        'cleanup' , 0 , ...
        'workingdir' , 'tmpjags' );
    toc
end;

%% Analysis

% Draw Joint Posterior and Marginal
figure(5);clf;hold on;
jointsize=.6;
joint=1e2;
marklabs=char('''s''','''o''','''^''','''+''');
%collabs=char('''k''','''r''','''g''','''b''');
collabs=char('''k''','''k''','''k''','''k''');
linelabs=char('''-''','''--''','''-.''',''':''');
subplot(221);hold on;h10=gca;
eps=.02;
bins1=[0:eps:1];
bins2=[0:eps:1];
bins1c=[eps/2:eps:1-eps/2];
bins2c=[eps/2:eps:1-eps/2];
axis([bins1(1) bins1(end) bins2(1) bins2(end)]);
set(h10,'yaxislocation','right','box','on','fontsize',13);
set(h10,'xtick',[],'ytick',[]);
subplot(222);hold on; h11=gca;
ylabel('Baseline','fontsize',14);
axis([0 1 bins2(1) bins2(end)]);
set(h11,'yaxislocation','right','ytick',[bins2(1):bins2(end)],...
	'box','on','xtick',[],'ticklength',[0 0],'fontsize',13);
subplot(223);hold on; h12=gca;
th=xlabel('Decay Rate','fontsize',14,'rot',0,'hor','left');
axis([bins1(1) bins1(end) 0 1]);
set(h12,'xtick',[bins1(1):bins1(end)],'box','on','ytick',[],...
	'ticklength',[0 0],'fontsize',14);
set(h10,'units','normalized','position',...
	[.1 1-jointsize-.1 jointsize jointsize]);
set(h11,'units','normalized','position',...
	[jointsize+.1+.05 1-jointsize-.1 1-.25-jointsize jointsize]);
set(h12,'units','normalized','position',...
	[.1 .1 jointsize 1-.25-jointsize]);
for i=1:4
    subplot(h10);hold on;
    keep=ceil(rand(joint,1)*nsamples);
    ph=plot(samples.alpha(1,keep,i),samples.beta(1,keep,i),'ko');
    eval(['set(ph,''markeredgecolor'',' deblank(collabs(i,:)) ...
		',''markersize'',6,''markerfacecolor'',''w'');']);
      eval(['set(ph,''marker'',' deblank(marklabs(i,:)) ');']);
 subplot(h11);hold on;
    count=histc(reshape(samples.beta(:,:,i),1,[]),bins2);
    count=count(1:end-1);
    count=count/sum(count);
    ph=plot(1-count,bins2c,'k-');
      set(ph,'linewidth',2);
  eval(['set(ph,''color'',' deblank(collabs(i,:)) ');']);
       eval(['set(ph,''linestyle'',' deblank(linelabs(i,:)) ');']);
  subplot(h12);hold on;
    count=histc(reshape(samples.alpha(:,:,i),1,[]),bins1);
        count=count(1:end-1);
count=count/sum(count);
    ph=plot(bins1c,count,'k-');
    set(ph,'linewidth',2);
    eval(['set(ph,''color'',' deblank(collabs(i,:)) ');']);
       eval(['set(ph,''linestyle'',' deblank(linelabs(i,:)) ');']);
end;
set(h11,'xlim',[1-.25 1-0]);
set(h12,'ylim',[0 .35]);

% Draw Posterior Predictive Analysis
figure(6);clf;
sc=20; % Scaling Constant for Drawing Boxes
for i=1:ns
	% Subplots for Subjects
	subplot(2,2,i);hold on;
	% Plot Subject Data
	ph=plot([1:nt],k(i,:),'k-');
	set(ph,'linewidth',2);
	% Plot Posterior Predictive
	for j=1:nt
		count=hist(samples.predk(1,:,i,j),[0:n]);
		count=count/sum(count);
		for x=0:n
			if count(x+1)>0
				ph=plot(j,x,'ks');
				set(ph,'markersize',sc*sqrt(count(x+1)));
				if k(i,j)==x
					set(ph,'markerfacecolor','k');
				end;
			end;
		end;
	end;
	% Set the Axes
	axis([0 nt+1 -1 19]);
	% Title the Subplot
	th=title(['Subject ' int2str(i)]);
	set(th,'fontsize',12,'verticalalignment','mid');
	xlabel('Time Lags','fontsize',12);
	ylabel('Retention Count','fontsize',12);
	% Tidy Up the Subplot
	set(gca,'box','on','xtick',[1:nt],'xticklabel',t,'ytick',[0 18]);
end;

