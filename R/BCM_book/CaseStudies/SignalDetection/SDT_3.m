% Hierarchical Signal Detection Theory, with Parameter Expansion

clear;

sampler = 0; % Choose 0=WinBUGS, 1=JAGS

%% Data
load Heit_Rotello sdt_d sdt_i; % sdt_d has the deduction data, sdt_i the induction

%% Sampling
% MCMC Parameters
nchains = 2; % How Many Chains?
nburnin = 1e3; % How Many Burn-in Samples?
nsamples = 1e4;  %How Many Recorded Samples?
nthin = 1; % How Often is a Sample Recorded?
doparallel = 0; % Parallel Option

% Do Analysis for Both Experimental Conditions
for dataset = 1:2
	switch dataset
		case 1, data = sdt_i;
		case 2, data = sdt_d;
	end;

	% Number of subjects
	[k,~] = size(data);

	% Hit, False-Alarm, Miss, and Correct-Rejection Counts for Conditions A and B
h = data(:,1);
f = data(:,2);
m = data(:,3);
c = data(:,4);
s = h+m;
n = f+c;
s = s(1); n = n(1); % Each subject gets same number of signal and noise trials

% Assign Matlab Variables to the Observed Nodes
datastruct = struct('h',h,'f',f,'s',s,'n',n,'k', k);

	% Initial Values to Supply to WinBugs
	for i=1:nchains
		S.mud = 0;
		S.muc = 0;
		S.lambdad = 1;
		S.lambdac = 1;
		init0(i) = S;
	end
    
    % Sampling
    if ~sampler
        % Use WinBUGS to Sample
        tic
        [samples, stats] = matbugs(datastruct, ...
            fullfile(pwd, 'SDT_3.txt'), ...
            'init', init0, ...
            'nChains', nchains, ...
            'view', 1, 'nburnin', nburnin, 'nsamples', nsamples, ...
            'thin', nthin, 'DICstatus', 0, 'refreshrate',100, ...
            'monitorParams', {'mud','muc','sigmad','sigmac'}, ...
            'Bugdir', 'C:/Program Files/WinBUGS14');
        
        toc
    else
        % Use JAGS to Sample
        tic
        doparallel = 0;
        fprintf( 'Running JAGS with chains serially...\n' );
        [samples, stats] = matjags( ...
            datastruct, ...
            fullfile(pwd, 'SDT_3.txt'), ...
            init0, ...
            'doparallel' , doparallel, ...
            'nchains', nchains,...
            'nburnin', nburnin,...
            'nsamples', nsamples, ...
            'thin', nthin, ...
            'monitorparams',   {'mud','muc','sigmad','sigmac'},...
            'savejagsoutput' , 1 , ...
            'verbosity' , 1 , ...
            'cleanup' , 0 , ...
            'workingdir' , 'tmpjags' );
        toc
    end;
    
    if dataset == 1
        samplesind = samples;
        statsind = stats;
    end;
end;

%% Analysis
% Draw Joint Posterior With Marginals
figure(1);clf;
clabs=char('k-','k--'); % d'/c contour line properties
collabs=char('''k''','[.8 .8 .8]'); % joint contour color
linelabs=char('-','--');
% number of joint plot
joint = 1e3;
% proportionate size of joint window
jointsize = .6;
subplot(221);hold on;h10 = gca;
bins1 = -1:.1:6;
bins2 = -3:.1:3;
axis([bins1(1) bins1(end) bins2(1) bins2(end)]);
set(h10,'yaxislocation','right','box','on','fontsize',13);
set(h10,'xtick',[],'ytick',[]);
subplot(222);hold on; h11 = gca;
ylabel('\mu_c','fontsize',14);
axis([0 1 bins2(1) bins2(end)]);
set(h11,'yaxislocation','right','ytick',[bins2(1):bins2(end)],...
	'box','on','xtick',[],'ticklength',[0 0],'fontsize',13);
subplot(223);hold on; h12 = gca;
th=xlabel('\mu_d','fontsize',14,'rot',0,'hor','left');
axis([bins1(1) bins1(end) 0 1]);
set(h12,'xtick',[bins1(1):bins1(end)],'box','on','ytick',[],...
	'ticklength',[0 0],'fontsize',14);
set(h10,'units','normalized','position',...
	[.1 1-jointsize-.1 jointsize jointsize]);
set(h11,'units','normalized','position',...
	[jointsize+.1+.05 1-jointsize-.1 1-.25-jointsize jointsize]);
set(h12,'units','normalized','position',...
	[.1 .1 jointsize 1-.25-jointsize]);
for i = 1:2
	subplot(h10);hold on;
	keep = ceil(rand(joint,1)*nsamples);
	switch i
		case 1, ph = plot(samplesind.mud(keep),samplesind.muc(keep),'ko');
		case 2, ph = plot(samples.mud(keep),samples.muc(keep),'ko');
	end;
	eval(['set(ph,''markeredgecolor'',' deblank(collabs(i,:)) ...
		',''markersize'',2,''markerfacecolor'',' deblank(collabs(i,:)) ');']);
	subplot(h11);hold on;
	switch i
		case 1, count = hist(reshape(samplesind.muc,1,[]),bins2);
		case 2, count = hist(reshape(samples.muc,1,[]),bins2);
	end;
	count = count/max(count);
	ph=plot(1-count,bins2,'k-');
	set(ph,'linewidth',2);
	eval(['set(ph,''color'',' deblank(collabs(i,:)) ');']);
	subplot(h12);hold on;
	switch i
		case 1, count = hist(reshape(samplesind.mud,1,[]),bins1);
		case 2, count = hist(reshape(samples.mud,1,[]),bins1);
	end;
	count=count/max(count);
	ph=plot(bins1,count,'k-');
	set(ph,'linewidth',2);
	eval(['set(ph,''color'',' deblank(collabs(i,:)) ');']);
end;
set(h11,'xlim',[1-1.2 1-0]);
set(h12,'ylim',[0 1.2]);