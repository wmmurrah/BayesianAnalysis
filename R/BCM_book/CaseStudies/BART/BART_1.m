% BART Model of Risky Decision Making

clear;

sampler = 0; % Choose 0=WinBUGS, 1=JAGS
run_model = 1; % set 0 to load samples, or 1 to run WinBUGS

%%  Data
p=.15;
ntrials=90;
tmp=importdata('GeorgeSober.txt','\t',1);
% tmp=importdata('BillSober.txt','\t',1);
cash=tmp.data(1:ntrials,7)~=0;
npumps=tmp.data(1:ntrials,6);
options=cash+npumps;
for j=1:ntrials
    if npumps(j)>0
        d(j,1:npumps(j)) = zeros(1,npumps(j));
    end;
	if cash(j)==1
        d(j,npumps(j)+1) = 1;
    end;
end;

return;

%% Sampling
% MCMC Parameters
nchains = 2; % How Many Chains?
nburnin = 2e3; % How Many Burn-in Samples?
nsamples = 5e3;  %How Many Recorded Samples?
nthin = 1; % How Often is a Sample Recorded?
doparallel = 0; % Parallel Option

% Assign Matlab Variables to the Observed Nodes
datastruct = struct('d',d,'ntrials',ntrials,'p',p,'options',options);

% Initial Values 
for i=1:nchains
    S.gplus = 1.2;
    S.beta = 0.5;
    init0(i) = S;
end

if ~run_model
    load BART_1 samples stats
else
    if ~sampler
        % Use WinBUGS to Sample
        tic
        [samples, stats] = matbugs(datastruct, ...
            fullfile(pwd, 'BART_1.txt'), ...
            'init', init0, ...
            'nChains', nchains, ...
            'view', 1, 'nburnin', nburnin, 'nsamples', nsamples, ...
            'thin', nthin, 'DICstatus', 0, 'refreshrate',100, ...
            'monitorParams',{'gplus','beta'}, ...
            'Bugdir', 'C:/Program Files/WinBUGS14');
        toc
    else
        % Use JAGS to Sample
        tic
        fprintf( 'Running JAGS ...\n' );
        [samples, stats] = matjags( ...
            datastruct, ...
            fullfile(pwd, 'BART_1.txt'), ...
            init0, ...
            'doparallel' , doparallel, ...
            'nchains', nchains,...
            'nburnin', nburnin,...
            'nsamples', nsamples, ...
            'thin', nthin, ...
            'monitorparams',{'gplus','beta'},...
            'savejagsoutput' , 1 , ...
            'verbosity' , 1 , ...
            'cleanup' , 0 , ...
            'workingdir' , 'tmpjags' );
        toc
    end;
    save BART_1 samples stats
end;

%% Analysis
figure(2);clf;hold on;
set(gcf,'units','norm','pos',[.2 .2 .6 .4],'paperpositionmode','auto');
subplot(131);cla;hold on;
count=hist(npumps,[1:max(npumps)]);
ph=bar([1:max(npumps)],count,'k');
set(gca,'xtick',[1:max(npumps)],'ytick',[1 max(count)],'box','on','fontsize',14,'xlim',[0 max(npumps)+1]);
xlabel('Number of Pumps','fontsize',16);
ylabel('Frequency','fontsize',16);
subplot(132);cla;hold on;
data=reshape(samples.gplus,1,[]);
eps=.02;
bins=[min(data)-eps:eps:max(data)+eps];
count=hist(data,bins);
count=count/sum(count)/eps;
ph=plot(bins,count,'k-');
set(gca,'box','on','fontsize',14,'ytick',[],'xlim',[bins(1) bins(end)]);
xlabel('\gamma^+','fontsize',16);
ylabel('Posterior Density','fontsize',16);
subplot(133);cla;hold on;
data=reshape(samples.beta,1,[]);
eps=.04;
bins=[min(data)-eps:eps:max(data)+eps];
count=hist(data,bins);
count=count/sum(count)/eps;
ph=plot(bins,count,'k-');
set(gca,'box','on','fontsize',14,'ytick',[],'xlim',[bins(1) bins(end)]);
xlabel('\beta','fontsize',16);
ylabel('Posterior Density','fontsize',16);

