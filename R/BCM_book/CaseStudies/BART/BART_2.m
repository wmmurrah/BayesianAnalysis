% BART Model of Risky Decision Making

clear;

sampler = 0; % Choose 0=WinBUGS, 1=JAGS
run_model = 1; % set 0 to load samples, or 1 to run WinBUGS

%% Load Data
p = .15;
ntrials = 90;
nconds = 3;
options = zeros(nconds,ntrials);
npumps = zeros(nconds,ntrials);
for i = 1:nconds
    switch i,
        case 1, tmp = importdata('GeorgeSober.txt','\t',1);
        case 2, tmp = importdata('GeorgeTipsy.txt','\t',1);
        case 3, tmp = importdata('GeorgeDrunk.txt','\t',1);
    end;
    cash = tmp.data(1:ntrials,7) ~= 0;
    npumps(i,:) = tmp.data(1:ntrials,6);
    options(i,:) = cash'+npumps(i,:);
    for j = 1:ntrials
        if npumps(j) > 0
            d(i,j,1:npumps(i,j)) = zeros(1,npumps(i,j));
        end;
        if cash(j) == 1
            npumps(j) + 1;
            d(i,j,npumps(i,j)+1) = 1;
        end;
    end;
end;


%% Sampling
% MCMC Parameters
nchains = 2; % How Many Chains?
nburnin = 2e3; % How Many Burn-in Samples?
nsamples = 5e3;  %How Many Recorded Samples?
nthin = 1; % How Often is a Sample Recorded?
doparallel = 0; % Parallel Option

% Assign Matlab Variables to the Observed Nodes
datastruct = struct('d',d,'ntrials',ntrials,'p',p,'options',options,'nconds',nconds);

% Initial Values
for i=1:nchains
    S.mug = 1.2;
    S.sigmag = 0.1;
    S.mub = 0.8;
    S.sigmab = 0.8;
    init0(i) = S;
end

if ~run_model
    load BART_2 samples stats
else
    if ~sampler
        % Use WinBUGS to Sample
        tic
        [samples, stats] = matbugs(datastruct, ...
            fullfile(pwd, 'BART_2.txt'), ...
            'init', init0, ...
            'nChains', nchains, ...
            'view', 1, 'nburnin', nburnin, 'nsamples', nsamples, ...
            'thin', nthin, 'DICstatus', 0, 'refreshrate',100, ...
            'monitorparams',{'gplus','beta','mug','sigmag','mub','sigmab'},...
            'Bugdir', 'C:/Program Files/WinBUGS14');
        toc
    else
        % Use JAGS to Sample
        tic
        fprintf( 'Running JAGS ...\n' );
        [samples, stats] = matjags( ...
            datastruct, ...
            fullfile(pwd, 'BART_2J.txt'), ...
            init0, ...
            'doparallel' , doparallel, ...
            'nchains', nchains,...
            'nburnin', nburnin,...
            'nsamples', nsamples, ...
            'thin', nthin, ...
            'monitorparams',{'gplus','beta','mug','sigmag','mub','sigmab'},...
            'savejagsoutput' , 1 , ...
            'verbosity' , 1 , ...
            'cleanup' , 0 , ...
            'workingdir' , 'tmpjags' );
        toc
    end;
    save BART_2 samples stats
end;


%% Analysis
figure(1);clf;hold on;
condlabs=char('Sober','Tipsy','Drunk');
set(gcf,'units','norm','pos',[.2 .2 .6 .6],'paperpositionmode','auto');
for i=1:nconds
subplot(3,3,1+(3*(i-1)));cla;hold on;
count=hist(npumps(i,:),[1:max(npumps(:))]);
ph=bar([1:max(npumps(:))],count,'k');
set(gca,'xtick',[1:max(npumps(:))],'ytick',[1 max(count)],'box','on','fontsize',14,'xlim',[0 max(npumps(:))+1]);
axis([0.5 9.5 0 34]);
if i==nconds
    xlabel('Number of Pumps','fontsize',16);
end;
th=ylabel(deblank(condlabs(i,:)),'fontsize',16,'rot',0,'vert','mid','hor','right');

subplot(3,3,2+(3*(i-1)));cla;hold on;
data=reshape(samples.gplus(:,:,i),1,[]);
eps=.02;
bins=[min(samples.gplus(:))-eps:eps:max(samples.gplus(:))+eps];
count=hist(data,bins);
count=count/sum(count)/eps;
ph=plot(bins,count,'k-');
set(gca,'box','on','fontsize',14,'ytick',[],'xlim',[bins(1) bins(end)]);
if i==nconds
    xlabel('\gamma^+','fontsize',16);
end;
subplot(3,3,3+(3*(i-1)));cla;hold on;
data=reshape(samples.beta(:,:,i),1,[]);
eps=.04;
bins=[min(samples.beta(:))-eps:eps:max(samples.beta(:))+eps];
count=hist(data,bins);
count=count/sum(count)/eps;
ph=plot(bins,count,'k-');
set(gca,'box','on','fontsize',14,'ytick',[],'xlim',[bins(1) bins(end)]);
if i==nconds
    xlabel('\beta','fontsize',16);
end;
end;

