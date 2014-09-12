% GCM with Individual Differences on Kruschke Data

clear;

sampler = 0; % Choose 0=WinBUGS, 1=JAGS
run_model = 1; % set 0 to load samples, or 1 to run WinBUGS

%% Load Data
load KruschkeData y d1 d2 n nstim nsubj a x;

%% Sampling
% MCMC Parameters
nchains = 2; % How Many Chains?
nburnin = 1e3; % How Many Burn-in Samples?
nsamples = 1e3;  %How Many Recorded Samples?
nthin = 1; % How Often is a Sample Recorded?
doparallel = 0; % Parallel Option

% Assign Matlab Variables to the Observed Nodes
datastruct = struct('y',x','nstim',nstim,'n',n,'nsubj',nsubj,'a',a,'d1',d1,'d2',d2);

% Initial Values to Supply to WinBugs
for i=1:nchains
    S.c = ones(1,nsubj);
    S.w = 0.5*ones(1,nsubj);
    init0(i) = S;
end

if ~run_model
    load GCM_2 samples stats
else
    if ~sampler
        % Use WinBUGS to Sample
        tic
        [samples, stats] = matbugs(datastruct, ...
            fullfile(pwd, 'GCM_2.txt'), ...
            'init', init0, ...
            'nChains', nchains, ...
            'view', 1, 'nburnin', nburnin, 'nsamples', nsamples, ...
            'thin', nthin, 'DICstatus', 0, 'refreshrate',100, ...
            'monitorParams', {'c','w','predy'}, ...
            'Bugdir', 'C:/Program Files/WinBUGS14');
        toc
    else
        % Use JAGS to Sample
        tic
        fprintf( 'Running JAGS ...\n' );
        [samples, stats] = matjags( ...
            datastruct, ...
            fullfile(pwd, 'GCM_2.txt'), ...
            init0, ...
            'doparallel' , doparallel, ...
            'nchains', nchains,...
            'nburnin', nburnin,...
            'nsamples', nsamples, ...
            'thin', nthin, ...
            'monitorparams', {'c','w','predy'},...
            'savejagsoutput' , 1 , ...
            'verbosity' , 1 , ...
            'cleanup' , 0 , ...
            'workingdir' , 'tmpjags' );
        toc
    end;
    save GCM_2 samples stats
end;

%% Analysis
% Joint Posterior
figure(3);clf;hold on;
set(gcf,'units','norm','pos',[.2 .2 .6 .6],'paperpositionmode','auto');
rng(542);
tmp=randperm(nsamples);
keep=tmp(1:20);
for i=1:nsubj
    for j=1:length(keep)
        ph=plot([stats.mean.c(i) samples.c(1,keep(j),i)],[stats.mean.w(i) samples.w(1,keep(j),i)],'k-');
        set(ph,'color',.7*ones(1,3));
    end;
end;
for i=1:nsubj
    ph=plot(stats.mean.c(i),stats.mean.w(i),'ko');
    set(ph,'markersize',6,'markerfacecolor','k');
    if ismember(i,[33,3,31])
        th=text(stats.mean.c(i),stats.mean.w(i),[' ' int2str(i)]);
        set(th,'vert','mid','hor','left','fontsize',14,'fontw','b');
    end;
end;
axis([0 4 0 1]);
xlabel('Generalization','fontsize',16);
ylabel('Attention Weight','fontsize',16);
set(gca,'fontsize',14,'xtick',[0:5],'ytick',[0:.2:1],'box','on');

% Individual Subjects Data Only
figure(4);clf;hold on;
set(gcf,'units','norm','pos',[.2 .2 .6 .6],'paperpositionmode','auto');
    sc=20;
    bins=[0:n];
for s=1:nsubj
    subplot(5,8,s);hold on;
    th=title(int2str(s));
    set(th,'vert','mid','hor','cen','fontsize',12);
    axis([.5 nstim+.5 -.5 n+.5]);
    set(gca,'ticklength',[0 0],'fontsize',14,'xtick',[],'ytick',[],'box','on');
        ph=plot([1:nstim],x(s,:),'k-');
    set(ph,'color','k','linewidth',1.5);
    if s==33
        set(gca,'xtick',[1:nstim],'ytick',[0 n],'yticklabel',{'B','A'});
        xlabel('Stimulus','fontsize',16);
         ylabel('Category','fontsize',16);
   end;
end;
