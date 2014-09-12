% GCM on Kruschke Data

clear;

    run_model=1; % set 0 to load samples, or 1 to run WinBUGS
    sampler=1; % Choose 0=WinBUGS, 1=JAGS

%% Load Data
load KruschkeData y d1 d2 n nstim nsubj a x;

%% Sampling
% Sampling Parameters
nchains=2; % number of chains
nburnin=1e3; % number of burn-in samples
nsamples=1e3; % number of samples

% Data to Supply to WinBugs
datastruct = struct('y',y,'nstim',nstim,'t',n*nsubj,'a',a,'d1',d1,'d2',d2);

% Initial Values to Supply to WinBugs
for i=1:nchains
	S.w = 0.5;
    S.b = 0.5;
	S.c = 1;
	init0(i) = S;
end

if ~run_model
    load GCM_1_question samples stats
else
    if ~sampler
        % Use WinBUGS to Sample
        tic
[samples, stats] = matbugs(datastruct, ...
	fullfile(pwd, 'GCM_1_question.txt'), ...
	'init', init0, ...
	'nChains', nchains, ...
	'view', 1, 'nburnin', nburnin, 'nsamples', nsamples, ...
	'thin', 1, 'DICstatus', 0, 'refreshrate',100, ...
	'monitorParams', {'c','w','b','predy'}, ...
	'Bugdir', 'C:/Program Files/WinBUGS14');
 toc
    else
        % Use JAGS to Sample
        tic
        doparallel = 1;
        fprintf( 'Running JAGS ...\n' );
        [samples, stats] = matjags( ...
            datastruct, ...
            fullfile(pwd, 'GCM_1_question.txt'), ...
            init0, ...
            'doparallel' , doparallel, ...
            'nchains', nchains,...
            'nburnin', nburnin,...
            'nsamples', nsamples, ...
            'thin', 1, ...
            'monitorparams', {'c','w','b','predy'},...
            'savejagsoutput' , 1 , ...
            'verbosity' , 1 , ...
            'cleanup' , 0 , ...
            'workingdir' , 'tmpjags' );
        toc
    end;
    save GCM_1_question samples stats
end;

%% Analysis
stats.mean.b
%% Analysis
% Joint Posterior
figure(1);clf;hold on;
ph=plot(reshape(samples.c,1,[]),reshape(samples.w,1,[]),'kx');
set(ph,'markersize',4);
axis([0 5 0 1]);
xlabel('Generalization','fontsize',16);
ylabel('Attention Weight','fontsize',16);
set(gca,'fontsize',14,'xtick',[0:5],'ytick',[0:.2:1],'box','on');
% Posterior Predictive and Individual Subjects
figure(2);clf;hold on;
set(gcf,'units','norm','pos',[.2 .2 .6 .45],'paperpositionmode','auto');
axis([.5 nstim+.5 -.5 nt+.5]);
xlabel('Stimulus','fontsize',16);
ylabel('Category Decision','fontsize',16);
set(gca,'ticklength',[0 0],'fontsize',14,'xtick',[1:nstim],'ytick',[0 nt],'yticklabel',{'B','A'},'box','on');
ph=plot([1:nstim],x,'k:');
 set(ph,'color',0*ones(1,3),'linewidth',1);
sc=2.5;
bins=[0:3:nt*nsubj];
for i=1:nstim
        count=hist(reshape(samples.predy(:,:,i),1,[]),bins);
        count=count/sum(count);
        for j=1:length(count)
            if count(j)>0
                ph=plot([i-sc*count(j) i+sc*count(j)],ones(1,2)*bins(j)/nsubj,'k-');
                set(ph,'color',.5*ones(1,3),'linewidth',3);
            end;
        end;
end;
ph=plot([1:nstim],y/nsubj,'k-');
set(ph,'linewidth',3);
