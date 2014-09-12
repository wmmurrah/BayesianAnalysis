% GCM on Kruschke Data

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
datastruct = struct('y',y,'nstim',nstim,'t',n*nsubj,'a',a,'d1',d1,'d2',d2);

% Initial Values to Supply to WinBugs
for i=1:nchains
	S.w = 0.5;
	S.c = 1;
	init0(i) = S;
end

if ~run_model
    load GCM_1 samples stats
else
    if ~sampler
        % Use WinBUGS to Sample
        tic
[samples, stats] = matbugs(datastruct, ...
	fullfile(pwd, 'GCM_1.txt'), ...
	'init', init0, ...
	'nChains', nchains, ...
	'view', 0, 'nburnin', nburnin, 'nsamples', nsamples, ...
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
            fullfile(pwd, 'GCM_1.txt'), ...
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
    save GCM_1 samples stats
end;

%% Analysis
% Joint Posterior
figure(1);clf;hold on;
ph=plot(reshape(samples.c,1,[]),reshape(samples.w,1,[]),'kx');
set(ph,'markersize',4);
axis([0 5 0 1]);
xlabel('Generalization','fontsize',18);
ylabel('Attention Weight','fontsize',18);
set(gca,'fontsize',14,'xtick',[0:5],'ytick',[0:.2:1],'box','on');

% Posterior Predictive and Individual Subjects
figure(2);clf;hold on;
set(gcf,'units','norm','pos',[.2 .2 .6 .45],'paperpositionmode','auto');
axis([.5 nstim+.5 -.5 nt+.5]);
xlabel('Stimulus','fontsize',16);
ylabel('Category Decision','fontsize',16);
set(gca,'ticklength',[0 0],'fontsize',14,'xtick',[1:nstim],'ytick',[0 nt],'yticklabel',{'B','A'},'box','on');
%ph=plot([1:nstim],x,'k:');
 set(ph,'color',0*ones(1,3),'linewidth',1);
sc=3.5;
bins=[0:2:nt*nsubj];
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

