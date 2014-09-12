% Number Concept Using Knower Levels On WOTC Data

clear;

sampler = 0; % Choose 0=WinBUGS, 1=JAGS
run_model = 1; % set 0 to load samples, or 1 to run WinBUGS

%% Load Data
load fc_given ns nz gn gnq gq ga fnq fq fa
% ns = number of subjects (children)
% nz = number of knower-level stages
% gn = maximum give-n response
% gnq = number questions for each child on give-n
% gq = matrix of questions asked to each child on give-n
% ga = matrix of answers given by each child on give-n
% fnq = number questions for each child on fast cards
% fq = matrix of questions asked to each child on fast cards
% fa = matrix of answers given by each child on fast cards

% Maximum answer in data
fn=max(fa(:));

% Hide Give-N behavior for half subjects
subjlist=[38 17 45 40 44 53];
ga(subjlist,:)=nan;

%% Sampling
% MCMC Parameters
nchains = 2; % How Many Chains?
nburnin = 1e3; % How Many Burn-in Samples?
nsamples = 1e3;  %How Many Recorded Samples?
nthin = 1; % How Often is a Sample Recorded?
doparallel = 0; % Parallel Option

% Assign Matlab Variables to the Observed Nodes
datastruct = struct('ns',ns,'fnq',fnq,'fn',fn,'fa',fa,'fq',fq,'nz',nz,'gnq',gnq,'gn',gn,'ga',ga,'gq',gq);

% Initialize Values to Supply to WinBUGS
for i=1:nchains
    S0.gv=10;
    S0.fv=10;
    S0.z=floor(rand(1,ns)*nz)+1;
    S0.pitmp=1/gn*ones(1,gn);
    S0.fpitmp=1/fn*ones(1,fn);
    init0(i) = S0;
end

if ~run_model
    load NumberConcept_3answer samples stats
else
    if ~sampler
        % Use WinBUGS to Sample
        tic
        % Use WinBUGS to Sample
        [samples, stats, structarray] = matbugs(datastruct, ...
            fullfile(pwd, 'NumberConcept_3.txt'), ...
            'init', init0, ...
            'nChains', nchains, ...
            'view', 1, 'nburnin', nburnin, 'nsamples', nsamples, ...
            'thin', nthin, 'DICstatus', 0, 'refreshrate',1, ...
            'monitorParams', {'predgaz'}, ...
            'Bugdir', 'C:/Program Files/WinBUGS14');
        toc
    else
        % Use JAGS to Sample
        tic
        fprintf( 'Running JAGS ...\n' );
        [samples, stats] = matjags( ...
            datastruct, ...
            fullfile(pwd, 'NumberConcept_3.txt'), ...
            init0, ...
            'doparallel' , doparallel, ...
            'nchains', nchains,...
            'nburnin', nburnin,...
            'nsamples', nsamples, ...
            'thin', nthin, ...
            'monitorparams', {'predgaz'},...
            'savejagsoutput' , 1 , ...
            'verbosity' , 1 , ...
            'cleanup' , 0 , ...
            'workingdir' , 'tmpjags' );
        toc
    end;
    save NumberConcept_3answer samples stats
end;

%% Analysis


% Posterior Prediction For Six Missing Data Subjects
figure(1);clf;
set(gcf,'units','norm','pos',[.1 .1 .8 .8],'paperpositionmode','auto');
bins=[1:gn];threshold=0;
sc=9;sc2=1;
eps=.000;
% Posterior Predictive Shaded Underlay
cm=-.5;cb=.99;
for z=1:length(subjlist)
    subj=subjlist(z);
    subplot(2,3,z);hold on;
    for i=1:gn
     count=hist(squeeze(samples.predgaz(1,:,subj,i)),bins);
      count=count/max(count);
        for j=1:gn
            ph=plot(i,j,'ks');
            set(ph,'markersize',10);
            set(ph,'markeredgecolor',min(1,max(0,(cm*count(j)+cb)))*ones(1,3),'markerfacecolor',min(1,max(0,(cm*count(j)+cb)))*ones(1,3));
        end;
        tmp=gq(subj,:);
        match=find(tmp==i);
        if ~isempty(match)
            for j=1:gn
                count=length(find(ga(subj,match)==j));
                count=count/length(match);
                if count>0
                    ph=plot(i,j,'ks');
                    set(ph,'markersize',sc*sqrt(count),'linewidth',2);
                end;
            end;
        end;
    end;
    set(ph,'markersize',8,'linewidth',2);
    axis([0 gn+1 0 gn+1]);
    set(gca,'box','on','fontsize',14);
    set(gca,'xtick',setdiff(unique(gq),0),'ytick',[1:gn]);
    th=title(['Child ' deblank(int2str(subj))]);
    set(th,'fontsize',16,'vert','mid');
end;
[ax,th]=suplabel('Question','x');
set(th,'fontsize',20,'vert','mid');
[ax,th]=suplabel('Answer','y');
set(th,'fontsize',20,'vert','mid');


