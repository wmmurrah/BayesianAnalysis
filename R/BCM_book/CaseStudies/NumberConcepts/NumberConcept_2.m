% Number Concept Using Knower Levels On WOTC Data

clear;

sampler = 0; % Choose 0=WinBUGS, 1=JAGS
run_model = 1; % set 0 to load samples, or 1 to run WinBUGS

%% Load Data
load fc_given ns nz gn fnq fq fa
% ns = number of subjects (children)
% nz = number of knower-level stages
% fnq = number questions for each child on fast cards
% fq = matrix of questions asked to each child on fast cards
% fa = matrix of answers given by each child on fast cards

% Answers in data
wlist=setdiff(unique(fa),0);
fn=length(wlist);
fqm=zeros(size(fa));
fam=zeros(size(fa));
for i=1:ns
    for j=1:max(fnq)
        [tf loc]=ismember(fa(i,j),wlist);
        fam(i,j)=loc;
            [tf loc]=ismember(fq(i,j),wlist);
        fqm(i,j)=loc;
end;
end;

%% Sampling
% MCMC Parameters
nchains = 2; % How Many Chains?
nburnin = 1e3; % How Many Burn-in Samples?
nsamples = 1e3;  %How Many Recorded Samples?
nthin = 1; % How Often is a Sample Recorded?
doparallel = 0; % Parallel Option

% Assign Matlab Variables to the Observed Nodes
datastruct = struct('ns',ns,'fnq',fnq,'fn',fn,'fa',fam,'fq',fqm,'nz',nz,'gn',gn);

% Initialize Values to Supply to WinBUGS
for i=1:nchains
    S0.v=2;
    S0.z=floor(rand(1,ns)*nz)+1;
    S0.fpitmp=1/fn*ones(1,fn);
    init0(i) = S0;
end

if ~run_model
    load NumberConcept_2 samples stats
else
    if ~sampler
        % Use WinBUGS to Sample
        tic
        % Use WinBUGS to Sample
        [samples, stats, structarray] = matbugs(datastruct, ...
            fullfile(pwd, 'NumberConcept_2.txt'), ...
            'init', init0, ...
            'nChains', nchains, ...
            'view', 1, 'nburnin', nburnin, 'nsamples', nsamples, ...
            'thin', nthin, 'DICstatus', 0, 'refreshrate',1, ...
            'monitorParams', {'predfa','predz','predfpi','v','z'}, ...
            'Bugdir', 'C:/Program Files/WinBUGS14');
        toc
    else
        % Use JAGS to Sample
        tic
        fprintf( 'Running JAGS ...\n' );
        [samples, stats] = matjags( ...
            datastruct, ...
            fullfile(pwd, 'NumberConcept_2.txt'), ...
            init0, ...
            'doparallel' , doparallel, ...
            'nchains', nchains,...
            'nburnin', nburnin,...
            'nsamples', nsamples, ...
            'thin', nthin, ...
            'monitorparams', {'predfa','predz','predfpi','v','z'},...
            'savejagsoutput' , 1 , ...
            'verbosity' , 1 , ...
            'cleanup' , 0 , ...
            'workingdir' , 'tmpjags' );
        toc
    end;
    save NumberConcept_2 samples stats
end;

%% Analysis

% Constants
labs=char('PN','One','Two','Three','Four','CP');
labs2=char('PN','1','2','3','4','CP');

% Posterior for BaseRate
figure(1);clf;hold on;
set(gcf,'units','norm','pos',[.1 .1 .8 .4],'paperpositionmode','auto');
count=hist(samples.predfpi,[1:fn]);
ph=bar(wlist,count/sum(count));
set(ph,'facecolor','k');
set(gca,'xlim',[0 max(wlist)+1],'xtick',[1 10:10:max(wlist)],'ytick',[],'box','on','fontsize',20,'ticklength',[0 0]);
xlabel('Number','fontsize',24);
ylabel('Probability','fontsize',24);

% Posterior Over Knower Levels
figure(2);clf;hold on;
set(gcf,'units','norm','pos',[.1 .1 .8 .8],'paperpositionmode','auto');
m=stats.mean.z;
[val ind]=sort(m);
ind=[1:ns];
bins=[1:nz];
for i=1:ns
    subplot(4,5,i);hold on;
    subj=ind(i);
    count=hist(samples.z(1,:,subj),bins);
    count=count/sum(count);
    ph=bar(bins,count);
    set(ph,'facecolor','k');
    axis([0 nz+1 0 1]);
    set(gca,'xtick',[1:6],'xticklabel',[],'ytick',[],...
        'ticklength',[.05 0],'box','on','fontsize',20);
    th=text(3,1,['Child ' int2str(subj)]);
    set(th,'fontsize',20,'vert','bot','hor','cen');
    if i==16
        xlabel({'Knower','Level'},'fontsize',20);
        ylabel({'Posterior','Mass'},'fontsize',20);
        set(gca,'xticklabel',labs2(:,1));
    end;
end;

% Posterior Predictive for Knower-Levels
figure(3);clf;
set(gcf,'units','norm','pos',[.1 .1 .8 .8],'paperpositionmode','auto');
bins=[1:fn];
sc=9;sc2=1;sc3=2;
mv=zeros(ns,1);
for i=1:ns
    mv(i)=mode(reshape(samples.z(:,:,i),1,[]));
end;
% Posteriore Predictive Shaded Underlay
cm=-.5;cb=.99;
for z=1:nz
    subplot(2,3,z);hold on;
    for i=1:gn
        count=hist(squeeze(samples.predz(:,:,z,i)),bins);
        count=count/max(count);
        for j=1:fn
            ph=plot(i,wlist(j),'ks');
            set(ph,'markersize',10);
            set(ph,'markeredgecolor',min(1,max(0,(cm*count(j)+cb)))*ones(1,3),'markerfacecolor',min(1,max(0,(cm*count(j)+cb)))*ones(1,3));
            %end;
        end;
    end;
    axis([0 gn+1 0 20.75]);
    set(gca,'box','on','fontsize',16);
    set(gca,'xtick',setdiff(unique(fq),0),'ytick',[1 10 20],'ticklength',[0 0]);
    th=title([deblank(labs(z,:)) '-Knower']);
    set(th,'fontsize',20,'vert','mid');
end;
% DataOverlay
dm=zeros(gn,fn,nz);
for z=1:nz
    match=find(mv==z);
    for i=1:length(match)
        for j=1:fnq(match(i))
            dm(fqm(match(i),j),fam(match(i),j),z)=dm(fqm(match(i),j),fam(match(i),j),z)+1;
        end;
    end;
end;
for z=1:nz
    subplot(2,3,z);hold on;
    for i=1:gn
        count=dm(i,:,z)/sum(dm(i,:,z));
        for j=1:fn
            if count(j)>0
                ph=plot(i,wlist(j),'ks');
                set(ph,'markersize',sc*sqrt(count(j)),'linewidth',sc3*count(j));
            end;
        end;
    end;
end;
[ax,th]=suplabel('Question','x');
set(th,'fontsize',24,'vert','mid');
[ax,th]=suplabel('Answer','y');
set(th,'fontsize',24,'vert','mid');

