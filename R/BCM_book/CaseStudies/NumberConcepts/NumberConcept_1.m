% Number Concept Using Knower Levels On Give-N Data

clear;

sampler = 0; % Choose 0=WinBUGS, 1=JAGS
run_model = 1; % set 0 to load samples, or 1 to run WinBUGS

%% Load Data
load fc_given ns nz gn gnq gq ga
% ns = number of subjects (children)
% nz = number of knower-level stages
% gn = maximum give-n response
% gnq = number questions for each child on give-n
% gq = matrix of questions asked to each child on give-n
% ga = matrix of answers given by each child on give-n

%% Sampling
% MCMC Parameters
nchains = 2; % How Many Chains?
nburnin = 1e3; % How Many Burn-in Samples?
nsamples = 1e3;  %How Many Recorded Samples?
nthin = 1; % How Often is a Sample Recorded?
doparallel = 0; % Parallel Option

% Assign Matlab Variables to the Observed Nodes
datastruct = struct('ns',ns,'gnq',gnq,'gn',gn,'ga',ga,'gq',gq,'nz',nz);

% Initialize Values to Supply to WinBUGS
for i=1:nchains
    S0.v=2;
    S0.z=floor(rand(1,ns)*nz)+1;
    S0.pitmp=1/gn*ones(1,gn);
    init0(i) = S0;
end

if ~run_model
    load NumberConcept_1 samples stats
else
    if ~sampler
        % Use WinBUGS to Sample
        tic
        % Use WinBUGS to Sample
        [samples, stats, structarray] = matbugs(datastruct, ...
            fullfile(pwd, 'NumberConcept_1.txt'), ...
            'init', init0, ...
            'nChains', nchains, ...
            'view', 1, 'nburnin', nburnin, 'nsamples', nsamples, ...
            'thin', nthin, 'DICstatus', 0, 'refreshrate',1, ...
            'monitorParams', {'predga','predz','predpi','v','z'}, ...
            'Bugdir', 'C:/Program Files/WinBUGS14');
        toc
    else
        % Use JAGS to Sample
        tic
        fprintf( 'Running JAGS ...\n' );
        [samples, stats] = matjags( ...
            datastruct, ...
            fullfile(pwd, 'NumberConcept_1.txt'), ...
            init0, ...
            'doparallel' , doparallel, ...
            'nchains', nchains,...
            'nburnin', nburnin,...
            'nsamples', nsamples, ...
            'thin', nthin, ...
            'monitorparams', {'predga','predz','predpi','v','z'},...
            'savejagsoutput' , 1 , ...
            'verbosity' , 1 , ...
            'cleanup' , 0 , ...
            'workingdir' , 'tmpjags' );
        toc
    end;
    save NumberConcept_1 samples stats
end;

%% Analysis

% Constants
labs=char('PN','One','Two','Three','Four','CP');
labs2=char('PN','1','2','3','4','CP');

% Posterior for BaseRate
figure(1);clf;hold on;
set(gcf,'units','norm','pos',[.1 .1 .8 .5],'paperpositionmode','auto');
bins=[1:gn];
count=hist(samples.predpi,bins);
ph=bar(bins,count/sum(count));
set(ph,'facecolor','k');
set(gca,'xlim',[0 gn+1],'xtick',[1:gn],'ytick',[],'box','on','fontsize',20);
xlabel('Number','fontsize',24);
ylabel('Probability','fontsize',24);

% Posterior Over Knower Levels
figure(2);clf;hold on;
set(gcf,'units','norm','pos',[.1 .1 .8 .8],'paperpositionmode','auto');
m=stats.mean.z;
[val ind]=sort(m);
ind=[1:ns];
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
        xlabel({'Knower','Level'},'fontsize',24);
        ylabel({'Posterior','Mass'},'fontsize',24);
        set(gca,'xticklabel',labs2(:,1));
    end;
end;

% Posterior Predictive for Knower-Levels
figure(3);clf;
set(gcf,'units','norm','pos',[.1 .05 .8 .9],'paperpositionmode','auto');
bins=[1:gn];
sc=9;sc2=1;sc3=2;
eps=.000;
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
        for j=1:gn
            ph=plot(i,j,'ks');
            set(ph,'markersize',10);
            set(ph,'markeredgecolor',min(1,max(0,(cm*count(j)+cb)))*ones(1,3),'markerfacecolor',min(1,max(0,(cm*count(j)+cb)))*ones(1,3));
            %end;
        end;
    end;
    axis([0 gn+1 0 gn+1]);
    set(gca,'box','on','fontsize',16);
    set(gca,'xtick',setdiff(unique(gq),0),'ytick',[1:gn]);
    th=title([deblank(labs(z,:)) '-Knower']);
    set(th,'fontsize',20,'vert','mid');
end;
% DataOverlay
dm=zeros(gn,gn,nz);
for z=1:nz
    match=find(mv==z);
    for i=1:length(match)
        for j=1:gnq(match(i))
            dm(gq(match(i),j),ga(match(i),j),z)=dm(gq(match(i),j),ga(match(i),j),z)+1;
        end;
    end;
end;
for z=1:nz
    subplot(2,3,z);hold on;
    for i=1:gn
        count=dm(i,:,z)/sum(dm(i,:,z));
        for j=1:gn
            if count(j)>0
                ph=plot(i,j,'ks');
                set(ph,'markersize',sc*sqrt(count(j)),'linewidth',sc3*count(j));
            end;
        end;
    end;
end;
[ax,th]=suplabel('Question','x');
set(th,'fontsize',24,'vert','mid');
[ax,th]=suplabel('Answer','y');
set(th,'fontsize',24,'vert','mid');

% Posterior Prediction For Six Subjects
subjlist=[15 2 4 3 10 20];
figure(4);clf;
set(gcf,'units','norm','pos',[.1 .1 .8 .8],'paperpositionmode','auto');
bins=[1:gn];threshold=0;
sc=9;sc2=1;sc3=2;
eps=.000;
% Posterior Predictive Shaded Underlay
cm=-.5;cb=.99;
for z=1:length(subjlist)
    subj=subjlist(z);
    subplot(2,3,z);hold on;
    for i=1:gn
        count=hist(squeeze(samples.predga(1,:,subj,i)),bins);
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
                     set(ph,'markersize',sc*sqrt(count),'linewidth',sc3*count);
                end;
            end;
        end;
    end;
	set(ph,'markersize',8,'linewidth',2);
    axis([0 gn+1 0 gn+1]);
    set(gca,'box','on','fontsize',20);
    set(gca,'xtick',setdiff(unique(gq),0),'ytick',[1:gn]);
    th=title(['Child ' deblank(int2str(subj))]);
    set(th,'fontsize',24,'vert','mid');
end;
[ax,th]=suplabel('Question','x');
set(th,'fontsize',24,'vert','mid');
[ax,th]=suplabel('Answer','y');
set(th,'fontsize',24,'vert','mid');

