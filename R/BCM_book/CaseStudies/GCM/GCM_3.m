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
nsamples = 1e4;  %How Many Recorded Samples?
nthin = 1; % How Often is a Sample Recorded?
doparallel = 0; % Parallel Option

% Assign Matlab Variables to the Observed Nodes
datastruct = struct('y',x','nsubj',nsubj,'nstim',nstim,'n',n,'a',a,'d1',d1,'d2',d2);

% Initial Values to Supply to WinBugs
for i=1:nchains
    S.muctmp = 0.5;
    S.muwtmp = 0.25;
    S.delta = 0.5;
    S.sigmactmp = 0.5;
    S.sigmawtmp = 0.5;
    S.zg = floor(rand(nsubj,1)*2);
    S.zc = floor(rand(nsubj,1)*2);
    S.phic = 1/2;
    S.phig = 1/2;
    init0(i) = S;
end

if ~run_model
    load GCM_3 samples stats
else
    if ~sampler
        % Use WinBUGS to Sample
        tic
        [samples, stats] = matbugs(datastruct, ...
            fullfile(pwd, 'GCM_3.txt'), ...
            'init', init0, ...
            'nChains', nchains, ...
            'view', 1, 'nburnin', nburnin, 'nsamples', nsamples, ...
            'thin', nthin, 'DICstatus', 0, 'refreshrate',100, ...
            'monitorParams', {'z','phig','phic','muc','sigmac','muw','sigmaw','c','w','predy','wpredg','cpredg','delta','predyg'}, ...
            'Bugdir', 'C:/Program Files/WinBUGS14');
        toc
    else
        % Use JAGS to Sample
        tic
        fprintf( 'Running JAGS ...\n' );
        [samples, stats] = matjags( ...
            datastruct, ...
            fullfile(pwd, 'GCM_3J.txt'), ...
            init0, ...
            'doparallel' , doparallel, ...
            'nchains', nchains,...
            'nburnin', nburnin,...
            'nsamples', nsamples, ...
            'thin', nthin, ...
            'monitorparams', {'z','phig','phic','muc','sigmac','muw','sigmaw','c','w','predy','wpredg','cpredg','delta','predyg'}, ...
            'savejagsoutput' , 1 , ...
            'verbosity' , 1 , ...
            'cleanup' , 0 , ...
            'workingdir' , 'tmpjags' );
        toc
    end;
    save GCM_3 samples stats
end;

%% Analysis

% Latent Assignment
figure(5);clf;hold on;
set(gcf,'units','norm','pos',[.2 .2 .6 .4],'paperpositionmode','auto');
count=[];
for i=1:nsubj
    count=[count;hist(reshape(samples.z(:,:,i),1,[]),[1 2 3])/nsamples/nchains];
end;
[val ind]=max(count');
[val2 ind2]=sort(ind+val,'descend');
ph=plot(count(ind2,3),'k^-');
set(ph,'linewidth',1,'markersize',10,'markerfacecolor','w');
ph=plot(count(ind2,2),'ko-');
set(ph,'linewidth',1,'markersize',8,'markerfacecolor','w');
ph=plot(count(ind2,1),'ks-');
set(ph,'linewidth',1,'markersize',5,'markerfacecolor','w');
[lh oh]=legend('Contaminant','Attend Position','Attend Height','Location','east');
axis([0 nsubj+1 -.05 1.05]);
xlabel('Subject','fontsize',18);
ylabel('Membership Probability','fontsize',18);
set(gca,'fontsize',16,'xtick',[1 10:10:30 40],'ytick',[0:.2:1],'box','on','ticklength',[ 0 0]);
set(lh,'position',[0.4155 0.4185 0.1667 0.2241]);

% Posterior Predictive for Groups Subjects
figure(6);clf;hold on;
set(gcf,'units','norm','pos',[.2 .2 .6 .4],'paperpositionmode','auto');
sc=50;
bins=[0:n];
for i=1:nsubj
    grp(i)=mode(reshape(samples.z(:,:,i),1,[]));
end;
labs=char('Contaminant','Attend Position','Attend Height');
clabs=char('k','k','k');
mix=[3 2 1];
for g=1:3
    subplot(1,3,mix(g));hold on;
    th=title(deblank(labs(mix(g),:)));
    set(th,'vert','mid','hor','cen','fontsize',16);
    axis([.5 nstim+.5 -.5 n+.5]);
    set(gca,'ticklength',[0 0],'fontsize',14,'xtick',[],'ytick',[],'box','on');
         ph=plot([1:nstim],x(find(grp==g),:)','k-');
     set(ph,'color',.6*ones(1,3),'linewidth',1);
for i=1:nstim
        count=hist(reshape(samples.predyg(:,:,g,i),1,[]),bins);
        count=count/sum(count);
        for j=1:length(count)
            if count(j)>0
                ph=plot(i,bins(j),'ks');
                set(ph,'markersize',sc*(count(j)),'markerfacecolor','w','linewidth',1.5);
                            end;
        end;
    end;
             ph=plot([1:nstim],mean(x(find(grp==g),:)),'k-');
     set(ph,'color',clabs(mix(g)),'linewidth',3);
if mix(g)==1
        set(gca,'xtick',[1:nstim],'ytick',[0 n],'yticklabel',{'B','A'});
        xlabel('Stimulus','fontsize',18);
        ylabel('Category','fontsize',18);
    end;
end;
