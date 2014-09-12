% Psychometric Function 2, including contaminant model

clear;

sampler = 0; % Choose 0=WinBUGS, 1=JAGS
run_model = 1; % set 0 to load samples, or 1 to run WinBUGS

%% Data
datarow = 1;
fid = fopen('data_x.txt','r');
while ~feof(fid)
    dataline=fgetl(fid);
    dataline=regexp(dataline,'\t','split');
    for c=1:28
        x(datarow,c)=str2double(dataline{c});
    end
    datarow=datarow+1;
end
fclose(fid);
datarow = 1;
fid = fopen('data_n.txt','r');
while ~feof(fid)
    dataline = fgetl(fid);
    dataline = regexp(dataline,'\t','split');
    for c=1:28
        n(datarow,c) = str2double(dataline{c});
    end
    datarow = datarow+1;
end
fclose(fid);
datarow = 1;
fid=fopen('data_r.txt','r');
while ~feof(fid)
    dataline = fgetl(fid);
    dataline = regexp(dataline,'\t','split');
    for c = 1:28
        r(datarow,c) = str2double(dataline{c});
    end
    datarow = datarow+1;
end
fclose(fid);
datarow = 1;
fid = fopen('data_rprop.txt','r');
while ~feof(fid)
    dataline = fgetl(fid);
    dataline = regexp(dataline,'\t','split');
    for c = 1:28
        rprop(datarow,c) = str2double(dataline{c});
    end
    datarow = datarow+1;
end
fclose(fid);

xmean = [318.888,311.0417,284.4444,301.5909,296.2000,305.7692,294.6429,280.3571];
nstim = [27, 24, 27, 22, 25, 26, 28, 28];
nsubjs = 8;

%% Sampling
% MCMC Parameters
nchains = 2; % How Many Chains?
nburnin = 5e3; % How Many Burn-in Samples?
nsamples = 1e4;  %How Many Recorded Samples?
nthin = 1; % How Often is a Sample Recorded?
doparallel = 0; % Parallel Option

% Assign Matlab Variables to the Observed Nodes
datastruct = struct('x',x,'n',n,'r',r,'xmean',xmean,'nstim',nstim,'nsubjs',nsubjs);

% Initial Values to Supply to WinBugs
for i=1:nchains
    S.mua = 0;
    S.mub = 0;
    S.sigmaa = 1;
    S.sigmab = 1;
    S.alpha = -2 + 4.*rand(1,nsubjs);
    S.beta = 0.5.*rand(1,nsubjs);
    init0(i) = S;
end

if ~run_model
    load Psychophysical_2 samples stats
else
    if ~sampler
        % Use WinBUGS to Sample
        tic
        [samples, stats] = matbugs(datastruct, ...
            fullfile(pwd, 'Psychophysical_2.txt'), ...
            'init', init0, ...
            'nChains', nchains, ...
            'view', 1, 'nburnin', nburnin, 'nsamples', nsamples, ...
            'thin', nthin, 'DICstatus', 0, 'refreshrate',100, ...
            'monitorparams',{'alpha','beta','z'},...
            'Bugdir', 'C:/Program Files/WinBUGS14');
        toc
    else
        % Use JAGS to Sample
        tic
        fprintf( 'Running JAGS ...\n' );
        [samples, stats] = matjags( ...
            datastruct, ...
            fullfile(pwd, 'Psychophysical_2.txt'), ...
            init0, ...
            'doparallel' , doparallel, ...
            'nchains', nchains,...
            'nburnin', nburnin,...
            'nsamples', nsamples, ...
            'thin', nthin, ...
            'monitorparams',{'z','alpha','beta'},...
            'savejagsoutput' , 1 , ...
            'verbosity' , 1 , ...
            'cleanup' , 0 , ...
            'workingdir' , 'tmpjags' );
        toc
    end;
    save Psychophysical_2 samples stats
end;

%% Analysis
% Basic Model
load Psychophysical_1 samples stats
alpha_all = squeeze([samples.alpha(1,:,:) samples.alpha(2,:,:)]);
beta_all = squeeze([samples.beta(1,:,:) samples.beta(2,:,:)]);
alpha_avg = stats.mean.alpha;
beta_avg = stats.mean.beta;
alpha_range = zeros(nsubjs,20);
beta_range= zeros(nsubjs,20);
for s=1:nsubjs
    alpha_range(s,:)=randsample(alpha_all(:,s),20);
    beta_range(s,:)=randsample(beta_all(:,s),20);
end
% Contaminant Model
load Psychophysical_2 samples stats
alpha_all2 = squeeze([samples.alpha(1,:,:) samples.alpha(2,:,:)]);
beta_all2 = squeeze([samples.beta(1,:,:) samples.beta(2,:,:)]);
alpha_avg2 = squeeze(mean(alpha_all2,1));
beta_avg2 = squeeze(mean(beta_all2,1));
alpha_range2 = zeros(nsubjs,20);
beta_range2 = zeros(nsubjs,20);
for s=1:nsubjs
    alpha_range2(s,:)=randsample(alpha_all2(:,s),20);
    beta_range2(s,:)=randsample(beta_all2(:,s),20);
end
zmean=squeeze(mean(mean(samples.z,1),2));
% Construct JNDs
for s=1:nsubjs
    JND(:,s) = psychfunc_inv(0.84,xmean(s),alpha_all(:,s),beta_all(:,s)) - psychfunc_inv(0.5,xmean(s),alpha_all(:,s),beta_all(:,s));
    JND2(:,s) = psychfunc_inv(0.84,xmean(s),alpha_all2(:,s),beta_all2(:,s)) - psychfunc_inv(0.5,xmean(s),alpha_all2(:,s),beta_all2(:,s));
end
% Plots with psychometric functions are called from psychfunc.m and psychfun_inv.m function-files
figure(2);clf;hold on;
set(gcf,'units','norm','position',[.1 .1 .85 .7],'paperpositionmode','auto');
col=zeros(max(nstim),3,nsubjs);
for s=1:nsubjs
    subplot(2,4,s);hold on;
    for stim=1:nstim(s)
        col(stim,:,s)=zmean(s,stim);
    end
    set(gca,'xlim',[190 410],'ylim',[-0.1 1.1],'fontsize',18,'box','on');
    hold on
    xscale = x(s,1):0.1:x(s,nstim(s));
    plot(xscale,psychfunc(xscale,xmean(s),alpha_avg(s),beta_avg(s)),'k--','linewidth',1.5); hold on
    plot(xscale,psychfunc(xscale,xmean(s),alpha_avg2(s),beta_avg2(s)),'k');
       for dat=1:nstim(s)
        plot(x(s,dat),rprop(s,dat),'s','markeredgecolor','k','markerfacecolor',1-col(dat,:,s),'markersize',5); hold on
    end;
  th=title(['Subject ',num2str(s)]);
    set(th,'fontsize',18,'vert','mid');
end;
ah=gca;
axes('position',[0,0,1,1],'visible','off');
th=text(0.5,0.02,'Test Interval (ms)','fontsize',20,'horizontalalignment','center');
th=text(0.075,0.5,'Proportion of Long Responses','fontsize',20,'rotation',90,'horizontalalignment','center');

print -depsc ../../../Content/CaseStudies/PsychophysicalFunctions/Figures/Psychophysical_1.eps

figure(3);clf;hold on;
set(gcf,'units','norm','position',[.1 .1 .85 .785],'paperpositionmode','auto','color','w');
for s=1:nsubjs
    subplot(2,4,s);hold on
    [f xi]=ksdensity(JND(:,s));
    %area(xi,f,'facecolor',[0.6 0.6 0.6],'linestyle','none'); 
    ph=plot(xi,f,'k-');
        set(ph,'linewidth',1.5);
 [f xi]=ksdensity(JND2(:,s));
    %area(xi,f,'facecolor',[0.8 0.8 0.8],'linestyle','none');
     ph=plot(xi,f,'k--');
     set(ph,'linewidth',1.5);
     set(gca,'xlim',[0 100],'ylim',[0 0.12],'box','on','xtick',[0:20:100],'ytick',[],'fontsize',14);
    th=title(['Subject ',num2str(s)]);
    set(th,'fontsize',18,'vert','mid');
    if s==1
        [lh oh]=legend('without','with','location','northwest');
        set(lh,'box','off','fontsize',14,'fontweight','normal',...
            'pos',get(lh,'pos')+[0 0 0 0]);
    end
end
ah=gca;
axes('position',[0,0,1,1],'visible','off');
text(0.5,0.02,'JND (ms)','fontsize',20,'horizontalalignment','center');
text(0.1,0.5,'Posterior Density','fontsize',20,'rotation',90,'horizontalalignment','center');
axes(ah);

