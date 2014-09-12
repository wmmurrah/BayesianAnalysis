%% Inferring Number of Surveys Distributed and Return Rate Simultataneously

clear;

sampler = 0; % Choose 0=WinBUGS, 1=JAGS

% JAGS 3.2.0+ seems to work, but 3.1.0 did not

%% Data
dataset = 1;

nmax = 500;
switch dataset
    case 1, k = [16 18 22 25 27];
    case 2, k = [16 18 22 25 28];
end;
m = length(k);

%% Sampling
% MCMC Parameters
nchains = 2; % How Many Chains?
nburnin = 1e3; % How Many Burn-in Samples?
nsamples = 5e3;  %How Many Recorded Samples?
nthin = 1; % How Often is a Sample Recorded?
doparallel = 0; % Parallel Option

% Assign Matlab Variables to the Observed Variables
datastruct = struct('k',k,'m',m,'nmax',nmax);

% Initialize Unobserved Variables
for i=1:nchains
    S.theta = 0.5;
    S.n = round(nmax/2);
    init0(i) = S;
end

if ~sampler
    % Use WinBUGS to Sample
    tic
    [samples, stats] = matbugs(datastruct, ...
        fullfile(pwd, 'Survey.txt'), ...
        'init', init0, ...
        'nChains', nchains, ...
        'view', 1, 'nburnin', nburnin, 'nsamples', nsamples, ...
        'thin', nthin, 'DICstatus', 0, 'refreshrate',100, ...
        'monitorParams', {'theta','n'}, ...
        'Bugdir', 'C:/Program Files/WinBUGS14');
    toc
else
    % Use JAGS to Sample
    tic
    fprintf( 'Running JAGS with chains serially...\n' );
    [samples, stats] = matjags( ...
        datastruct, ...
        fullfile(pwd, 'Survey.txt'), ...
        init0, ...
        'doparallel' , doparallel, ...
        'nchains', nchains,...
        'nburnin', nburnin,...
        'nsamples', nsamples, ...
        'thin', nthin, ...
        'monitorparams', {'theta','n'}, ...
        'savejagsoutput' , 1 , ...
        'verbosity' , 1 , ...
        'cleanup' , 0 , ...
        'workingdir' , 'tmpjags' );
    toc
end;

%% Analysis
n=reshape(samples.n,1,[]);
theta=reshape(samples.theta,1,[]);
figure(50+dataset);clf;hold on;
set(gcf,'units','norm','pos',[.2 .2 .45 .55],'paperpositionmode','auto');
jointsize=.6;
subplot(221);hold on;h10=gca;
bins1=[0:15:nmax];
bins2=[0:.03:1];
axis([bins1(1) bins1(end) bins2(1) bins2(end)]);
set(h10,'yaxislocation','right','box','on','fontsize',13);
set(h10,'xtick',[],'ytick',[]);
subplot(222);hold on; h11=gca;
ylabel('Rate of Return','fontsize',16);
axis([0 1 bins2(1) bins2(end)]);
set(h11,'yaxislocation','right','ytick',[bins2(1):.2:bins2(end)],...
    'box','on','xtick',[],'ticklength',[0 0],'fontsize',13);
subplot(223);hold on; h12=gca;
th=xlabel('Number of Surveys','fontsize',16);
set(th,'rot',0,'hor','left');
axis([bins1(1) bins1(end) 0 1]);
set(h12,'xtick',[bins1(1):100:bins1(end)],'box','on',...
    'ytick',[],'ticklength',[0 0],'fontsize',14);
set(h10,'units','normalized','position',...
    [.1 1-jointsize-.1 jointsize jointsize]);
set(h11,'units','normalized','position',...
    [jointsize+.1+.05 1-jointsize-.1 1-.25-jointsize jointsize]);
set(h12,'units','normalized','position',...
    [.1 .1 jointsize 1-.25-jointsize]);
subplot(h10);hold on;
ph=plot(samples.n,samples.theta,'k.');
subplot(h11);hold on;
count=hist(theta,bins2);
count=count/max(count);
ph2=barh(bins2,1-count);
set(ph2,'facecolor','k','basevalue',1);
subplot(h12);hold on;
count=hist(n,bins1);
count=count/max(count);
ph=bar(bins1,count);
set(ph,'facecolor','k');
set(h11,'xlim',[1-1.2 1-0]);
set(h12,'ylim',[0 1.2]);ph=get(gcf,'children');axes(ph(3));hold on;
set(gca,'fontsize',14);

% Expectation
ph=plot(mean(n),mean(theta),'rx');set(ph,'markersize',12,'linewidth',2);

% Maximum Likelihood
cc=-inf;
ind=0;
for i=1:nsamples
    logL=0;
    for j=1:m
        logL=logL+gammaln(n(i)+1)-gammaln(k(j)+1)-gammaln(n(i)-k(j)+1);
        logL=logL+k(j)*log(theta(i))+(n(i)-k(j))*log(1-theta(i));
    end;
    if logL>cc
        ind=i;
        cc=logL;
    end;
end;
ph=plot(n(ind),theta(ind),'go');set(ph,'markerfacecolor','w');
