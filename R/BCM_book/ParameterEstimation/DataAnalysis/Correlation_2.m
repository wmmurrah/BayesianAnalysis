% Correlation Coefficient With Measurement Error

clear;

sampler = 0; % Choose 0=WinBUGS, 1=JAGS

%% Data
% Choose dataset
dataset = 3;

x = [.8,102; 1,98; .5,100; 0.9,105; .7,103;...
    0.4,110; 1.2,99; 1.4,87; 0.6,113; 1.1,89; 1.3,93];
switch dataset
    case 3,
        sigmaerror = [.03 1];
   case 4,
        sigmaerror = [.03 10];
end;

% Constants
[n,~] = size(x);
lambdaerror = 1./sigmaerror.^2; % Precision of Measurements


%% Sampling
% MCMC Parameters
nchains = 2; % How Many Chains?
nburnin = 5e3; % How Many Burn-in Samples?
nsamples = 5e3;  %How Many Recorded Samples?
nthin = 1; % How Often is a Sample Recorded?
doparallel = 0; % Parallel Option

% Assign Matlab Variables to the Observed Nodes
datastruct = struct('x',x,'n',n,'lambdaerror',1./sigmaerror.^2);

% Initial Values to Supply to WinBugs
for i=1:nchains
    S.r = 0;
    S.mu = zeros(1,2);
    S.lambda = ones(1,2);
    init0(i) = S;
end

if ~sampler
    % Use WinBUGS to Sample
    tic
[samples, stats] = matbugs(datastruct, ...
    fullfile(pwd, 'Correlation_2.txt'), ...
    'init', init0, ...
    'nChains', nchains, ...
    'view', 1, 'nburnin', nburnin, 'nsamples', nsamples, ...
    'thin', nthin, 'DICstatus', 0, 'refreshrate',100, ...
    'monitorParams', {'r','mu','sigma'}, ...
    'Bugdir', 'C:/Program Files/WinBUGS14');
else
    % Use JAGS to Sample
    tic
    fprintf( 'Running JAGS ...\n' );
    [samples, stats] = matjags( ...
        datastruct, ...
        fullfile(pwd, 'Correlation_2.txt'), ...
        init0, ...
        'doparallel' , doparallel, ...
        'nchains', nchains,...
        'nburnin', nburnin,...
        'nsamples', nsamples, ...
        'thin', nthin, ...
        'monitorparams', {'r','mu','sigma'}, ...
        'savejagsoutput' , 1 , ...
        'verbosity' , 1 , ...
        'cleanup' , 0 , ...
        'workingdir' , 'tmpjags' );
    toc
end;

%% Analysis
figure(20+dataset);clf;hold on;
set(gcf,'units','norm','pos',[.2 .2 .6 .4],'paperpositionmode','auto');
eps = .02; binsc = -1+eps/2:eps:1-eps/2; binse = -1:eps:1;
subplot(121);hold on;
ph = plot(x(:,1),x(:,2),'ko');
set(ph,'markersize',4,'markerfacecolor','k');
for i = 1:n
    ph = plot([x(i,1)-sigmaerror(1) x(i,1)+sigmaerror(1)],ones(1,2)*x(i,2),'k-');
    ph = plot(ones(1,2)*x(i,1),[x(i,2)-sigmaerror(2) x(i,2)+sigmaerror(2)],'k-');
end;
axis([0 1.5 85 115]);
set(gca,'box','on','fontsize',14,'xtick',0:.25:1.5,'ytick',85:5:115);
xlabel('Response Time (sec)','fontsize',16);
ylabel('IQ','fontsize',16);
subplot(122);hold on;
count = histc(reshape(samples.r,1,[]),binse);
count = count(1:end-1);
count = count/sum(count)/eps;
ph = plot(binsc,count,'k-');
set(gca,'box','on','fontsize',14,'xtick',[-1:.5:1],'ytick',[]);
tmp = corrcoef(x);
ph = plot(ones(1,2)*tmp(1,2),[0 max(get(gca,'ylim'))],'k--');
xlabel('Correlation','fontsize',16);
ylabel('Posterior Density','fontsize',16);



