%% ChasaSoon

clear;

sampler = 0; % Choose 0=WinBUGS, 1=JAGS

%% Data
nattempts = 950;
nfails = 949;
n = 50; % questions
y = [ones(nfails,1);0]; % Indicate which scores are censored
z = [nan*ones(nfails,1);30]'; % All scores except the last are unknown

%% Sampling
% MCMC Parameters
nchains = 2; % How Many Chains?
nburnin = 1e3; % How Many Burn-in Samples?
nsamples = 2e3;  %How Many Recorded Samples?
nthin = 1; % How Often is a Sample Recorded?
doparallel = 0; % Parallel Option

% Assign Matlab Variables to the Observed Nodes
datastruct = struct('nattempts',nattempts,'n',n,'z',z,'y',y);

% Initialize Unobserved Variables
for i=1:nchains
    S.theta = .8;
    if sampler == 1
        S.z = [20*ones(nfails,1);30];
    end;
    init0(i) = S;
end

if ~sampler
    % Use WinBUGS to Sample
    tic
[samples, stats] = matbugs(datastruct, ...
    fullfile(pwd, 'ChaSaSoon.txt'), ...
    'init', init0, ...
    'nChains', nchains, ...
    'view', 1, 'nburnin', nburnin, 'nsamples', nsamples, ...
    'thin', nthin, 'DICstatus', 0, 'refreshrate',1, ...
    'monitorParams', {'theta','z'}, ...
    'Bugdir', 'C:/Program Files/WinBUGS14');
else
    % Use JAGS to Sample
    tic
    fprintf( 'Running JAGS ...\n' );
    [samples, stats] = matjags( ...
        datastruct, ...
        fullfile(pwd, 'ChaSaSoonJ.txt'), ...
        init0, ...
        'doparallel' , doparallel, ...
        'nchains', nchains,...
        'nburnin', nburnin,...
        'nsamples', nsamples, ...
        'thin', nthin, ...
        'monitorparams', {'theta','z'}, ...
        'savejagsoutput' , 1 , ...
        'verbosity' , 1 , ...
        'cleanup' , 0 , ...
        'workingdir' , 'tmpjags' );
    toc
end;

%% Analysis
% Credible Interval
cred = 0.95;
b1 = (1-cred)/2; b2 = 1-b1;
val = sort(reshape(samples.theta,1,[]));
lo = val(round(b1*nsamples*nchains));
hi = val(round(b2*nsamples*nchains));
disp(sprintf('%d percent credible interval is [%1.3f, %1.3f]',cred*100,lo,hi));
% Posterior Over Rate
figure(1);clf;hold on;
eps=.004; binsc = eps/2:eps:1-eps/2; binse = 0:eps:1;
count = histc(reshape(samples.theta,1,[]),binse);
count = count(1:end-1);
count = count/sum(count)/eps;
ph = plot(binsc,count,'k-');
[val ind] = max(count);
th = text(binsc(ind),val*1.1,sprintf('%1.2f - %1.2f',lo,hi));
set(th,'hor','cen','fontsize',14);
th = text(binsc(ind),val*1.2,sprintf('%d%%',cred*100));
set(th,'hor','cen','fontsize',14);
set(gca,'box','on','fontsize',14);
set(gca,'ylim',get(gca,'ylim')*1.4);
xlabel('Rate','fontsize',16);
ylabel('Posterior Density','fontsize',16);

