% Stop
clear;

sampler = 0; % Choose 0=WinBUGS, 1=JAGS
run_model = 1; % set 0 to load samples, or 1 to run WinBUGS

%% Load Data
load StopSearchData y m p v x;
% y is 20x50 binary decision data for subjects by problems
% m is 83x9 stimulus by cue matrix for stimuli used in test questions
% p is a 50x2 problem matrix, giving the stimulus numbers for test questions
% v is 9x1 vector of cue validities
% x is 9x1 vector of validities on log-odds evidence scale

% Constants
[n nc]=size(m); % number of stimuli and cues
[nq junk]=size(p); % number of questions
[ns junk]=size(y); % number of subjects


%% Sampling
% MCMC Parameters
nchains = 4; % How Many Chains?
nburnin = 1e3; % How Many Burn-in Samples?
nsamples = 1e3;  %How Many Recorded Samples?
nthin = 1; % How Often is a Sample Recorded?
doparallel = 0; % Parallel Option

% Assign Matlab Variables to the Observed Nodes
datastruct = struct('ns',ns,'nc',nc,'nq',nq,'m',m,'p',p,'y',y,'v',v,'x',x);

% Initial Values to Supply to WinBugs
for i=1:nchains
    S.gamma = 0.75; 
	init0(i) = S;
end

if ~run_model
    load Stop samples stats
else
% Sampling
if ~sampler
    % Use WinBUGS to Sample
    tic
[samples, stats] = matbugs(datastruct, ...
    fullfile(pwd, 'Stop.txt'), ...
    'init', init0, ...
    'nChains', nchains, ...
    'view', 1, 'nburnin', nburnin, 'nsamples', nsamples, ...
    'thin', nthin, 'DICstatus', 0, 'refreshrate',1, ...
    'monitorParams', {'gamma','z','phi','ypred'}, ...
    'Bugdir', 'C:/Program Files/WinBUGS14');
toc
else
    % Use JAGS to Sample
    tic
    fprintf( 'Running JAGS ...\n' );
    [samples, stats] = matjags( ...
        datastruct, ...
        fullfile(pwd, 'StopJ.txt'), ...
        init0, ...
        'doparallel' , doparallel, ...
        'nchains', nchains,...
        'nburnin', nburnin,...
        'nsamples', nsamples, ...
        'thin', nthin, ...
        'monitorparams',{'gamma','z','phi','ypred'},...
        'savejagsoutput' , 1 , ...
        'verbosity' , 1 , ...
        'cleanup' , 0 , ...
        'workingdir' , 'tmpjags' );
    toc
end;
save Stop samples stats
end;

%% Analysis
ind1 = 1:nq;
[val2 ind2] = sort(sum(y,2),'descend');
ind2 = 1:ns;
% Plot
figure(2);clf;hold on;
sc=7;sc2=7;
set(gcf,'units','norm','pos',[.2 .2 .7 .6],'paperpositionmode','auto','color','w');
axis([0 nq+1 0 ns+1]);
for i=1:nq
    for j=1:ns
            if stats.mean.ypred(ind2(j),ind1(i))>0
                ph=plot(i,j,'k+');
                set(ph,'markersize',sc*stats.mean.ypred(ind2(j),ind1(i)));
            end;
       
    end;
end;
for i=1:nq
    for j=1:ns
        if y(ind2(j),ind1(i))==1
            ph=plot(i,j,'kx');
            set(ph,'markersize',sc2,'linewidth',1);
        end;
    end
end
set(gca,'fontsize',16,'box','on','ticklength',[.0 .0],'xtick',[1 nq-6+1 nq],'ytick',[1 5:5:ns]);
xlabel('Question','fontsize',18);
ylabel('Subject','fontsize',18);

ph = plot(ones(1,2)*nq-6+.5,[0 21],'k--');

correspondence = mean(mean(y.*stats.mean.ypred+(1-y).*(1-stats.mean.ypred)))

figure(3);clf;hold on;
set(gcf,'units','norm','pos',[.2 .05 .7 .3],'paperpositionmode','auto');
axis([0 ns+1 0 1]);
ph = bar(1:ns,stats.mean.z(ind2));
set(ph,'facecolor','k');
set(gca,'fontsize',16,'box','on','ticklength',[.0 .0],'xtick',[1:ns],'ytick',[0 1],'yticklabel',{'TTB','WADD'});
xlabel('Subject','fontsize',18);
ylabel('Group','fontsize',18);

