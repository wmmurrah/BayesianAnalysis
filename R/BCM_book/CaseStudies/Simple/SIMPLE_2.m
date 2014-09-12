% Hierarchical SIMPLE Model

clear;

sampler = 0; % Choose 0=WinBUGS, 1=JAGS
run_model = 1; % set 0 to load samples, or 1 to run sampler

%% Load Data
load Murdock1962 y n nsets listlength pc labs

% Set Dataset to Use
m=zeros(size(y));
dsets=6;
gsets=9;
for dset=1:gsets
    switch dset,
        case 1, nwords = 10; lag = 2; offset = 15;
        case 2, nwords = 15; lag = 2; offset = 20;
        case 3, nwords = 20; lag = 2; offset = 25;
        case 4, nwords = 20; lag = 1; offset = 10;
        case 5, nwords = 30; lag = 1; offset = 15;
        case 6, nwords = 40; lag = 1; offset = 20;
        case 7, nwords = 10; lag = 1; offset = 5; %(Generalization)
        case 8, nwords = 25; lag = 1; offset = 12.5; %(Generalization)
        case 9, nwords = 50; lag = 1; offset = 25; %(Generalization)
    end; % switch
    % Temporal Offset For Free Recall
    m(dset,1:nwords) = offset+[(nwords-1)*lag:-lag:0];
    w(dset) = nwords; l(dset) = lag;
    listlength(dset) = nwords;
end;
n(dsets+1:gsets)=1200;
labs = char(labs,'10-1','25-1','50-1');
y = y'; m = m';

%% Sampling
% MCMC Parameters
nchains = 2; % How Many Chains?
nburnin = 1e3; % How Many Burn-in Samples?
nsamples = 1e3;  %How Many Recorded Samples?
nthin = 1; % How Often is a Sample Recorded?
doparallel = 0; % Parallel Option

% Assign Matlab Variables to the Observed Nodes
datastruct = struct('gsets',gsets,'y',y,'n',n,'listlength',listlength,'m',m,'dsets',dsets,'w',w);

% Initial Values to Supply to WinBugs
for i=1:nchains
    S.c = 10;
    S.s = 10;
    S.a = [0 .5];
    init0(i) = S;
end

if ~run_model
    load SIMPLE_2 samples stats
else
    % Sampling
    if ~sampler
        % Use WinBUGS to Sample
        tic
        [samples, stats] = matbugs(datastruct, ...
            fullfile(pwd, 'SIMPLE_2.txt'), ...
            'init', init0, ...
            'nChains', nchains, ...
            'view', 0, 'nburnin', nburnin, 'nsamples', nsamples, ...
            'thin', nthin, 'DICstatus', 0, 'refreshrate',nsamples, ...
            'monitorParams',    {'c','s','t','a'...
            'predpc[1:10,1]','predpc[1:15,2]','predpc[1:20,3]',...
            'predpc[1:20,4]','predpc[1:30,5]','predpc[1:40,6]',...
            'predpc[1:10,7]','predpc[1:25,8]','predpc[1:50,9]'},...
            'Bugdir', 'C:/Program Files/WinBUGS14');
        toc
    else
        % Use JAGS to Sample
        tic
        fprintf( 'Running JAGS ...\n' );
        [samples, stats] = matjags( ...
            datastruct, ...
            fullfile(pwd, 'SIMPLE_2.txt'), ...
            init0, ...
            'doparallel' , doparallel, ...
            'nchains', nchains,...
            'nburnin', nburnin,...
            'nsamples', nsamples, ...
            'thin', nthin, ...
            'monitorparams',   {'c','s','t','a'...
            'predpc[1:10,1]','predpc[1:15,2]','predpc[1:20,3]',...
            'predpc[1:20,4]','predpc[1:30,5]','predpc[1:40,6]',...
            'predpc[1:10,7]','predpc[1:25,8]','predpc[1:50,9]'},...
            'savejagsoutput' , 1 , ...
            'verbosity' , 1 , ...
            'cleanup' , 0 , ...
            'workingdir' , 'tmpjags' );
        toc
    end;
    save SIMPLE_2 samples stats
end;

%% Analysis

% Posterior Predictive
figure(3);clf;	hold on;
set(gcf,'units','norm','pos',[.1 .1 .6 .6],'paperpositionmode','auto');
% Drawing Constants
mark=char('o','o','o','o','o','o');
ccc=char('k','k','k','k','k','k');
hm=20;
for dset=1:gsets
    subplot(3,3,dset);hold on;
        for i=1:hm
        r1 = ceil(rand*nchains);
        r2 = ceil(rand*nsamples);
        data = squeeze(samples.predpc(r1,r2,:,dset));
        ph = plot(1:listlength(dset),data(1:listlength(dset)),'k-');
        set(ph,'linewidth',2,'Color',[0.7 0.7 0.7]);
        end;
        % Draw the Probability Correct Data
    if dset<=6
        ph=plot([1:listlength(dset)],pc(dset,1:listlength(dset)),'ko');
        set(ph,'linewidth',.75,'markeredgecolor','k','markerfacecolor','w','markersize',5);
        switch dset,
            case {3,4}
                set(ph,'markersize',4);
            case {5,6}
                set(ph,'markersize',3);
        end;
    end;
    axis([0 51 0 1]);
    set(gca,'xtick',[0:10:50],'ytick',[0:.2:1],'box','on','fontsize',14);
    %[lh oh]=legend(labs,'location','eastoutside');
    th=text(51,1,deblank(labs(dset,:)));
    set(th,'fontsize',14,'hor','right','vert','top');
end;
[ax,th]=suplabel('Serial Position','x');set(th,'fontsize',16);
[ax,th]=suplabel('Probability Correct','y');set(th,'fontsize',16);

% Posteriors
figure(4);clf;hold on;
set(gcf,'units','norm','pos',[.2 .2 .7 .3],'paperpositionmode','auto');
% Drawing Constants
epsc=.17;
cx=[18:epsc:24];
cxe=[18+epsc/2:epsc:24-epsc/2];
epss=.1;
sx=[9:epss:11];
sxe=[9+epss/2:epss:11-epss/2];
set(gcf,'paperpositionmode','auto','units','norm','pos',[.1 .1 .5 .3]);
% Threshold Noise Parameter
subplot(131);hold on;
count=hist(reshape(samples.s,1,[]),sx);
count=count(1:end-1);count=count/sum(count)*epss;
keep=find(count>1e-12);
ph=bar(sxe(keep),count(keep),'k');
set(gca,'ytick',[],'box','on','ticklength',[0 0],'fontsize',14);
set(gca,'xlim',sx([1 end]));
xlabel('Threshold Noise (s)','fontsize',14);
ylabel('Posterior Density','fontsize',14);
% Distinctiveness Parameter
subplot(132);hold on;
count=hist(reshape(samples.c,1,[]),cx);
count=count(1:end-1);count=count/sum(count)*epsc;
keep=find(count>1e-12);
ph=bar(cxe(keep),count(keep),'k');
set(gca,'ytick',[],'box','on','ticklength',[0 0],'fontsize',14);
set(gca,'xlim',cx([1 end]));
xlabel('Distinctiveness (c)','fontsize',14);
ylabel('Posterior Density','fontsize',14);
% Threshold Parameter as a Function Of List Length
subplot(133);hold on;
howmany=200;
keep=ceil(rand(howmany,1)*nsamples);
wdom=[1:50];
for i=1:howmany
    predt=samples.a(1,keep(i),1)*wdom+samples.a(1,keep(i),2);
    ph=plot(wdom,predt,'ko');
    set(ph,'markersize',2,'markerfacecolor',[.5 .5 .5],...
        'markeredgecolor',[.5 .5 .5]);
end;
set(gca,'box','on','ytick',[0:.2:1],'ylim',[0 1],...
    'xtick',[1 10:10:50],'xlim',[0 51],'fontsize',14);
xlabel('Item List Length (W)','fontsize',14);
ylabel('Threshold (t)','fontsize',14);