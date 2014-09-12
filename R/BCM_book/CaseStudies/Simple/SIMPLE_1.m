% SIMPLE Model

clear;

sampler = 0; % Choose 0=WinBUGS, 1=JAGS
run_model = 1; % set 0 to load samples, or 1 to run sampler

%% Load Data
load Murdock1962 y n nsets listlength pc labs

% Set Dataset to Use
dsets=6;
m=zeros(size(y));
for dset=1:dsets
    switch dset,
        case 1, nwords = 10; lag = 2; offset = 15;
        case 2, nwords = 15; lag = 2; offset = 20;
        case 3, nwords = 20; lag = 2; offset = 25;
        case 4, nwords = 20; lag = 1; offset = 10;
        case 5, nwords = 30; lag = 1; offset = 15;
        case 6, nwords = 40; lag = 1; offset = 20;
    end; % switch
    % Temporal Offset For Free Recall
    m(dset,1:nwords) = offset+[(nwords-1)*lag:-lag:0];
end;
y= y'; m = m';

%% Sampling
% MCMC Parameters
nchains = 2; % How Many Chains?
nburnin = 1e3; % How Many Burn-in Samples?
nsamples = 1e3;  %How Many Recorded Samples?
nthin = 1; % How Often is a Sample Recorded?
doparallel = 0; % Parallel Option

% Assign Matlab Variables to the Observed Nodes
datastruct = struct('y',y,'n',n,'listlength',listlength,'m',m,'dsets',dsets);

% Initial Values to Supply to WinBugs
for i=1:nchains
    S.c = 10*ones(dset,1);
    S.t = 0.5*ones(dset,1);
    S.s = 10*ones(dset,1);
    init0(i) = S;
end

if ~run_model
    load SIMPLE_1 samples stats
else
    % Sampling
    if ~sampler
        % Use WinBUGS to Sample
        tic
        [samples, stats] = matbugs(datastruct, ...
            fullfile(pwd, 'SIMPLE_1.txt'), ...
            'init', init0, ...
            'nChains', nchains, ...
            'view', 0, 'nburnin', nburnin, 'nsamples', nsamples, ...
            'thin', nthin, 'DICstatus', 0, 'refreshrate',nsamples, ...
            'monitorParams', {'c','s','t',...
            'predpc[1:10,1]','predpc[1:15,2]','predpc[1:20,3]',...
            'predpc[1:20,4]','predpc[1:30,5]','predpc[1:40,6]'}, ...
            'Bugdir', 'C:/Program Files/WinBUGS14');
        toc
    else
        % Use JAGS to Sample
        tic
        fprintf( 'Running JAGS ...\n' );
        [samples, stats] = matjags( ...
            datastruct, ...
            fullfile(pwd, 'SIMPLE_1.txt'), ...
            init0, ...
            'doparallel' , doparallel, ...
            'nchains', nchains,...
            'nburnin', nburnin,...
            'nsamples', nsamples, ...
            'thin', nthin, ...
            'monitorparams',   {'c','s','t',...
            'predpc[1:10,1]','predpc[1:15,2]','predpc[1:20,3]',...
            'predpc[1:20,4]','predpc[1:30,5]','predpc[1:40,6]'},...
            'savejagsoutput' , 1 , ...
            'verbosity' , 1 , ...
            'cleanup' , 0 , ...
            'workingdir' , 'tmpjags' );
        toc
    end;
    save SIMPLE_1 samples stats
end;

%% Analysis
% Posterior Predictive
figure(1);clf;hold on;
set(gcf,'units','norm','pos',[.2 .2 .6 .5],'paperpositionmode','auto');
% Drawing Constants
mark=char('o','o','o','o','o','o');
ccc=char('k','k','k','k','k','k');
hm=20; % how many posterior predictive samples
for dset=1:dsets
    subplot(2,3,dset);hold on;
    for i=1:hm
        r1 = ceil(rand*nchains);
        r2 = ceil(rand*nsamples);
        data = squeeze(samples.predpc(r1,r2,:,dset));
        ph = plot(1:listlength(dset),data(1:listlength(dset)),'k-');
        set(ph,'linewidth',2,'Color',[0.7 0.7 0.7]);
    end;
    % Draw the Probability Correct Data
    ph=plot([1:listlength(dset)],pc(dset,1:listlength(dset)),'ko');
    set(ph,'linewidth',.75,'markeredgecolor','k','markerfacecolor','w','markersize',5);
    switch dset,
        case {3,4}
            set(ph,'markersize',4);
        case {5,6}
            set(ph,'markersize',3);
    end;
    axis([0 41 0 1]);
    set(gca,'xtick',[0:10:40],'ytick',[0:.2:1],'box','on','fontsize',14);
    th=text(41,1,deblank(labs(dset,:)));
    set(th,'fontsize',14,'hor','right','vert','top');
end;
[ax,th]=suplabel('Serial Position','x');set(th,'fontsize',16);
set(th,'pos',get(th,'pos')+[0 .02 0]);
[ax,th]=suplabel('Probability Correct','y');set(th,'fontsize',16);
set(th,'pos',get(th,'pos')+[.02 0 0]);

% Joint Posterior
figure(2);clf;hold on;
% Drawing Constants
mark=char('o','s','d','^','v','<');
ccc=char('k','k','k','k','k','k');
epsc=.1;cx=[5:epsc:25];
cxe=[0+epsc/2:epsc:25-epsc/2];
epss=.1;sx=[7:epss:15];
sxe=[5+epss/2:epss:16-epss/2];
epst=.004;tx=[.45:epst:.65];
txe=[.4+epst/2:epst:.7-epst/2];
sc=2;joint=20;
% Draw Joint Posterior and Marginal
grid on;view([-25 30]);
for i=1:dsets
    ph=plot(-10,-10,'ko');
    set(ph,'marker',mark(i));
end;
for i=1:dsets
    keep=ceil(rand(joint,1)*nsamples);
    ph=plot3(samples.s(1,keep,i),samples.t(1,keep,i),samples.c(1,keep,i),'kp');
    set(ph,'marker',mark(i),'markersize',7,'markerfacecolor','w');
    ph=plot3(samples.s(1,keep,i),samples.t(1,keep,i),zeros(1,joint),'k.');
    ph=plot3(20*ones(1,joint),samples.t(1,keep,i),samples.c(1,keep,i),'k.');
    ph=plot3(samples.s(1,keep,i),.7*ones(1,joint),samples.c(1,keep,i),'k.');
    count=hist(reshape(samples.s,[],1),sx);
    count=count/sum(count);
    ph=plot3(sx,.3+sc*count,zeros(size(sx)),'k-');
    set(ph,'linewidth',1);
    count=hist(reshape(samples.t,[],1),tx);
    count=count/sum(count);
    ph=plot3(sc*20/.4*count,tx,zeros(size(tx)),'k-');
    set(ph,'linewidth',1);
    count=hist(reshape(samples.c,[],1),cx);
    count=count/sum(count);
    ph=plot3(sc*30/.4*count,.7*ones(size(cx)),cx,'k-');
    set(ph,'linewidth',1);
end;
set(gca,'fontsize',14,'gridline','--','linewidth',.25);
set(gca,'xlim',[0 20],'xtick',[0:5:20]);
set(gca,'ylim',[.3 .7],'ytick',[.3:.1:.7]);
set(gca,'zlim',[0 30],'ztick',[0:10:30]);
xlabel('Threshold Noise (s)','fontsize',16);
ylabel('Threshold (t)','fontsize',16);
zlabel('Distinctiveness (c)','fontsize',16);
[lh oh]=legend(labs,'location','bestoutside');
set(lh,'box','off');
