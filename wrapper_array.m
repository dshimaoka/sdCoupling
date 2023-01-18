%rmpath('C:\Users\dshi0006\git\dsbox\Stacked_Plot');
if ispc
    saveServer = ['X:' filesep 'Massive' filesep 'sdCoupling'];
    addpath('C:\Users\dshi0006\git\dsbox\');
else
    %saveServer = '/tmp/$(id -u)/gvfs/smb-share:server=storage.erc.monash.edu.au,share=shares/MNHS-dshi0006/Massive/sdCoupling';
    saveServer = '/home/dshi0006/tmpData';
    addpath('/home/dshi0006/git/dsbox');
end
close all

if getenv('SLURM_ARRAY_TASK_ID')
    idStr = getenv('SLURM_ARRAY_TASK_ID'); % get task ID string 
    fprintf('ID %s\n', idStr);
    arrayTaskID = str2num(idStr); % conv string to number
    pen = arrayTaskID;
else
    pen = 1;
end

%KNa for soma
gKNaPer = [100 60 40 20];

%GABA decay constant
tauGABAPer = [100 60 40 20];

%synaptic transmission delay
delay = 0; %[ms]

%soma-dendrite coupling
plrPer = [15];%5
gsdPer = [0 15];%15

%synaptic conductance
gEEPer = [13];%15%cFac * gEEPer/100 = 1
gIIPer = [20]; %[15 20 25 30];% 20 for SW?
gEIPer = [20];%[100];% 15 20 30 40 50]; %10< for E-I balance
gIEPer = [30];%[22];

sz = [numel(plrPer) numel(gsdPer) numel(gEEPer) numel(gIIPer) numel(gEIPer) ...
    numel(gIEPer) numel(gKNaPer) numel(tauGABAPer)];
%total jobs: prod(sz)
[plrPers,gsdPers,gEEPers,gIIPers,gEIPers,gIEPers,gKNaPers,tauGABAPers] = ...
    ind2sub(sz, pen);


%dt = 0.06; %ms Compte 2003
dt = 0.25; %ms Mainen 1996
tspan = [0:dt:2000];%ms


%default parameters
saveFig = 1;
lognWeight = 0;
doSingle = 0;
run('param.m');
p0 = p;

%% reconfigure Es > Ed connections in pLR fraction of excitatory neurons
p0.WEEd = zeros(p.Ne, p.Ne);
p = getLRconnectivity(p0, plrPer(plrPers)/100);


%% log-normal weight of WEEs
if lognWeight
    %find all connections
    theseConnections = find(p.WEEs(:,:)>0);
    %f = @(x)exp(-(x - tgtPosition).^2/2/sigma^2);
    sigma = 1;
    mu = log(0.2)+sigma^2;
    newWeights = lognrnd(mu, sigma, numel(theseConnections),1);
    p.WEEs(theseConnections) = newWeights;
end

%spike transmission delay
p.lags = delay;%[ms]

%% cell ID to visualize
thisCell = 1;



%% assign non-zero gsd to long-range projection neurons
LRidx = find(sum(p.WEEd,1)>0);
p.gsd = zeros(p.Ne,1);

p.gKNa = gKNaPer(gKNaPers)/100*p0.gKNa;
p.tauGABA = tauGABAPer(tauGABAPers)/100*p0.tauGABA;

p.gEEsAMPA = gEEPer(gEEPers)/100*p0.gEEsAMPA;
p.gEEsNMDA = gEEPer(gEEPers)/100*p0.gEEsNMDA;
p.gEEdAMPA = gEEPer(gEEPers)/100*p0.gEEdAMPA;
p.gEEdNMDA = gEEPer(gEEPers)/100*p0.gEEdNMDA;

p.gIEs = gIEPer(gIEPers)/100*p0.gIEs;
p.gIEd = gIEPer(gIEPers)/100*p0.gIEd;

p.gII = gIIPer(gIIPers)/100*p0.gII;
p.gEIAMPA = gEIPer(gEIPers)/100*p0.gEIAMPA;
p.gEINMDA = gEIPer(gEIPers)/100*p0.gEINMDA;
p.gsd = gsdPer(gsdPers)/100*p0.gsd;

suffix = ['_gsd' num2str(gsdPer(gsdPers)) '_rho500' '_pLR' num2str(plrPer(plrPers)) ...
    '_gIIper' num2str(gIIPer(gIIPers)) '_gEIper' num2str(gEIPer(gEIPers)) ...
    '_gEEper' num2str(gEEPer(gEEPers)) '_gIEper' num2str(gIEPer(gIEPers)) ...
    '_gKNaper' num2str(gKNaPer(gKNaPers)) '_tauGABAper' num2str(tauGABAPer(tauGABAPers)) ...
    '_lags' num2str(p.lags)];
if lognWeight
    suffix = [suffix '_lognWeight'];
end

saveDir = [saveServer filesep suffix(2:end)];
[status, msg, msgID] = mkdir(saveDir);

disp(status);
disp(msg);
disp(msgID);

disp(suffix);


%other modified parameters
p.Ad = 500/100 * p0.As;
p.Rpump = 1.5*p0.Rpump;
%p.gEEsAMPA = 0.7 * p0.gEEsAMPA;
%p.gEEsNMDA = 0.7 * p0.gEEsNMDA;
%p.gIEs = 1.1*p0.gIE; %necessary to reain E-I balance (Compte 2003 fig6)

if doSingle
    % quench variability across cells
    VLSD = 0;
    p.gLi = 0.1025;
    p.VLi = -63.8;
else
    VLSD = 12;
end
p.VL = -60.95+VLSD*randn(p.Ne,1); %necessary for spontaneous firing

%p.gNad = 1.5*p0.gNad;
%p.VCa = 145;%does not change spike amplitude
%p.gLd = 0;%;
% p.gNad = 10*p0.gNad;
% p.gCa = 10*p0.gCa;
% p.gKCa = 10*p0.gKCa;
% p.gKm = 10*p0.gKm;
% p.gKv = 10*p0.gKv;

init = @(t)getInitialValues(t,p);

[taxis, tcourse, spikeTimes] = ...
    simulatorWClass(p,tspan, init);
%[taxis,tcourse,spikeTime, cellID] = simulatorWClassWInput(p,tspan_c,initialValue_c, I);


%% spike statistics
for icell = 1:p.Ne
    
    if ~isempty(spikeTimes{1}{icell})
        ISI = diff(spikeTimes{1}{icell});
        CV(icell) = sqrt(mean(ISI.^2)-(mean(ISI))^2)/mean(ISI);
        nrSpikes(icell) = size(spikeTimes{1}{icell},1);
    else
        CV(icell) = 0;
        nrSpikes(icell) = 0;
    end
end

mCV = mean(CV);
mfrRate = mean(nrSpikes)/((taxis(end)-taxis(1))*1e-3); %[Hz]

%% obtain all variables
o = compte_ds_mainen(p, tcourse'); %compte_ds_mainen is fine


if saveFig
    %% Intrinsic conductance for dendrite
    figure('position',[0 0 1900 1000]);
    [INad,condInt_Nad] = o.INad;
    [IKm,condInt_Km] = o.IKm;
    [IKCa,condInt_KCa] = o.IKCa;
    [ICa,condInt_Ca] = o.ICa;
    Ca = o.Ca;%expected o(uM)
    
    varNames_c = ["Ca","condInt_Nad (55)",...
        "condInt_Km (-90)","condInt_Ca (140)","condInt_KCa (-90)"];
    thisTable = array2timetable([...
        Ca(thisCell,:);condInt_Nad(thisCell,:);condInt_Km(thisCell,:);condInt_Ca(thisCell,:);...
        condInt_KCa(thisCell,:)]',...
        'TimeStep',seconds(1e-3*dt),'variableNames',varNames_c);
    stackedplot(thisTable);
    saveas(gcf,[saveDir filesep 'conductance_dendrite' suffix '.png']);close;
    
    
    %% Intrinsic conductance for soma
    figure('position',[0 0 1900 1000]);
    Na = o.Na;
    [INa,condInt_Na] = o.INa;
    [IK,condInt_K] = o.IK;
    [IA,condInt_A] = o.IA;
    [IKs,condInt_KS] = o.IKS;
    [IKNa,condInt_KNa] = o.IKNa;
    
    varNames_c = ["Na","condInt_Na",...
        "condInt_K","condInt_A","condInt_Ks","condInt_KNa"];
    thisTable = array2timetable([Na(thisCell,:);...
        condInt_Na(thisCell,:);condInt_K(thisCell,:);condInt_A(thisCell,:);...
        condInt_KS(thisCell,:);condInt_KNa(thisCell,:)]',...
        'TimeStep',seconds(1e-3*dt),'variableNames',varNames_c);
    set(gcf,'position',[0 0 1900 1000]);
    stackedplot(thisTable);
    saveas(gcf, [saveDir filesep 'conductance_soma' suffix '.png']);close;
    
    % %% figure for intrinsic current for dendrite
    % figure
    % plot(taxis, -o.p.Ad.*o.INad, taxis, -o.p.Ad.*o.ICa,taxis, -o.p.Ad.*o.IKCa,taxis, ...
    %     -o.p.Ad.*o.IKm )
    % legend('INad','Ica','IKCa','IKm');
    
    
    
    %% mean membrane potential
    if ~doSingle
        mVs = mean(tcourse(:,1:p.Ne)');
        mVd = mean(tcourse(:,1+p.Ne:2*p.Ne)');
        mVi = mean(tcourse(:,1+p.Netot:1+p.Ni+p.Netot)');
    else
        mVs = tcourse(:,1)';
        mVd = tcourse(:,2)';
        mVi = tcourse(:,1+Netot)';
    end
    
    figure('position',[0 0 1900 1000]);
    %[pspec_s, axisPspec] = pmtm(mVs-mean(mVs),3,numel(taxis),1e3/dt);
    varNames_e = ["mVs","mVd","mVi"];
    thisTable = array2timetable([mVs' mVd' mVi'], ...
        'TimeStep',seconds(1e-3*dt),'variableNames',varNames_e);
    stackedplot(thisTable);
    saveas(gcf, [saveDir filesep 'mV' suffix '.png']);close;
    
    %% powerspectrum
    [pspec, axisPspec] = pmtm(mVs - mean(mVs), 3, numel(taxis), 1/dt*1e3);


    if ~doSingle
        %% raster plot of all neurons
        figure('position',[0 0 1900 1000]);
        showRasterEI(spikeTimes,p,taxis);
      
        saveas(gcf,[saveDir filesep 'spikes' suffix '.png']);close;
        
        
        %% variables for an excitatory neuron
        %rmpath('C:\Users\dshi0006\git\dsbox\Stacked_Plot');
        figure('position',[0 0 1900 1000]);
        idx_e = thisCell:p.Ne:thisCell+14*p.Ne; %excitatory
        varNames_e = ["Vs","Vd","Ca","Na","ssGABA","sdGABA","sdAMPA","sdNMDA",...
            "ssAMPA","ssNMDA","h","hd","n","ha","mks"];
        thisTable = array2timetable(tcourse(:,idx_e),'TimeStep',...
            seconds(1e-3*dt),'variableNames',varNames_e);
        stackedplot(thisTable);
        saveas(gcf, [saveDir filesep 'exc' suffix '.png']);close;
        
        
        %% variables for an inhibitory neuron
        figure('position',[0 0 1900 1000]);
        
        idx_i = p.Netot+thisCell:p.Ni:p.Netot+6*p.Ni; %inhibitory
        varNames_i = ["Vi","siAMPA","siNMDA","siGABA","hi","ni"];
        thisTable = array2timetable(tcourse(:,idx_i),'TimeStep',...
            seconds(1e-3*dt),'variableNames',varNames_i);
        s=stackedplot(thisTable);
        saveas(gcf,[saveDir filesep 'inh' suffix '.png']);close;
        
        
        %% synaptic conductance
        condSyn_exc_d=mean(o.condSyn_exc_d);
        condSyn_exc_s=mean(o.condSyn_exc_s);
        condSyn_exc_i=mean(o.condSyn_exc_i);
        condSyn_inh_d=mean(o.condSyn_inh_d);
        condSyn_inh_s=mean(o.condSyn_inh_s);
        condSyn_inh_i=mean(o.condSyn_inh_i);
        
        figure('position',[0 0 1900 1000]);
        subplot(311);
        yyaxis left; plot(taxis, condSyn_exc_s);
        yyaxis right; plot(taxis, condSyn_inh_s);
        title('excitatory soma');
        legend('exc synpatictic conductance','inh synpatictic conductance');
        subplot(312);
        yyaxis left; plot(taxis, condSyn_exc_d);
        yyaxis right; plot(taxis, condSyn_inh_d);
        title('excitatory dendrite');
        subplot(313);
        yyaxis left; plot(taxis, condSyn_exc_i);
        yyaxis right; plot(taxis, condSyn_inh_i);
        title('inhibitory neuron');
        saveas(gcf,[saveDir filesep 'synpaticConductance' suffix '.png']);close;
    end
end

save([saveDir filesep 'stats' suffix],'p','mCV','mfrRate','spikeTimes','tspan',...
    'axisPspec','pspec','mVs','mVd','mVi');




%% TO DO
% DONE introduce connectivity from dendrite to soma
% introduce mechanism for spike generation in dendrite (potassium and sodium conductance?)
% < DONE replace INaP (does not inactivate) with INa (does inactivate) in dendrite?
% < replace dendritic component with Pinsky 1994 or Mainen 96?
% change soma-dendrite coupling
% apply external perturbation

%factor F = faraday constant? > YES
%dimension consistency check
%calcium dynamics is too slow compared to Pinsky 94 fig2
%> replace w calcium dynamics used in Compte 2003
%check implementation of Mainen 1996 by others (neuron simulator)
%> https://senselab.med.yale.edu/modeldb/ShowModel?model=2488&file=/cells/demofig2.hoc#tabs-2
%gKv is absent in dendrite > YES
%should mutiply 1e-2 in 1st term of dCa? > NO
%value of gLd > INVERSE OF RESISTIVITY
%super fast subthreshold oscillation only to Km when gLd>0
%whatever the actual resting potential, dendrite does not fire (Voltage)
%>no voltage spike mechanism?
%<triple check Na and K channel kinetics, compare to soma
%<simulate without calcium-related channels

%DONE where to apply temp factor?
%DONE replace Nad channel according to neurodb?
%INad effect to Vd is too small
%scale factor for conductance? 10 or 0.1?
%scale factor for gLd
%revise parameters according to neurodb implementation
%DONE Nad (vshift?)
%DONE ca (mca, hca)
%DONE cad
%DONE kca
%DONE km

%Nad decay is too fast
%DONE make the function more readable (md and hd)
%DONE apply vshift to all dendritic components
%make the calcium spike larger
%introduce long-range connectivity Weed
%introduce p(proportion of cell area)
%increasing synI > reduce gEEs

%find parameter that achieves the following:
%s-d is only weakly coupled
%slow osc in soma but not in dendrite
%E-I balance is achieved

%explore rho so it can perform dendritic amplification
%> do this with transient pulse input
%
%apply pulse to network
%> to probe complexity
%
%assign rho values randomly?


%to achieve asynchronous state,
%assign log-normal gEEs(d) to each SYNAPSE? (Teramae 2012)

