cd('Z:\Shared\Daisuke\sandbox\Compte2003');
run('param_ds_single.m');
rmpath('C:\Users\dshi0006\git\dsbox\Stacked_Plot');
close all

%dt = 0.06; %ms Compte 2003
dt = 0.25; %ms Mainen 1996
% tspan = [0:dt:250];%ms
tspan_c = [0:dt:1000];%ms


%default parameters
doSingle = 0;
run('param.m');
p0 = p;

%% reconfigure Es > Ed connections in pLR fraction of excitatory neurons
plrPer = [0];%[0 10 20 40];

gEEPer = [13];%cFac * gEEPer/100 = 1
gIIPer = [20]; %[15 20 25 30];% 20 for SW?
gEIPer = [100];% [10 15 20 30 40 50]; %10< for E-I balance
gIEPer = [20 30 50 70 90];%[22];
gsdPer = [15];

p0.WEEd = zeros(p.Ne, p.Ne);

mCV = [];
mnrSpikes = [];
for nn = 1:numel(plrPer)
    %% define WEEd
    p = getLRconnectivity(p0, plrPer(nn)/100);
    
    %% assign non-zero gsd to long-range projection neurons
    LRidx = find(sum(p.WEEd,1)>0);
    p.gsd = zeros(p.Ne,1);
    
    for mm = 1:numel(gIEPer)
        for ll = 1:numel(gEEPer)
            for kk = 1:numel(gsdPer)
                    p.gsd(LRidx) = gsdPer(kk);
                for ii=1:numel(gIIPer)
                    for jj = 1:numel(gEIPer)
                        
                        p.gEEsAMPA = gEEPer(ll)/100*p0.gEEsAMPA;
                        p.gEEsNMDA = gEEPer(ll)/100*p0.gEEsNMDA;
                        p.gEEdAMPA = gEEPer(ll)/100*p0.gEEdAMPA;
                        p.gEEdNMDA = gEEPer(ll)/100*p0.gEEdNMDA;
                        
                        p.gIEs = gIEPer(mm)/100*p0.gIEs;
                        p.gIEd = gIEPer(mm)/100*p0.gIEd;
                        
                        p.gII = gIIPer(ii)/100*p0.gII;
                        p.gEIAMPA = gEIPer(jj)/100*p0.gEIAMPA;
                        p.gEINMDA = gEIPer(jj)/100*p0.gEINMDA;
                        %p.gsd = gsdPer(kk)/100*p0.gsd;
                        
                        suffix = ['_gsd0-15' '_rho500' '_pLR' num2str(plrPer(nn)) ...
                            '_gIIper' num2str(gIIPer(ii)) '_gEIper' num2str(gEIPer(jj)) ...
                            '_gEEper' num2str(gEEPer(ll)) '_gIEper' num2str(gIEPer(mm))];
                        
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
                        
                        
                        Netot = 20*p.Ne;
                        Nitot = 6*p.Ni;
                        
                        initialValue_c = zeros(Netot+Nitot,1);
                        initialValue_c(1:2*p.Ne,1) = -80*rand(2*p.Ne,1); %Vs, Vd
                        initialValue_c(Netot+1:Netot+p.Ni,1) = -80* rand(p.Ni,1); %Vi
                        initialValue_c(3*p.Ne+1:4*p.Ne,1) = p.Naeq; %Na+
                        initialValue_c(2*p.Ne+1:3*p.Ne,:) = 1e-3; %Ca2+ concentration
                        
                        [taxis,tcourse,spikeTime, cellID] = simulatorWClass(p,tspan_c,initialValue_c);
                        %[taxis,tcourse,spikeTime, cellID] = simulatorWClassWInput(p,tspan_c,initialValue_c, I);
                        
                        %% spike statistics
                        spikeTimes = [];
                        CV = [];
                        nrSpikes = [];
                        for icell = 1:p.Ne
                            spikeTimes{icell} = trace2Event(tcourse(:,icell)>0);
                            
                            if ~isempty(spikeTimes{icell})
                                ISI = diff(spikeTimes{icell}(:,1));
                                CV(icell) = sqrt(mean(ISI.^2)-(mean(ISI))^2)/mean(ISI);
                                nrSpikes(icell) = size(spikeTimes{icell},1);
                            else
                                CV(icell) = 0;
                                nrSpikes(icell) = 0;
                            end
                        end
                        
                        mCV(kk,ii,jj) = mean(CV);
                        mfrRate(kk,ii,jj) = mean(nrSpikes)/((taxis(end)-taxis(1))*1e-3);
                        
                        %% conductances
                        o = compte_ds_mainen(p, tcourse');
                        
                        if doSingle
                            figure(1);
                            plot(taxis,o.Vs,taxis,o.Vd)
                        end
                        icell = 1;
                        
                        figure(4);
                        %Intrinsic conductance for dendrite
                        [INad,condInt_Nad] = o.INad;
                        [IKm,condInt_Km] = o.IKm;
                        [IKCa,condInt_KCa] = o.IKCa;
                        [ICa,condInt_Ca] = o.ICa;
                        Ca = o.Ca;%expected o(uM)
                        
                        varNames_c = ["Ca","condInt_Nad (55)",...
                            "condInt_Km (-90)","condInt_Ca (140)","condInt_KCa (-90)"];
                        thisTable = array2timetable([...
                            Ca(icell,:);condInt_Nad(icell,:);condInt_Km(icell,:);condInt_Ca(icell,:);...
                            condInt_KCa(icell,:)]',...
                            'TimeStep',seconds(1e-3*dt),'variableNames',varNames_c);
                        set(gcf,'position',[0 0 1900 1000]);
                        stackedplot(thisTable);
                        screen2png(['conductance_dendrite' suffix]);close;
                        
                        figure(5);
                        %Intrinsic conductance for soma
                        Na = o.Na;
                        [INa,condInt_Na] = o.INa;
                        [IK,condInt_K] = o.IK;
                        [IA,condInt_A] = o.IA;
                        [IKs,condInt_KS] = o.IKS;
                        [IKNa,condInt_KNa] = o.IKNa;
                        
                        varNames_c = ["Na","condInt_Na",...
                            "condInt_K","condInt_A","condInt_Ks","condInt_KNa"];
                        thisTable = array2timetable([Na(icell,:);...
                            condInt_Na(icell,:);condInt_K(icell,:);condInt_A(icell,:);...
                            condInt_KS(icell,:);condInt_KNa(icell,:)]',...
                            'TimeStep',seconds(1e-3*dt),'variableNames',varNames_c);
                        set(gcf,'position',[0 0 1900 1000]);
                        stackedplot(thisTable);
                        screen2png(['conductance_soma' suffix]);close;
                        
                        % %% figure for intrinsic current for dendrite
                        % figure
                        % plot(taxis, -o.p.Ad.*o.INad, taxis, -o.p.Ad.*o.ICa,taxis, -o.p.Ad.*o.IKCa,taxis, ...
                        %     -o.p.Ad.*o.IKm )
                        % legend('INad','Ica','IKCa','IKm');
                        
                        
                        
                        %%
                        if ~doSingle
                            mVs = mean(tcourse(:,1:p.Ne)');
                            mVd = mean(tcourse(:,1+p.Ne:2*p.Ne)');
                            mVi = mean(tcourse(:,1+Netot:1+p.Ni+Netot)');
                        else
                            mVs = tcourse(:,1)';
                            mVd = tcourse(:,2)';
                            mVi = tcourse(:,1+Netot)';
                        end
                        %[pspec_s, axisPspec] = pmtm(mVs-mean(mVs),3,numel(taxis),1e3/dt);
                        varNames_e = ["mVs","mVd","mVi"];
                        thisTable = array2timetable([mVs' mVd' mVi'], ...
                            'TimeStep',seconds(1e-3*dt),'variableNames',varNames_e);
                        stackedplot(thisTable);
                        screen2png(['mV' suffix]);close;
                        
                        
                        if ~doSingle
                            figure(1); %raster plot of all excitatory neurons
                            set(gcf,'position',[0 0 1900 1000]);
                            plot(spikeTime(cellID<=p.Ne),cellID(cellID<=p.Ne),'r.');
                            hold on
                            plot(spikeTime(find((cellID<=2*p.Ne) .* (cellID>p.Ne)==1)),...
                                cellID(find((cellID<=2*p.Ne) .* (cellID>p.Ne)==1))-p.Ne,'b.');
                            legend('soma','dendrite');
                            xlabel('time [ms');
                            ylabel('exc cell ID');
                            screen2png(['spikes' suffix]);close;
                            
                            rmpath('C:\Users\dshi0006\git\dsbox\Stacked_Plot');
                            figure(2);
                            set(gcf,'position',[0 0 1900 1000]);
                            idx_e = icell:p.Ne:icell+14*p.Ne; %excitatory
                            varNames_e = ["Vs","Vd","Ca","Na","ssGABA","sdGABA","sdAMPA","sdNMDA",...
                                "ssAMPA","ssNMDA","h","hd","n","ha","mks"];
                            thisTable = array2timetable(tcourse(:,idx_e),'TimeStep',...
                                seconds(1e-3*dt),'variableNames',varNames_e);
                            stackedplot(thisTable);
                            screen2png(['exc' suffix]);close;
                            
                            
                            figure(3);
                            set(gcf,'position',[0 0 1900 1000]);
                            
                            idx_i = Netot+icell:p.Ni:Netot+6*p.Ni; %inhibitory
                            varNames_i = ["Vi","siAMPA","siNMDA","siGABA","hi","ni"];
                            thisTable = array2timetable(tcourse(:,idx_i),'TimeStep',...
                                seconds(1e-3*dt),'variableNames',varNames_i);
                            s=stackedplot(thisTable);
                            screen2png(['inh' suffix]);close;
                            
                            %% synaptic conductance
                            condSyn_exc_d=mean(o.condSyn_exc_d);
                            condSyn_exc_s=mean(o.condSyn_exc_s);
                            condSyn_exc_i=mean(o.condSyn_exc_i);
                            condSyn_inh_d=mean(o.condSyn_inh_d);
                            condSyn_inh_s=mean(o.condSyn_inh_s);
                            condSyn_inh_i=mean(o.condSyn_inh_i);
                            
                            subplot(311);
                            plot(taxis, condSyn_exc_s,taxis, condSyn_inh_s);
                            title('excitatory soma');
                            legend('exc synpatictic conductance','inh synpatictic conductance');
                            subplot(312);
                            plot(taxis, condSyn_exc_d,taxis, condSyn_inh_d);
                            title('excitatory dendrite');
                            subplot(313);
                            plot(taxis, condSyn_exc_i,taxis, condSyn_inh_i);
                            title('inhibitory neuron');
                            screen2png(['synpaticConductance' suffix]);close;
                        end
                        save(['param' suffix],'p');
                    end
                end
            end
        end
    end
end


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

