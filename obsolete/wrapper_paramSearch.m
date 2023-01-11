cd('Z:\Shared\Daisuke\sandbox\Compte2003');
run('param_ds_single.m');
rmpath('C:\Users\dshi0006\git\dsbox\Stacked_Plot');

dt = 0.06; %ms
% tspan = [0:dt:250];%ms
tspan_c = [-4000:dt:0];%ms

%default parameters
%run('param_ds_single.m');
run('param_ds.m');
p0 = p;

%modified parameters
gKNa = [0.33 1];
p.Rpump = 1.5*p.Rpump;
rgsd = [0.25 1];
rgEEs = [0.6 0.7 0.8];
p.gLs = 0.0667;%067*randn(p.Ne,1);%mS/cm^2
%p.gLd = 0.667;%67*randn(p.Ne,1);%mS/cm^2
p.VL = -60.95+12*randn(p.Ne,1); %threshold:12 for soma


for igkna = 1:numel(gKNa)
    for irgEE = 1:numel(rgEEs)
        
        p.gsd = 0;
        p.gEEsAMPA = rgEEs(irgEE) * p0.gEEsAMPA;
        p.gEEsNMDA = rgEEs(irgEE) * p0.gEEsNMDA;
        p.gKNa = gKNa(igkna) * p0.gKNa;
        
        suffix = ['_gKNa' num2str(100*gKNa(igkna)) '_gEEs' num2str(100*rgEEs(irgEE))];
        
        Netot = 14*p.Ne;
        Nitot = 7*p.Ni;
        
        initialValue_c = zeros(Netot+Nitot,1);
        initialValue_c(1:2*p.Ne,1) = -80*rand(2*p.Ne,1); %Vs, Vd
        initialValue_c(Netot+1:Netot+p.Ni,1) = -80* rand(p.Ni,1); %Vi
        initialValue_c(3*p.Ne+1:4*p.Ne,1) = p.Naeq; %Na+
        
        [taxis,tcourse,spikeTime, cellID] = simulatorWClass(p,tspan_c,initialValue_c);
        
         %% conductances
         icell = 1;
        o = compte_ds(p, tcourse');
        plot(taxis,o.Vs,taxis,o.Vd)
        
        %total excitatory synaptic conductance of soma[nS]
        condSyn_exc = o.condSyn_exc_s;
        %total inhibtory synaptic conductance of dendrite[nS]
        condSyn_inh = o.condSyn_inh;
        %Intrinsic conductance
        [~,condInt_NaP] = o.INaP;
        [~,condInt_KS] = o.IKS;
        [~,condInt_AR] = o.IAR;
        [~,condInt_KCa] = o.IKCa;
        [~,condInt_KNa] = o.IKNa;
        
        varNames_c = ["condSyn_exc","condSyn_inh","condInt_NaP",...
            "condInt_KS","condInt_AR","condInt_KCa","condInt_KNa"];
        thisTable = array2timetable([condSyn_exc(icell,:); condSyn_inh(icell,:);...
            condInt_NaP(icell,:);condInt_KS(icell,:);condInt_AR(icell,:);condInt_KCa(icell,:);condInt_KNa(icell,:)]',...
            'TimeStep',seconds(1e-3*dt),'variableNames',varNames_c);
        figure(4);
        set(gcf,'position',[0 0 1900 1000]);
        stackedplot(thisTable);
        screen2png(['conductance' suffix]);close;
        
        
        mVs = mean(tcourse(:,1:p.Ne)');
        mVd = mean(tcourse(:,1+p.Ne:2*p.Ne)');
        mVi = mean(tcourse(:,1+Netot:1+p.Ni+Netot)');
        %[pspec_s, axisPspec] = pmtm(mVs-mean(mVs),3,numel(taxis),1e3/dt);
        varNames_e = ["mVs","mVd","mVi"];
        thisTable = array2timetable([mVs' mVd' mVi'],'TimeStep',seconds(1e-3*dt),'variableNames',varNames_e);
        stackedplot(thisTable);
        screen2png(['mV' suffix]);close;
        
        
        figure(1); %raster plot of all excitatory neurons
        set(gcf,'position',[0 0 1900 1000]);
        plot(spikeTime(cellID<=p.Ne),cellID(cellID<=p.Ne),'r.');
        hold on
        plot(spikeTime(find((cellID<=2*p.Ne) .* (cellID>p.Ne)==1)),cellID(find((cellID<=2*p.Ne) .* (cellID>p.Ne)==1))-p.Ne,'b.');
        legend('soma','dendrite');
        xlabel('time [ms');
        ylabel('exc cell ID');
        screen2png(['spikes' suffix]);close;
        

        icell = 10;

        rmpath('C:\Users\dshi0006\git\dsbox\Stacked_Plot');
        figure(2);
        set(gcf,'position',[0 0 1900 1000]);
        idx_e = icell:p.Ne:icell+13*p.Ne; %excitatory
        varNames_e = ["Vs","Vd","Ca","Na","ssGABA","sdAMPA","sdNMDA","ssAMPA","ssNMDA","h","hd","n","ha","mks"];
        thisTable = array2timetable(tcourse(:,idx_e),'TimeStep',seconds(1e-3*dt),'variableNames',varNames_e);
        stackedplot(thisTable);
        screen2png(['exc' suffix]);close;
        
        
        figure(3);
        set(gcf,'position',[0 0 1900 1000]);
        icell = 1;
        
        idx_i = Netot+icell:p.Ni:Netot+7*p.Ni; %inhibitory
        varNames_i = ["Vi","siAMPA","xiNMDA","siNMDA","siGABA","hi","ni"];
        thisTable = array2timetable(tcourse(:,idx_i),'TimeStep',seconds(1e-3*dt),'variableNames',varNames_i);
        s=stackedplot(thisTable);
        screen2png(['inh' suffix]);close;
        
       
    end
end
    
    
%% TO DO
% DONE introduce connectivity from dendrite to soma
% introduce mechanism for spike generation in dendrite (potassium and sodium conductance?)
% < DONE replace INaP (does not inactivate) with INa (does inactivate) in dendrite?
% < replace dendritic component with Pinsky 1994 or Mainen 96?
% change soma-dendrite coupling
% apply external perturbation
    
    
