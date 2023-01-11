cd('Z:\Shared\Daisuke\sandbox\Compte2003');

dt = 0.06; %ms
 % tspan = [0:dt:250];%ms
tspan_c = [-3000:dt:0];%ms

run('param_ds.m');
%run('param_single.m');
p.gKNa = 0.33; %fig5
p.Rpump = 1.5*p.Rpump;
p.WIE = cell2mat(arrayfun(@(x)getConnectivity(p.pos_e, x, p.sigma_i, p.Cie), p.pos_i, 'UniformOutput', false));
p.WEEs = cell2mat(arrayfun(@(x)getConnectivity(p.pos_e, x, p.sigma_e, p.Cees), p.pos_e, 'UniformOutput', false));
p.WEEd = cell2mat(arrayfun(@(x)getConnectivity(p.pos_e, x, p.sigma_e, p.Ceed), p.pos_e, 'UniformOutput', false));
p.WEI = cell2mat(arrayfun(@(x)getConnectivity(p.pos_i, x, p.sigma_e, p.Cei), p.pos_e, 'UniformOutput', false));
p.WII = cell2mat(arrayfun(@(x)getConnectivity(p.pos_i, x, p.sigma_i, p.Cii), p.pos_i, 'UniformOutput', false));


Netot = 13*p.Ne;
Nitot = 7*p.Ni;


initialValue_c = zeros(Netot+Nitot,1);
initialValue_c(1:2*p.Ne,1) = -80*rand(2*p.Ne,1); %Vs, Vd
initialValue_c(Netot+1:Netot+p.Ni,1) = -80* rand(p.Ni,1); %Vi
initialValue_c(3*p.Ne+1:4*p.Ne,1) = p.Naeq; %Na+
% initialValue_c(Netot+5*p.Ni+1:Netot+6*p.Ni,1) = 1; %hi
% load('initValue.mat');
% initialValue_c = initialValue;

[taxis,tcourse,spikeTime, cellID] = simulatorWClass(p,tspan_c,initialValue_c);


figure(1); %raster plot of all excitatory neurons
set(gcf,'position',[0 0 1900 1000]);
plot(spikeTime(cellID<=p.Ne),cellID(cellID<=p.Ne),'r.');
hold on
plot(spikeTime(find((cellID<=2*p.Ne) .* (cellID>p.Ne)==1)),cellID(find((cellID<=2*p.Ne) .* (cellID>p.Ne)==1))-p.Ne,'b.');
legend('soma','dendrite');
xlabel('time [ms');
ylabel('exc cell ID');
screen2png('spikes');close;


rmpath('C:\Users\dshi0006\git\dsbox\Stacked_Plot');
figure(2);
set(gcf,'position',[0 0 1900 1000]);
icell = 2;
idx_e = icell:p.Ne:icell+11*p.Ne; %excitatory 
varNames_e = ["Vs","Vd","Ca","Na","ssGABA","sdAMPA","xdNMDA","sdNMDA","h","n","ha","mks"];
thisTable = array2timetable(tcourse(:,idx_e),'TimeStep',seconds(1e-3*dt),'variableNames',varNames_e);
stackedplot(thisTable);
screen2png('exc');close;


figure(3);
set(gcf,'position',[0 0 1900 1000]);
icell = 1;

idx_i = Netot+icell:p.Ni:Netot+7*p.Ni; %inhibitory
varNames_i = ["Vi","siAMPA","xiNMDA","siNMDA","siGABA","hi","ni"];
thisTable = array2timetable(tcourse(:,idx_i),'TimeStep',seconds(1e-3*dt),'variableNames',varNames_i);
s=stackedplot(thisTable);
screen2png('inh');close;

%Ca rise is too fast?
%Na is building up?

%% mean across neurons - Not working yet
%idx = {1:p.Ne, p.Ne+1:2*p.Ne, Netot+1:Netot+p.Ni};

% meanFiringRate = [];
% for ii = 1:numel(idx)
%     meanFiringRate(:,ii) = 1/numel(idx{ii})*event2Trace(taxis, spikeTime(intersect(cellID, idx{ii})));
% end
% plot(taxis, meanFiringRate);


%% conductances
o = compte(p, tcourse');
%total excitatory synaptic conductance of soma[nS]
condSyn_exc = o.condSyn_exc;
%total inhibtory synaptic conductance of dendrite[nS]
condSyn_inh = o.condSyn_inh;
%Intrinsic conductance
[~,condInt_NaP] = o.INaP;
[~,condInt_KS] = o.IKS;
[~,condInt_AR] = o.IAR;
[~,condInt_KCa] = o.IKCa;
[~,condInt_KNa] = o.IKNa;

icell = 100;
varNames_c = ["condSyn_exc","condSyn_inh","condInt_NaP",...
    "condInt_KS","condInt_AR","condInt_KCa","condInt_KNa"];
thisTable = array2timetable([condSyn_exc(icell,:); condSyn_inh(icell,:);...
    condInt_NaP(icell,:);condInt_KS(icell,:);condInt_AR(icell,:);condInt_KCa(icell,:);condInt_KNa(icell,:)]',...
    'TimeStep',seconds(1e-3*dt),'variableNames',varNames_c);
figure(4);
set(gcf,'position',[0 0 1900 1000]);
stackedplot(thisTable);
screen2png('conductance');close;

%condInt NaP,KS too large?
%condInt AR & KCa looks ok
%condInt KNa is too slow because Na is too slow


%% TO DO
% try increasing gEE
% DONE check if NaP current from dendrite is critical
% change soma-dendrite coupling
% introduce connectivity from dendrite to soma
% apply external perturbation


%CIE=CEE=CII=CEI=50 - wave propagations w multiple speed
%NMDA only 1st order (xNMDA ignored). C= all 20 - inhibitory conductance accumulation
%CEE20 CEI5 CII5 CIE20 E-I balance achieved
%pGIE 4.56 silent for 3s
%p.gKNa = 0.33 0.4Hz
%pGIE 110%, p.gKNa = 0.33 0.35Hz
%pGEI 90%, p.gKNa = 0.33 0.35Hz
%pGEI 90%, pGIE110% silent for 3s
%all connections doubled, wave in each spike in burst then silent for 3s
%p.gKNa = 0.33 Rpump 0.027  0.43Hz
%p.gKNa = 0.33 Rpump 0.036  0.43Hz
%gKNa = 0.33, Rpump 0.027 No dendritic contribution to Na > 3times larger inh syn conductance. 0.35Hz