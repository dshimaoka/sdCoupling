cd('Z:\Shared\Daisuke\sandbox\Compte2003');
run('param_ds_single.m');
rmpath('C:\Users\dshi0006\git\dsbox\Stacked_Plot');
close all

%dt = 0.06; %ms Compte 2003
dt = 0.25; %ms Mainen 1996
% tspan = [0:dt:250];%ms
tspan_c = [-1500:dt:200];%ms


%default parameters
doSingle = 1;
run('param_ds_mainen.m');
p0 = p;

%% variable parameters
gsdPer = sort([1 0:5:35]);%15; %[1 5 10 20 40 60 80 100]; %[%]
%15 is the best to evoke single-cell burst
%40 is too large
rhoPer = [1 100 200 300 400 500 600 700 800];

%% external curent
I.tstart = 0;
I.tend = 10;
I.VsExtCurrent = 1; %nA %0.5;Mainen96 (fig3)
I.VdExtCurrent = 1;
I.ViExtCurrent = 0;

nrSpikes = [];
for kk = 1:numel(rhoPer)
    for ii = 1:numel(gsdPer)
        suffix = ['_gsd' num2str(gsdPer(ii)) '_rho' num2str(rhoPer(kk)) ...
            '_VsCurrent' num2str(I.VsExtCurrent) '_VdCurrent' num2str(I.VdExtCurrent)];
        disp(suffix);
        
        p.gsd = gsdPer(ii)/100*p0.gsd;
        
        %other modified parameters
        p.Ad = rhoPer(kk)/100 * p0.As;
        p.Rpump = 1.5*p0.Rpump;
        p.gEEsAMPA = 0.7 * p0.gEEsAMPA;
        p.gEEsNMDA = 0.7 * p0.gEEsNMDA;
        p.gIEs = 1.1*p0.gIE; %necessary to reain E-I balance (Compte 2003 fig6)
        
        if doSingle
            % quench variability across cells
            VLSD = 0;
            p.gLi = 0.1025;
            p.VLi = -63.8;
        else
            VLSD = 12;
        end
        p.VL = -60.95+VLSD*randn(p.Ne,1); %necessary for spontaneous firing        
        
        Netot = 20*p.Ne;
        Nitot = 6*p.Ni;
        
        initialValue_c = zeros(Netot+Nitot,1);
        initialValue_c(1:2*p.Ne,1) = -80*rand(2*p.Ne,1); %Vs, Vd
        initialValue_c(Netot+1:Netot+p.Ni,1) = -80* rand(p.Ni,1); %Vi
        initialValue_c(3*p.Ne+1:4*p.Ne,1) = p.Naeq; %Na+
        initialValue_c(2*p.Ne+1:3*p.Ne,:) = 1e-3; %Ca2+ concentration
        
        %[taxis,tcourse,spikeTime, cellID] = simulatorWClass(p,tspan_c,initialValue_c);
        [taxis,tcourse,spikeTime, cellID] = simulatorWClassWInput(p,tspan_c,initialValue_c, I);
        
        %% count number of spikes
        nrSpikes(kk,ii) = size(trace2Event(tcourse(taxis>0,1)>0),1);
        %nrSpikes(kk,ii) = numel(find(tcourse(taxis>0,1) > 0)); %counts up and down
        
        %%
        mVs = tcourse(taxis>-50,1)';
        mVd = tcourse(taxis>-50,2)';
        
        taxis_c = taxis(taxis>-50);
        %[pspec_s, axisPspec] = pmtm(mVs-mean(mVs),3,numel(taxis),1e3/dt);
        varNames_e = ["mVs","mVd","mVi"];
        %         thisTable = array2timetable([mVs' mVd' mVi'], ...
        %             'TimeStep',seconds(1e-3*dt),'variableNames',varNames_e);
        %         stackedplot(thisTable);
        plot(taxis_c, mVs, taxis_c, mVd);
        vbox(I.tstart, I.tend, gca, [.9 .9 .9 .9]);
        legend(varNames_e);
        grid minor
        axis padded
        title(suffix(2:end));
        xlabel('time from current injection [ms]');
        ylabel('[mV]');
        screen2png(['mV' suffix]);close;
        
        %save(['param' suffix],'p');
    end
end
save(suffix(2:end),'p','nrSpikes','gsdPer','rhoPer','I');
imagesc(nrSpikes);
set(gca,'ytick',1:size(nrSpikes,2),'yticklabel',rhoPer,...
    'xtick',1:size(nrSpikes,1),'xticklabel',gsdPer);
ylabel('rho [%]');
xlabel('gsd [%]');

screen2png(gcf,'testBurstness_compte_ds_mainen');

