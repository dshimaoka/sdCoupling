
%% take iteration number
if ispc
    saveServer = ['X:' filesep 'Massive' filesep 'sdCoupling'];
    addpath(genpath('C:\Users\dshi0006\git\dsbox\'));
else
    %saveServer = '/tmp/$(id -u)/gvfs/smb-share:server=storage.erc.monash.edu.au,share=shares/MNHS-dshi0006/Massive/sdCoupling';
    saveServer = '/home/dshi0006/tmpData';
    addpath(genpath('/home/dshi0006/git/dsbox'));
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

nItr = 20; %total iteration per condition
dt = 0.25; %ms Mainen 1996
tspan = [0:dt:10000];%ms
ntgtCells = 100;%
extCurrent = 1.5; %nA
stimDur = 10; %ms
misi = 1000; %ms
jisi = 300; %ms actual isi ranges [misi-jisi misi+jisi]

%% load parameter
if pen < nItr
        suffix = '_gsd0_rho500_pLR15_gIIper20_gEIper20_gEEper13_gIEper30'; %synchronous
elseif pen < 2*nItr
        suffix = '_gsd15_rho500_pLR15_gIIper20_gEIper20_gEEper13_gIEper50';%asynchrnous 1
% else
%         suffix = '_gsd15_rho500_pLR15_gIIper20_gEIper20_gEEper10_gIEper30';%asynchrnous 2
end

pname = ['stats' suffix '.mat'];
load(pname, 'p');


%% set random number generator
rng('shuffle');
%otherwise the seed is identical across the jobs?


%% set initial values

initialValue_c = zeros(p.Netot+p.Nitot,1);
initialValue_c(1:2*p.Ne,1) = -80*rand(2*p.Ne,1); %Vs, Vd
initialValue_c(p.Netot+1:p.Netot+p.Ni,1) = -80* rand(p.Ni,1); %Vi
initialValue_c(3*p.Ne+1:4*p.Ne,1) = p.Naeq; %Na+
initialValue_c(2*p.Ne+1:3*p.Ne,:) = 1e-3; %Ca2+ concentration


%% set external input
tgtEcells = p.Ne/2-round(ntgtCells/2)+1:p.Ne/2+round(ntgtCells/2);
isis = 2*jisi*rand(round(tspan(end)/misi),1)+misi-jisi;
I.tstart = cumsum(isis);
I.tend = I.tstart + stimDur;
I.VsExtCurrent = zeros(p.Ne,1);
I.VsExtCurrent(tgtEcells) = extCurrent;
I.VdExtCurrent = 0;
I.ViExtCurrent = 0;


%% run the simulator
saveDir = [saveServer filesep suffix(2:end)];
[status, msg, msgID] = mkdir(saveDir);

[taxis,tcourse,spikeTimes] = simulatorWClassWInput(p,tspan,initialValue_c, I);


%% save the result
save([saveDir filesep 'stats' suffix '_I' num2str(1e3*extCurrent) 'pA_' num2str(pen)],...
    'p','spikeTimes','taxis','tspan','initialValue_c','I');


%% analysis
saveDir = 'X:\Massive\sdCoupling\20230117\gsd0_rho500_pLR15_gIIper20_gEIper20_gEEper13_gIEper30';
saveDir = 'X:\Massive\sdCoupling\20230117\gsd15_rho500_pLR15_gIIper20_gEIper20_gEEper13_gIEper50';
dd=dir(saveDir);
idx = strfind({dd.name},'.mat');

spikeTrace_all = [];
taxis_rs_all = [];
tstart_all = [];
for pen = 1:numel(dd)
    if ~isempty(idx{pen})
        load(fullfile(saveDir,dd(pen).name),...
            'p','spikeTimes','taxis','tspan','I');
        %         showRasterEI(spikeTimes,p,taxis,I);
        
        %% stimulus triggered avg
        dt_r = 5;%ms
        taxis_rs = tspan(1):dt_r:tspan(end);
        spikeTrace = [];
        for icell = 1:p.Ne
            spikeTrace(icell,:) = event2Trace(taxis_rs, spikeTimes{1}{icell});
        end
        
        spikeTrace_all = [spikeTrace_all spikeTrace];
        if isempty(taxis_rs_all)
            taxis_rs_all = taxis_rs;
            tstart_all = I.tstart;
        else
            taxis_rs_all = [taxis_rs_all max(taxis_rs_all)+dt_r+taxis_rs];
            tstart_all = [tstart_all max(taxis_rs_all)+dt_r+I.tstart];
        end
    end
end

[avgSpikeTrace, win_rs] = eventLockedAvg(single(spikeTrace_all), taxis_rs_all,...
    tstart_all, ones(numel(tstart_all),1),[-200 misi-jisi]);
figure('position',[0 0 1900 1000]);
subplot(3,1,[1 2]);imagesc(win_rs, 1:p.Ne, squeeze(avgSpikeTrace));
colormap(1-gray);
mcolorbar;
vline(0);
subplot(3,1,3);plot(win_rs, squeeze(avgSpikeTrace)');vline(0);
xlabel('time from stim onset [ms]');
saveas(gcf,[saveDir filesep 'sta' '.png']);close;



