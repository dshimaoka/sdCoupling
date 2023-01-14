
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


%% load parameter
%suffix = '_gsd0_rho500_pLR15_gIIper20_gEIper20_gEEper13_gIEper30';
suffix = '_gsd15_rho500_pLR15_gIIper20_gEIper20_gEEper13_gIEper50';
%suffix = '_gsd15_rho500_pLR15_gIIper20_gEIper20_gEEper10_gIEper30';

pname = ['stats' suffix '.mat'];
load(pname, 'p');

%% set initial values
dt = 0.25; %ms Mainen 1996
tspan = [0:dt:10000];%ms

initialValue_c = zeros(p.Netot+p.Nitot,1);
initialValue_c(1:2*p.Ne,1) = -80*rand(2*p.Ne,1); %Vs, Vd
initialValue_c(p.Netot+1:p.Netot+p.Ni,1) = -80* rand(p.Ni,1); %Vi
initialValue_c(3*p.Ne+1:4*p.Ne,1) = p.Naeq; %Na+
initialValue_c(2*p.Ne+1:3*p.Ne,:) = 1e-3; %Ca2+ concentration


%% set external input
tgtEcells = p.Ne/2-50+1:p.Ne/2+50;
extCurrent = 1; %nA
stimDur = 10; %ms
misi = 800; %ms
jisi = 300; %ms actual isi ranges [misi-jisi misi+jisi]

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
    'p','spikeTimes','taxis','tspan');
