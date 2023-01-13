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


%soma-dendrite coupling
plrPer = [0];%5
gsdPer = [0];%15

%synaptic conductance
gEEPer = [10 13 15];%15%cFac * gEEPer/100 = 1
gIIPer = [15 20 25]; %[15 20 25 30];% 20 for SW?
gEIPer = [15 20 25];%[100];% 15 20 30 40 50]; %10< for E-I balance
gIEPer = [15 20 25];% 30 50 70 90];%[22];

sz = [numel(plrPer) numel(gsdPer) numel(gEEPer) numel(gIIPer) numel(gEIPer) numel(gIEPer)];
%total jobs: prod(sz)



for pen = 1:81
    
    [plrPers,gsdPers,gEEPers,gIIPers,gEIPers,gIEPers] = ind2sub(sz, pen);
    
    
    suffix = ['_gsd' num2str(gsdPer(gsdPers)) '_rho500' '_pLR' num2str(plrPer(plrPers)) ...
        '_gIIper' num2str(gIIPer(gIIPers)) '_gEIper' num2str(gEIPer(gEIPers)) ...
        '_gEEper' num2str(gEEPer(gEEPers)) '_gIEper' num2str(gIEPer(gIEPers))];
    
    try

        saveDir = [saveServer filesep  '20230112_paramSearch' filesep suffix(2:end) ]; 
        load([saveDir filesep 'stats' suffix],'p','mCV','mfrRate','spikeTimes','tspan',...
            'axisPspec','pspec','mVs','mVd','mVi');
        mCVp(plrPers,gsdPers,gEEPers,gIIPers,gEIPers,gIEPers) = mCV;
        mfrRatep(plrPers,gsdPers,gEEPers,gIIPers,gEIPers,gIEPers) = mfrRate;
        mVsp(:,plrPers,gsdPers,gEEPers,gIIPers,gEIPers,gIEPers) = mVs;
        mpspecp(:,plrPers,gsdPers,gEEPers,gIIPers,gEIPers,gIEPers) = pspec;
        
    catch err
        disp([suffix ' NOT FOUND']);
    end
end

%% firing rate
for ii = 1:3
    for jj = 1:3
        subplot(3,3,ii+3*(jj-1));
        imagesc(squeeze(mfrRatep(1,1,:,ii,jj,:)));
        title(['giiPers' num2str(gIIPer(ii)) ', geiPers' num2str(gEIPer(jj))]);
        set(gca,'ytick',1:numel(gEEPer),'ytickLabel',gEEPer);ylabel('gEEper');
        set(gca,'xtick',1:numel(gIEPer),'xtickLabel',gIEPer);xlabel('gIEper');
        caxis(prctile(mfrRatep(:),[1 99]));
    end
end
mcolorbar;
saveas(gcf,'mfrRate.png');

%% CV
for ii = 1:3
    for jj = 1:3
        subplot(3,3,ii+3*(jj-1));
        imagesc(squeeze(mCVp(1,1,:,ii,jj,:)));
        title(['giiPers' num2str(gIIPer(ii)) ', geiPers' num2str(gEIPer(jj))]);
        set(gca,'ytick',1:numel(gEEPer),'ytickLabel',gEEPer);ylabel('gEEper');
        set(gca,'xtick',1:numel(gIEPer),'xtickLabel',gIEPer);xlabel('gIEper');
        caxis(prctile(mCVp(:),[1 99]));
    end
end
mcolorbar;
saveas(gcf,'mCV.png');
