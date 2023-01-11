[general] = basicParams150728(1);
run('C:\Users\daisuke\Documents\MATLAB\GolombAmitai97\param_default.m');
dt = 0.1;

% dont use cells near boundary
xrange = round(0.4*p.N):p.N-round(0.4*p.N);

%sampled points for time traces
xsample = round(linspace(xrange(1), xrange(end), 5));


%% filtering in time (same as Analysis_GainEqualization)
general.fmax =  general.fmax*10; %hack..
[b, a] = cheby2(general.order,general.attenuation,...
    [general.fmin*2/(1000/dt) general.fmax*2/(1000/dt)]);%Riedner 2007


gampa = 0.5:0.1:1.5;%[0.1:0.1:2.5];
gkslow = 2.0;%[0.5:0.1:2.0];
tWindow = 900;
%taxis_peak = dt*(-tWindow/2:tWindow/2);

% target parameters for showing timecourse
ii_l = 1;
jj_l = 1;

ii_h = 11;
jj_h = 1;


speed = []; width = [];
for ii = 1:length(gampa)
    for jj = 1:length(gkslow)
        
        sname = sprintf('gampa%1d_gkslow%1d',round(10*gampa(ii)), round(10*gkslow(jj)));
        load(sname, 'tcourse','dt');
        
        W = [];
        LOCS = [];
        for xx = xrange
            tcourse_f(:,xx) =  filtfilt(b,a,double(tcourse(:,xx)-tcourse(1,xx))) + tcourse(1,xx);
            %tcourse_r = resample(tcourse(:,xx)-tcourse(end,xx),  round(general.srate),round(1000/dt));
            
            % compute peak time and width for each cell
            [pks_c,locs_c,wc] = findpeaks(tcourse_f(:,xx));
            
            if ~isempty(wc)
                [~,idx] = max(pks_c);
                W(xx) = dt * wc(idx); %[ms]
                LOCS(xx) = dt * locs_c(idx); %[ms]
            end
        end
        
        
        %% propagation speed as mean(diff(pks)/dx)
        speed(ii,jj) = 1000*(p.L/p.N) / mean(diff(LOCS(xrange))); %mm/s
        
        % duration as mean(w)
        width(ii,jj) = mean(W(xrange)); %ms
        
        if jj == jj_l && ii == ii_l
            tcourse_l = tcourse;
            tcourse_lf = tcourse_f;
        elseif jj == jj_h && ii == ii_h
            tcourse_h = tcourse;
            tcourse_hf = tcourse_f;
        end
    end
end


%% visualization

subplot(313)
[hax]=plotyy(gampa, speed(:,1),gampa, width(:,1));
xlabel('gAMPA');
set(hax,'tickdir','out');
% ylim(hax(1),[5 20]);
ylim(hax(2),[20 35]);
ylabel(hax(1),'Speed [mm/s]');
ylabel(hax(2),'Avg duration across neurons [ms]');

subplot(311)%low conductance regime
[~,maxTidx_l] = max(mean(tcourse_lf(:,xsample),2));
tidx_peak = maxTidx_l-tWindow/2:maxTidx_l+tWindow/2;
tidx_peak = tidx_peak(tidx_peak>0);
for kkk = 1:length(xsample)
    plot(dt*(tidx_peak-maxTidx_l), tcourse_lf(tidx_peak,xsample(kkk)), 'color',...
        [1-kkk/length(xsample) 0.25 kkk/length(xsample)]);
    hold on
end
%ylim([-25 60]);
xlim([-25 25]);
xlabel('Time from wave peak [ms]');
ylabel('Membrane potential [Vm]')
marginplot;

subplot(312)%low conductance regime
[~,maxTidx_h] = max(mean(tcourse_hf(:,xsample),2));
tidx_peak = maxTidx_h-tWindow/2:maxTidx_h+tWindow/2;
tidx_peak = tidx_peak(tidx_peak>0);
for kkk = 1:length(xsample)
    plot(dt*(tidx_peak-maxTidx_h), tcourse_hf(tidx_peak,xsample(kkk)), 'color',...
        [1-kkk/length(xsample) 0.25 kkk/length(xsample)]);
    hold on
end
%ylim([-25 60]);
xlim([-25 25]);
xlabel('Time from wave peak [ms]');
ylabel('Membrane potential [Vm]')
marginplot;
legend('','','','','');

savePaperFigure(gcf,'GolombAmitaiModel')
