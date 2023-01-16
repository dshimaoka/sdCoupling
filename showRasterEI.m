function showRasterEI(spikeTimes, p, taxis, I)
%showRasterEI(spikeTimes, p)
%shows raster plot of all E and I cells


for icell = 1:p.Ne
    plot(spikeTimes{1}{icell},icell*ones(numel(spikeTimes{1}{icell}),1),'r.');
    hold on
end
for icell = 1:p.Ni
    plot(spikeTimes{2}{icell},(p.Ne+icell)*ones(numel(spikeTimes{2}{icell}),1),'b.');
    hold on
end
xlabel('time [ms]');
ylabel('cell ID (r:exc, b:inh)');
ylim([0 p.Ne+p.Ni]);

if nargin>=3
    xlim([taxis(1) taxis(end)]);
end

if nargin==4
    vbox(I.tstart, I.tend);
end

