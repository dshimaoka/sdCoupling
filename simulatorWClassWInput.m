function [taxis, tcourse, spikeTimes] = simulatorWClassWInput(p,tspan,initialValue, I)

%I.tstart
%I.tend
%I.VsExtCurrent
%I.VdExtCurrent
%I.ViExtCurrent

[tspanEv, polarity] = getTspanEvents(tspan, I.tstart, I.tend);

p0 = p;
p1 = p;
p1.VsExtCurrent = I.VsExtCurrent;
p1.VdExtCurrent = I.VdExtCurrent;
p1.ViExtCurrent = I.ViExtCurrent;

initialValue_c = initialValue;
taxis = [];
tcourse = [];
for it = 1:numel(tspanEv)
    
    tspan_c = tspanEv{it};
    
    switch polarity(it)
        case 0
            p = p0;
            [taxis_c,tcourse_c] = ode23(@Compte_ds_mainenModel_woInput, tspan_c, initialValue_c);
        case 1
            p = p1;
            [taxis_c,tcourse_c] = ode23(@Compte_ds_mainenModel_wInput, tspan_c, initialValue_c);
    end
    taxis = cat(1,taxis, taxis_c);
    tcourse = cat(1,tcourse, tcourse_c);
    initialValue_c = tcourse_c(end,:)';
end


if nargin>2
    for icell = 1:p.Ne
         cache = trace2Event(tcourse(:,icell)>0);
         spikeTimes{1}{icell} = taxis(cache(:,1)); %exc
     end
     for icell = p.Netot+1:p.Netot+p.Ni
         cache = trace2Event(tcourse(:,icell)>0);
         spikeTimes{2}{icell-p.Netot} = taxis(cache(:,1)); %inh
     end
end

    function dVar = Compte_ds_mainenModel_woInput(t, Var)
        if rem(t,1)==0
            disp(t);
        end
        o = compte_ds_mainen(p0, Var); %construct class
        dVar = o.dVar;
    end
    function dVar = Compte_ds_mainenModel_wInput(t, Var)
        if rem(t,1)==0
            disp(t);
        end
        o = compte_ds_mainen(p1, Var); %construct class
        dVar = o.dVar;
    end

    function [tspanEv, polarity] = getTspanEvents(tspan, tstart, tend)
        % [tspanEv, polarity] = getTspanEvents(tspan, tstart, tend)
        dt = mean(diff(tspan));
        tspan_trace = event2Trace(tspan,[tstart tend]);
        tspan_trace_i = tspan_trace==0;
        
        [evIdx] = sort(trace2Event(tspan_trace, (1:numel(tspan))));
        [evIdx_i] = sort(trace2Event(tspan_trace_i,(1:numel(tspan))));
        
        evIdx_all = union(evIdx,evIdx_i);
        polarity = tspan_trace(evIdx_all);
        
        tspanEv = cell(numel(evIdx_i)-1,1);
        for ii = 1:numel(evIdx_i)-1
            tspanEv{ii} = tspan(evIdx_i(ii)):dt:(tspan(evIdx_i(ii+1))-dt);
        end
        polarity = polarity(1:end-1);
    end

end
