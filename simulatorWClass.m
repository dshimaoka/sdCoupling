function [taxis, tcourse, spikeTimes] = simulatorWClass(p,tspan,initialValue)
%[taxis, tcourse, spikeTimes] = simulatorWClass(p,tspan,initialValue)
%
% 17/1/23 now compatible with p.delay

if p.lags==0
    [taxis,tcourse] = ode23(@Compte_ds_mainenModel, tspan, initialValue(1));
elseif p.lags>0
    options = ddeset('InitialStep',1);
    sol = dde23(@Compte_ds_mainen_delayModel, p.lags, initialValue, tspan, options);
    taxis = tspan';
    tcourse = deval(sol, taxis)';
end

if nargin>2
    %     [tidx, cellID] = find(tcourse > p.Vth);
    %     spikeTime = taxis(tidx);
     for icell = 1:p.Ne
         cache = trace2Event(tcourse(:,icell)>0);
         spikeTimes{1}{icell} = taxis(cache(:,1)); %exc
     end
     for icell = p.Netot+1:p.Netot+p.Ni
         cache = trace2Event(tcourse(:,icell)>0);
         spikeTimes{2}{icell-p.Netot} = taxis(cache(:,1)); %inh
     end
end
  
 function dVar = Compte_ds_mainen_delayModel(t, Var, VarL)
        if rem(t,1)==0
            disp(t);
        end
        
         o = compte_ds_mainen_delay(p, Var, VarL); %construct class
        dVar = o.dVar;
 end
function dVar = Compte_ds_mainenModel(t, Var)
        if rem(t,1)==0
            disp(t);
        end
        
        o = compte_ds_mainen(p, Var); %construct class
        dVar = o.dVar;
    end
end
