function [taxis, tcourse, spikeTime, cellID] = simulatorWClass(p,tspan,initialValue)

% [taxis,tcourse] = ode23(@GAmodel, tspan, initialValue);
% [taxis,tcourse] = ode45(@GAmodel, tspan, initialValue);
[taxis,tcourse] = ode23(@Compte_ds_mainenModel, tspan, initialValue);

if nargin>2
    [tidx, cellID] = find(tcourse > p.Vth);
    spikeTime = taxis(tidx);
end
  
    function dVar = CompteModel(t, Var)
        disp(t);
         o = compte(p, Var); %construct class
        dVar = o.dVar;
    end
    function dVar = Compte_dsModel(t, Var)
        disp(t);
         o = compte_ds(p, Var); %construct class
        dVar = o.dVar;
    end
 function dVar = Compte_ds_mainenModel(t, Var)
        disp(t);
         o = compte_ds_mainen(p, Var); %construct class
        dVar = o.dVar;
    end
end
