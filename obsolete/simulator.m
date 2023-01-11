function [taxis, tcourse, spikeTime, cellID] = simulator(p,tspan,initialValue)
%TODO
% option to produce only spikes 

% [taxis,tcourse] = ode23(@GAmodel, tspan, initialValue);
% [taxis,tcourse] = ode45(@GAmodel, tspan, initialValue);
[taxis,tcourse] = ode23(@CompteModel, tspan, initialValue);

if nargin>2
    [tidx, cellID] = find(tcourse > p.Vth);
    spikeTime = taxis(tidx);
end
    function dVar = GAmodel(t,Var)
        
        disp(t);
        %parameters
        %N: num of neurons
        %C: 1uF/cm2
        %gNa, VNa the_m sig_m the_h sig_h the_ht sig_ht
        %gNp the_p sig_p
        %gKdr, VK the_n sig_n the_nt sig_nt
        %gKA the_a sig_a the_b sig_b tau_b
        %gKslow the_z sig_z tau_z
        %gL VL
        %the_s sig_s kf kr kt kv gAMP VGlu
        %kfN krN thefN sigfN gNMDA VGlu
        %L
        %lambda
        
        
        %% decompose input variables
        VV = Var(1:p.N,1);
        hh = Var(p.N+1:2*p.N,1);
        nn = Var(2*p.N+1:3*p.N,1);
        bb = Var(3*p.N+1:4*p.N,1);
        zz = Var(4*p.N+1:5*p.N,1);
        
        sAMPA = Var(5*p.N+1:6*p.N,1);
        TGlu = Var(6*p.N+1:7*p.N,1);
        sNMDA = Var(7*p.N+1:8*p.N,1);
        
        
        %AMPA current
        sinf = 1./(1+exp(-(VV-p.the_s)./p.sig_s)); %A24
        dTGlu = -p.kt * sinf .* TGlu + p.kv * (1 - TGlu); %A23
        dsAMPA = p.kf * TGlu .* sinf .* (1 - sAMPA) - p.kr * sAMPA; %A22
        IAMPA = p.gAMPA * (VV - p.VGlu) .* (p.Wtilde * sAMPA);%A21
        
        %NMDA current
        fNMDA = 1./(1 + exp(-(VV - p.the_fN)/p.sig_fN)); %A27 %WRONG?
        dsNMDA = p.kfN * TGlu .* sinf .* (1 - sNMDA) - p.krN * sNMDA; %A26
        INMDA = p.gNMDA * fNMDA .* (VV - p.VGlu) .* (p.Wtilde * sNMDA); %A25
        
        %leak current
        IL = p.gL*(VV - p.VL);
        
        %slow potassium current
        zinf = 1./(1 + exp(-(VV-p.the_z)./p.sig_z)); %A19
        dzz = (zinf - zz) ./ p.tau_z; %A18
        IKslow = p.gKslow .* zz .* (VV - p.VK); %A17
        
        %A-type potassium current
        ainf = 1./(1 + exp(-(VV-p.the_a)/p.sig_a)); %A15
        binf = 1./(1 + exp(-(VV-p.the_b)/p.sig_b)); %A16
        dbb = (binf - bb)/p.tau_b; %A14
        IKA = p.gKA .* ainf.^3 .* bb .* (VV - p.VK); %A13
        
        %delayed rectified potassium current
        ninf = 1./(1+exp(-(VV-p.the_n)/p.sig_n)); %A12
        taun = 0.37 + 1.85 ./ (1 + exp(-(VV - p.the_nt)./p.sig_nt)); %A11
        dnn = (ninf - nn)./taun; %A10
        IKdr = p.gKdr * nn.^4 .* (VV - p.VK); %A9
        
        %persistent sodium current
        pinf = 1./(1+exp(-(VV-p.the_p)./p.sig_p)); %A8
        INaP = p.gNaP * pinf .* (VV - p.VNa); %A7
        
        %sodium current
        tauh = 0.37 + 2.78 ./(1 + exp(-(VV-p.the_ht)/p.sig_ht)); %A6
        hinf = 1./(1 + exp(-(VV - p.the_h)./p.sig_h)); %A5
        minf = 1./(1 + exp(-(VV - p.the_m)./p.sig_m)); %A4
        dhh = (hinf - hh) ./ tauh; %A3
        INa = p.gNa * minf.^3 .* hh .* (VV-p.VNa); %A2
        
        %current balance
        dVV = (-INa - INaP - IKdr - IKA - IKslow - IL - IAMPA - INMDA)./p.C;
        
        %% output
        dVar(1:p.N,1) = dVV;
        dVar(p.N+1:2*p.N,1) = dhh;
        dVar(2*p.N+1:3*p.N,1) = dnn;
        dVar(3*p.N+1:4*p.N,1) = dbb;
        dVar(4*p.N+1:5*p.N,1) = dzz;
        
        dVar(5*p.N+1:6*p.N,1) = dsAMPA;
        dVar(6*p.N+1:7*p.N,1) = dTGlu;
        dVar(7*p.N+1:8*p.N,1) = dsNMDA;
        
        %     function wtilde(idx)
        %     wtilde(idx) = tanh(p.L/(2*p.lambda*p.N)) * exp(-abs(idx)*p.L/(p.lambda*p.N));
        % end
    end

    function dVar = CompteModel(t, Var)
                disp(t);

        %% decompose input variables
        Netot = 12*p.Ne;
        Nitot = 7*p.Ni;
        
        Vs = Var(1:p.Ne,1); %membrane potential of excitatory neuron, soma mV
        Vd = Var(p.Ne+1:2*p.Ne,1); %membrane potential of excitatory neuron, dendrite mV
        Ca = Var(2*p.Ne+1:3*p.Ne,1);%intracellular calcium concentration, dendrite [mM]
        Na = Var(3*p.Ne+1:4*p.Ne,1);%intracellular sodium concentration, soma [mM]
        ssGABA = Var(4*p.Ne+1:5*p.Ne,1);%fraction of open GABA receptor channels on soma of excitatory neuron
        sdAMPA = Var(5*p.Ne+1:6*p.Ne,1);%fraction of open AMPA receptor channels on dendrite of excitatory neuron
        xdNMDA = Var(6*p.Ne+1:7*p.Ne,1);%second messenger for NMDA receptor channels on dendrite of excitatory neuron
        sdNMDA = Var(7*p.Ne+1:8*p.Ne,1);%fraction of open NMDA receptor channels on dendrite of excitatory neuron
        h = Var(8*p.Ne+1:9*p.Ne,1);%sodium channel for soma of excitatory neuron 
        n = Var(9*p.Ne+1:10*p.Ne,1);%delayed rectifier channel for soma of excitatory neuron
        ha = Var(10*p.Ne+1:11*p.Ne,1);%fast A-type sodium channel for soma of excitatory neuron
        mks = Var(11*p.Ne+1:12*p.Ne,1);%non-inactivating sodium channel for soma of excitatory neuron
        
        Vi = Var(Netot+1:Netot+p.Ni,1); %membrane potential of inhibitoroy neuron mV
        siAMPA  = Var(Netot+p.Ni+1:Netot+2*p.Ni,1); %fraction of open AMPA receptor channels on inhibitory neuron
        xiNMDA = Var(Netot+2*p.Ni+1:Netot+3*p.Ni,1); %second messenger for NMDA receptor channels on inhibory neuron
        siNMDA = Var(Netot+3*p.Ni+1:Netot+4*p.Ni,1); %fraction of open NMDA receptor channels on inhibitory neuron
        siGABA = Var(Netot+4*p.Ni+1:Netot+5*p.Ni,1); %fraction of open GABA receptor channels on inhibitory neuron
        hi = Var(Netot+5*p.Ni+1:Netot+6*p.Ni,1); %sodium channel for inhibitory neuron
        ni = Var(Netot+6*p.Ni+1:Netot+7*p.Ni,1); %sodium channel for inhibitory neuron
        
        %sodium current (soma)
        alpham = 0.1*(Vs+33)./(1 - exp(-(Vs+33)/10));
        betam = 4 * exp(-(Vs +53.7)/12);
        minf = alpham ./ (alpham + betam);
        alphah = 0.07 * exp(-(Vs + 50)/10);
        betah = 1./(1 + exp(-(Vs + 20)/10));
        INa = p.gNa .* minf.^3 .* h .* (Vs - p.VNa); %[mS/cm^2]*[mV] = [1e-2 * A/m^2]
        
        %sodium current (inhibitory neuron)
        alphami = 0.5*(Vi+35)./(1 - exp(-(Vi+35)/10));
        betami = 20 * exp(-(Vi + 60)/18);
        minfi = alphami ./ (alphami + betami);
        alphahi = 0.35 * exp(-(Vi + 58)/20);
        betahi = 5./(1 + exp(-(Vi + 28)/10));
        INai = p.gNai .* minfi.^3 .* hi .* (Vi - p.VNai); %[1e-2 * A/m^2]
        
        %potassium current (soma)
        alphan = 0.01 * (Vs + 34)./(1 - exp(-(Vs+34)/10));
        betan = 0.125 * exp(-(Vs + 44)/25);
        IK = p.gK .* n.^4 .* (Vs - p.VK); %[1e-2 * A/m^2]
        
        %potassium current (inhibitory neuron)
        alphani = 0.05 * (Vi + 34)./(1 - exp(-(Vi+34)/10));
        betani = 0.625 * exp(-(Vi + 44)/80);
        IKi = p.gKi .* ni.^4 .* (Vi - p.VKi); %[1e-2 * A/m^2]
        
        %leak current (soma)
        IL = p.gL .* (Vs - p.VL); %[1e-2 * A/m^2]
        
        %leak current (inhibitory neuron)
        ILi = p.gLi .* (Vi - p.VLi); %[1e-2 * A/m^2]
        
        %fast A-type potassium current (soma)
        mainf = 1./(1+exp(-(Vs+50)/20));
        hainf = 1./(1+exp((Vs+80)/6));
        IA = p.gA .* mainf.^3 .* ha .* (Vs - p.VK);
        
        % non-inactivating slow potassium current (soma)
        mksinf = 1./(1 + exp(-(Vs+34)/6.5));
        taumks = 8./(exp(-(Vs+55)/30) + exp((Vs+55)/30));
        IKS = p.gKS.*mks.*(Vs - p.VK); %[1e-2 * A/m^2]
        
        % sodium dependent potassium current (soma)
        winf = 0.37./(1 + (38.7./Na).^3.5);
        IKNa = p.gKNa .* winf .* (Vs - p.VK); %[1e-2 * A/m^2]
        
        % high threshold calcium current (dendrite)
        mcainf = 1./(1 + exp(-(Vd + 20)/9));
        ICa = p.gCa .* mcainf.^2 .* (Vd - p.VCa); %[1e-2 * A/m^2]
        
        % calcium dependent potassium current (dendrite)
        IKCa = p.gKCa .* Ca./(Ca + p.KD) .* (Vd - p.VK); %[1e-2 * A/m^2]
        
        % non-inactivating (persistent) sodium current  (dendrite)
        mnapinf = 1./(1+exp(-(Vd+55.7)/7.7));
        INaP = p.gNaP .* mnapinf.^3 .* (Vd - p.VNa); %[1e-2 * A/m^2]
        
        % inward rectifier (activated by hyperpolarization) (dendrite)
        % non-inactivating current
        harinf = 1./(1+exp((Vd+75)/4));
        IAR = p.gAR.*harinf.*(Vd - p.VK); %[1e-2 * A/m^2]
        
        % synaptic currents (not explicitly written in paper)
        Isyns = -1e4 * p.gIE .* ssGABA .* (Vs - p.VsynGABA); %(excitatory soma) %[1e-8 A]
        Isynd = 1e4 * (p.gEEAMPA .* sdAMPA .* (Vd - p.VsynAMPA) ... %excitatory dendrite 
              + p.gEENMDA .* sdNMDA .* (Vd - p.VsynNMDA)); %[1e-8 A]
        Isyni = 1e4 * (p.gEIAMPA .* siAMPA .* (Vi - p.VsynAMPA) ... %inhibitory
            + p.gEINMDA .* siNMDA .* (Vi - p.VsynNMDA) ...
            - p.gII .* siGABA .* (Vi - p.VsynGABA)); %[1e-8 A]
        
        %         %spikes
        %         Se = (Vs>p.Vth);
        %         Si = (Vi>p.Vth);
        %
        %         %gating variables. connectivity part is inferred from Wang 1999b
        %         dssGABA = p.alphaGABA .* ff(Vs) .* (p.WIE * Si) - ssGABA./p.tauGABA;
        %
        %         dsdAMPA = p.alphaAMPA .* ff(Vd) .* (p.WEE * Se) - sdAMPA./p.tauAMPA;
        %         dxdNMDA = p.alphaxNMDA .* ff(Vd)  .* (p.WEE * Se) - xdNMDA./p.tauxNMDA;
        %         dsdNMDA = p.alphasNMDA .* xdNMDA .* (1 - sdNMDA) - sdNMDA./p.tausNMDA; %sum over connected neurons to first term
        %
        %         dsiAMPA = p.alphaAMPA .* ff(Vi) .* (p.WEI * Se) - siAMPA./p.tauAMPA; %sum over connected neurons to first term
        %         dxiNMDA = p.alphaxNMDA .* ff(Vi) .* (p.WEI * Se) - xiNMDA./p.tauxNMDA;
        %         dsiNMDA = p.alphasNMDA .* xiNMDA .* (1 - siNMDA) - siNMDA./p.tausNMDA; %added xiNMDA & sum over connected neurons to first term
        %
        %         dsiGABA = p.alphaGABA .* ff(Vi) .* (p.WII * Si) - siGABA./p.tauGABA;
        
         %gating variables. 
         dssGABA = p.alphaGABA .* (p.WIE * ff(Vi)) - ssGABA./p.tauGABA;
         
         dsdAMPA = p.alphaAMPA .* (p.WEE * ff(Vs)) - sdAMPA./p.tauAMPA;
         dxdNMDA = p.alphaxNMDA .* (p.WEE * ff(Vs)) - xdNMDA./p.tauxNMDA;
         dsdNMDA = p.alphasNMDA .* xdNMDA .* (1 - sdNMDA) - sdNMDA./p.tausNMDA; %sum over connected neurons to first term
         
         dsiAMPA = p.alphaAMPA .* (p.WEI * ff(Vs)) - siAMPA./p.tauAMPA; %sum over connected neurons to first term
         dxiNMDA = p.alphaxNMDA .* (p.WEI * ff(Vs)) - xiNMDA./p.tauxNMDA;
         dsiNMDA = p.alphasNMDA .* xiNMDA .* (1 - siNMDA) - siNMDA./p.tausNMDA; %added xiNMDA & sum over connected neurons to first term
         
         dsiGABA = p.alphaGABA .* (p.WII * ff(Vi)) - siGABA./p.tauGABA;
         
        %channel
        %sodium channel excitatory soma
        dh = p.phih.*(alphah.*(1-h) - betah.*h);
        
        %sodium channel inhibitory
        dhi = p.phih.*(alphahi.*(1-hi) - betahi.*hi); %same phih used in exc soma?
        
        %potassium channel excitatory soma
        dn = p.phin.*(alphan.*(1-n) - betan.*n);
        
        %potassium channel inhibitory neuron
        dni = p.phin.*(alphani.*(1-ni) - betani.*ni);
        
        % fast A-type potassium channel excitatory soma
        dha = p.phiha .* (hainf - ha) ./ p.tauha;
        
        %
        dmks = p.phimks .* (mksinf - mks) ./ taumks;
        
        %current balance
        %soma of excitatory neuron
        dVs = (-p.As.*(IL+INa+IK+IA+IKS+IKNa) - Isyns - 10*p.gsd.*(Vs-Vd))./p.Cm./p.As;
        %dendrite of excitatory neuron
        dVd = (-p.Ad.*(ICa+IKCa+INaP+IAR) - Isynd - 10*p.gsd.*(Vd-Vs))./p.Cm./p.Ad;
        %inhibitory neuron
        dVi = (-p.Ai.*(ILi+INai+IKi) - Isyni)./p.Cm./p.Ai;
        
        %intracellular calcium concentration (dendrite)
        %dCa = -1e4 * p.alphaca .* p.Ad * ICa - Ca./p.tauca; %[mM/ms]?
        dCa = -1e2 * p.alphaca .* p.Ad * ICa - Ca./p.tauca; %[mM/ms]
        
        %intracellular sodium concentration (soma = dendrite?)
        %         dNa = -1e4*p.alphana .* (p.As .* INa + p.Ad .* INaP) - ...
        %             p.Rpump .* (Na.^3 ./(Na.^3 + 15^3) - p.Naeq.^3 ./ (p.Naeq.^3 + 15^3)); %[mM/ms]
        dNa = - 10 * p.alphana .* (p.As .* INa + p.Ad .* INaP) - ...
            p.Rpump .* (Na.^3 ./(Na.^3 + 15^3) - p.Naeq.^3 ./ (p.Naeq.^3 + 15^3)); %[mM/ms]

        
        % output
         dVar(1:p.Ne,1) = dVs; %membrane potential of excitatory neuron, soma
        dVar(p.Ne+1:2*p.Ne,1) = dVd; %membrane potential of excitatory neuron, dendrite
        dVar(2*p.Ne+1:3*p.Ne,1) = dCa;%intracellular calcium concentration, dendrite
        dVar(3*p.Ne+1:4*p.Ne,1) = dNa;%intracellular sodium concentration, soma 
        dVar(4*p.Ne+1:5*p.Ne,1) = dssGABA;%fraction of open GABA receptor channels on soma of excitatory neuron
        dVar(5*p.Ne+1:6*p.Ne,1) = dsdAMPA;%fraction of open AMPA receptor channels on dendrite of excitatory neuron
        dVar(6*p.Ne+1:7*p.Ne,1) = dxdNMDA;%second messenger for NMDA receptor channels on dendrite of excitatory neuron
        dVar(7*p.Ne+1:8*p.Ne,1) = dsdNMDA;%fraction of open NMDA receptor channels on dendrite of excitatory neuron
        dVar(8*p.Ne+1:9*p.Ne,1) = dh;%sodium channel for soma of excitatory neuron
        dVar(9*p.Ne+1:10*p.Ne,1) = dn;%delayed rectifier channel for soma of excitatory neuron
        dVar(10*p.Ne+1:11*p.Ne,1) = dha;%fast A-type sodium channel for soma of excitatory neuron
        dVar(11*p.Ne+1:12*p.Ne,1) = dmks;%non-inactivating sodium channel for soma of excitatory neuron
        
        dVar(Netot+1:Netot+p.Ni,1) = dVi; %membrane potential of inhibitoroy neuron
        dVar(Netot+p.Ni+1:Netot+2*p.Ni,1) = dsiAMPA; %fraction of open AMPA receptor channels on inhibitory neuron
        dVar(Netot+2*p.Ni+1:Netot+3*p.Ni,1) = dxiNMDA; %second messenger for NMDA receptor channels on inhibory neuron
        dVar(Netot+3*p.Ni+1:Netot+4*p.Ni,1) = dsiNMDA; %fraction of open NMDA receptor channels on inhibitory neuron
        dVar(Netot+4*p.Ni+1:Netot+5*p.Ni,1) = dsiGABA; %fraction of open GABA receptor channels on inhibitory neuron
        dVar(Netot+5*p.Ni+1:Netot+6*p.Ni,1) = dhi; %sodium channel for inhibitory neuron
        dVar(Netot+6*p.Ni+1:Netot+7*p.Ni,1) = dni; %potassium channel for inhibitory neuron
    end

    function spike = ff(Vpre) 
        %works as a spike detector
        spike = 1./(1 + exp(-(Vpre - 20)/2));
    end
end
