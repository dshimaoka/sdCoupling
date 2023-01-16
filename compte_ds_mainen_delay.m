classdef compte_ds_mainen_delay
    
    properties
        %excitatory soma (Compte 2003)
        Vs
        Na
        ssGABA
        ssAMPA %added
        ssNMDA %added
        h
        n
        ha
        mks
        
        %excitatory dendrite (Mainen 1996)
        Vd
        Ca
        sdGABA
        sdAMPA
        sdNMDA
        hd %added
        mca %added
        hca %added
        xi %added
        q %added
        md %added
        
        %inhibitory neuron (one compartment)
        Vi
        siAMPA
        siNMDA
        siGABA
        hi
        ni
        p
        
        %delay variables
        VsL
        VdL
        ViL
    end
    methods
        
        function o = compte_ds_mainen_delay(params,Var,VarL,varargin)
            o.p = params;
            o.Vs = Var(1:o.p.Ne,:); %membrane potential of excitatory neuron, soma mV
            o.Vd = Var(o.p.Ne+1:2*o.p.Ne,:); %membrane potential of excitatory neuron, dendrite mV
            o.Ca = Var(2*o.p.Ne+1:3*o.p.Ne,:);%intracellular calcium concentration, dendrite [mM]
            o.Na = Var(3*o.p.Ne+1:4*o.p.Ne,:);%intracellular sodium concentration, soma [mM]
            o.ssGABA = Var(4*o.p.Ne+1:5*o.p.Ne,:);%fraction of open GABA receptor channels on soma of excitatory neuron
            o.sdGABA = Var(5*o.p.Ne+1:6*o.p.Ne,:);%fraction of open GABA receptor channels on soma of excitatory neuron
            o.sdAMPA = Var(6*o.p.Ne+1:7*o.p.Ne,:);%fraction of open AMPA receptor channels on dendrite of excitatory neuron
            o.sdNMDA = Var(7*o.p.Ne+1:8*o.p.Ne,:);%fraction of open NMDA receptor channels on dendrite of excitatory neuron
            o.ssAMPA = Var(8*o.p.Ne+1:9*o.p.Ne,:);%fraction of open AMPA receptor channels on dendrite of excitatory neuron
            o.ssNMDA = Var(9*o.p.Ne+1:10*o.p.Ne,:);%fraction of open NMDA receptor channels on dendrite of excitatory neuron
            o.h = Var(10*o.p.Ne+1:11*o.p.Ne,:);%sodium channel for soma of excitatory neuron
            o.hd = Var(11*o.p.Ne+1:12*o.p.Ne,:);%sodium channel for dendrite of excitatory neuron
            o.n = Var(12*o.p.Ne+1:13*o.p.Ne,:);%delayed rectifier channel for soma of excitatory neuron
            o.ha = Var(13*o.p.Ne+1:14*o.p.Ne,:);%fast A-type sodium channel for soma of excitatory neuron
            o.mks = Var(14*o.p.Ne+1:15*o.p.Ne,:);%non-inactivating sodium channel for soma of excitatory neuron
            o.mca = Var(15*o.p.Ne+1:16*o.p.Ne,:);%calcium activation of excitatory dendrite Mainen 96
            o.hca = Var(16*o.p.Ne+1:17*o.p.Ne,:);%calcium inactivation of excitatory dendrite Mainen 96
            o.xi = Var(17*o.p.Ne+1:18*o.p.Ne,:);%KCa activation of excitatory dendrite Mainen 96
            o.q = Var(18*o.p.Ne+1:19*o.p.Ne,:);%KM activation of excitatory dendrite Mainen 96
            o.md = Var(19*o.p.Ne+1:20*o.p.Ne,:);%Nad activation of excitatory dendrite Mainen 96
            
            
            o.Vi = Var(o.p.Netot+1:o.p.Netot+o.p.Ni,:); %membrane potential of inhibitoroy neuron mV
            o.siAMPA  = Var(o.p.Netot+o.p.Ni+1:o.p.Netot+2*o.p.Ni,:); %fraction of open AMPA receptor channels on inhibitory neuron
            o.siNMDA = Var(o.p.Netot+2*o.p.Ni+1:o.p.Netot+3*o.p.Ni,:); %fraction of open NMDA receptor channels on inhibitory neuron
            o.siGABA = Var(o.p.Netot+3*o.p.Ni+1:o.p.Netot+4*o.p.Ni,:); %fraction of open GABA receptor channels on inhibitory neuron
            o.hi = Var(o.p.Netot+4*o.p.Ni+1:o.p.Netot+5*o.p.Ni,:); %sodium channel for inhibitory neuron
            o.ni = Var(o.p.Netot+5*o.p.Ni+1:o.p.Netot+6*o.p.Ni,:); %sodium channel for inhibitory neuron
            
            %delayed variables
            o.VsL = VarL(1:o.p.Ne,:); %membrane potential of excitatory neuron, soma mV
            o.VdL = VarL(o.p.Ne+1:2*o.p.Ne,:); %membrane potential of excitatory neuron, dendrite mV
            o.ViL = VarL(o.p.Netot+1:o.p.Netot+o.p.Ni,:);
        end
        
        
        function [currentDensity, conductance] = INa(o)
            %sodium current (soma)
            alpham = 0.1*(o.Vs+33)./(1 - exp(-(o.Vs+33)/10));
            betam = 4 * exp(-(o.Vs +53.7)/12);
            minf = alpham ./ (alpham + betam);
            currentDensity = o.p.gNa .* minf.^3 .* o.h .* (o.Vs - o.p.VNa); %[mS/cm^2]*[mV] = [1e-2 * A/m^2]
            conductance = 1e4* o.p.As * o.p.gNa .* minf.^3 .* o.h;
        end
        
        function [currentDensity, conductance] = INad(o)
            %sodium current (dendrite)
            %                 Zach Mainen created this particular model by adapting conductances
            %                 from lower temperature to run at higher temperature, and found it
            %                 necessary to reduce the temperature sensitivity of spike amplitude
            %                 and time course.  He accomplished this by increasing the net ionic
            %                 conductance through the heuristic of changing the standard HH
            %                 formula
            %                 g = gbar*product_of_gating_variables
            %                 to
            %                 g = tadj*gbar*product_of_gating_variables
            %                 where
            %                 tadj = q10^((celsius - temp)/10)
            %                 temp is the "reference temperature" (at which the gating variable
            %                 time constants were originally determined)
            %                 celsius is the "operating temperature"
            currentDensity = o.p.gNad .* o.p.phid .* o.md.^3 .* o.hd .* (o.Vd - o.p.VNad); %[mS/cm^2]*[mV] = [1e-2 * A/m^2]
            conductance = 1e4 * o.p.Ad *o.p.gNad .* o.p.phid .* o.md.^3 .* o.hd;
        end
        
        function [currentDensity, conductance] = INai(o)
            %sodium current (inhibitory neuron)
            alphami = 0.5*(o.Vi+35)./(1 - exp(-(o.Vi+35)/10));
            betami = 20 * exp(-(o.Vi + 60)/18);
            minfi = alphami ./ (alphami + betami);
            currentDensity = o.p.gNai .* minfi.^3 .* o.hi .* (o.Vi - o.p.VNai); %[1e-2 * A/m^2]
            conductance = 1e4*o.p.Ai * o.p.gNai .* minfi.^3 .* o.hi;
        end
        
        function [currentDensity, conductance] = IK(o)
            %potassium current (soma)
            currentDensity = o.p.gK .* o.n.^4 .* (o.Vs - o.p.VK); %[1e-2 * A/m^2]
            conductance = 1e4*o.p.As * o.p.gK .* o.n.^4;
        end
        
        function [currentDensity, conductance] = IKi(o)
            %potassium current (inhibitory neuron)
            currentDensity = o.p.gKi .* o.ni.^4 .* (o.Vi - o.p.VKi); %[1e-2 * A/m^2]
            conductance = 1e4 * o.p.Ai * o.p.gKi .* o.ni.^4;
        end
        
        function [currentDensity, conductance] = ILs(o)
            %leak current (soma)
            currentDensity = o.p.gLs .* (o.Vs - o.p.VL); %[1e-2 * A/m^2]
            conductance = 1e4*o.p.As * o.p.gLs;
        end
        
        function [currentDensity, conductance] = ILd(o)
            %leak current (dendrite)
            currentDensity = o.p.gLd .* (o.Vd - o.p.VLd); %[1e-2 * A/m^2]
            conductance = 1e4*o.p.Ad * o.p.gLd;
        end
        function [currentDensity, conductance] = ILi(o)
            %leak current (inhibitory neuron)
            currentDensity = o.p.gLi .* (o.Vi - o.p.VLi); %[1e-2 * A/m^2]
            conductance = 1e4 * o.p.Ai * o.p.gLi;
        end
        
        function [currentDensity, conductance] = IA(o)
            %fast A-type potassium current (soma)
            mainf = 1./(1+exp(-(o.Vs+50)/20));
            
            %should use main
            %should use ha NOT hainf
            currentDensity = o.p.gA .* mainf.^3 .* o.ha .* (o.Vs - o.p.VK);
            conductance = 1e4*o.p.As*o.p.gA .* mainf.^3 .* o.ha;
        end
        
        function [currentDensity, conductance] = IKS(o)
            % non-inactivating slow potassium current (soma)
            currentDensity = o.p.gKS.*o.mks.*(o.Vs - o.p.VK); %[1e-2 * A/m^2]
            conductance = 1e4 * o.p.As * o.p.gKS.*o.mks; %[nS]
        end
        
        function [currentDensity, conductance] = IKNa(o)
            % sodium dependent potassium current (soma)
            winf = 0.37./(1 + (38.7./o.Na).^3.5); %Birchoff 1998. Na in mM
            currentDensity = o.p.gKNa .* winf .* (o.Vs - o.p.VK); %[1e-2 * A/m^2]
            conductance = 1e4 * o.p.As * o.p.gKNa .* winf ;
        end
        
        function [currentDensity, conductance] = ICa(o)
            % high threshold calcium current (dendrite)
            %                 Zach Mainen created this particular model by adapting conductances
            %                 from lower temperature to run at higher temperature, and found it
            %                 necessary to reduce the temperature sensitivity of spike amplitude
            %                 and time course.  He accomplished this by increasing the net ionic
            %                 conductance through the heuristic of changing the standard HH
            %                 formula
            %                 g = gbar*product_of_gating_variables
            %                 to
            %                 g = tadj*gbar*product_of_gating_variables
            %                 where
            %                 tadj = q10^((celsius - temp)/10)
            %                 temp is the "reference temperature" (at which the gating variable
            %                 time constants were originally determined)
            %                 celsius is the "operating temperature"
            
            %mcainf = 1./(1 + exp(-(o.Vd + 20)/9));
            %currentDensity = o.p.gCa .* mcainf.^2 .* (o.Vd - o.p.VCa); %[1e-2 * A/m^2]%Compte 2003; Pinsky 1994
            %conductance = 1e4 * o.p.Ad * o.p.gCa .* mcainf.^2;
            
            currentDensity = o.p.gCa .* o.p.phid.* o.mca.^2 .* o.hca.* (o.Vd - o.p.VCa); %[1e-2 * A/m^2]%mainen 1993
            conductance = 1e4 * o.p.Ad * o.p.gCa .* o.p.phid .* o.mca.^2 .* o.hca;
        end
        
        function [currentDensity, conductance] = IKCa(o)
            % calcium dependent potassium current (dendrite)
            %Compte 2003:
            %                 currentDensity = o.p.gKCa .* o.Ca./(o.Ca + o.p.KD) .* (o.Vd - o.p.VK); %[1e-2 * A/m^2]
            %                 conductance = 1e4 * o.p.Ad * o.p.gKCa .* o.Ca./(o.Ca + o.p.KD);
            
            %                 Zach Mainen created this particular model by adapting conductances
            %                 from lower temperature to run at higher temperature, and found it
            %                 necessary to reduce the temperature sensitivity of spike amplitude
            %                 and time course.  He accomplished this by increasing the net ionic
            %                 conductance through the heuristic of changing the standard HH
            %                 formula
            %                 g = gbar*product_of_gating_variables
            %                 to
            %                 g = tadj*gbar*product_of_gating_variables
            %                 where
            %                 tadj = q10^((celsius - temp)/10)
            %                 temp is the "reference temperature" (at which the gating variable
            %                 time constants were originally determined)
            %                 celsius is the "operating temperature"
            currentDensity = o.p.gKCa .* o.p.phid.* o.xi .* (o.Vd - o.p.VKd); %[1e-2 * A/m^2] %Pinsky 1994
            conductance = 1e4 * o.p.Ad * o.p.gKCa .* o.p.phid.* o.xi; %Pinsky 1994
        end
        
        function [currentDensity, conductance] = IKm(o)
            % slow non-inactivating potassium current  (dendrite)
            %                 Zach Mainen created this particular model by adapting conductances
            %                 from lower temperature to run at higher temperature, and found it
            %                 necessary to reduce the temperature sensitivity of spike amplitude
            %                 and time course.  He accomplished this by increasing the net ionic
            %                 conductance through the heuristic of changing the standard HH
            %                 formula
            %                 g = gbar*product_of_gating_variables
            %                 to
            %                 g = tadj*gbar*product_of_gating_variables
            %                 where
            %                 tadj = q10^((celsius - temp)/10)
            %                 temp is the "reference temperature" (at which the gating variable
            %                 time constants were originally determined)
            %                 celsius is the "operating temperature"
            currentDensity = o.p.gKm .* o.p.phid.* o.q .* (o.Vd - o.p.VKd); %[1e-2 * A/m^2]
            conductance = 1e4 * o.p.Ad * o.p.gKm .* o.p.phid.* o.q; %[nS]
        end
        
        
        %% channel activation variables
        function dif = dh(o)
            %sodium channel excitatory soma
            alphah = 0.07 * exp(-(o.Vs + 50)/10);
            betah = 1./(1 + exp(-(o.Vs + 20)/10));
            dif = o.p.phih.*(alphah.*(1-o.h) - betah.*o.h);
        end
        
        function dif = dhi(o)
            %sodium channel inhibitory
            alphahi = 0.35 * exp(-(o.Vi + 58)/20);
            betahi = 5./(1 + exp(-(o.Vi + 28)/10));
            dif = o.p.phihi.*(alphahi.*(1-o.hi) - betahi.*o.hi);
        end
        
        function dif = dha(o)      % fast A-type potassium channel excitatory soma
            hainf = 1./(1+exp((o.Vs+80)/6));
            dif = o.p.phiha .* (hainf - o.ha) ./ o.p.tauha;
        end
        
        function dif = dmks(o)
            %excitatory soma
            mksinf = 1./(1 + exp(-(o.Vs+34)/6.5));
            taumks = 8./(exp(-(o.Vs+55)/30) + exp((o.Vs+55)/30));
            
            dif = o.p.phimks .* (mksinf - o.mks) ./ taumks;
        end
        
        function dif = dn(o) %persistent
            %potassium channel excitatory soma
            alphan = 0.01 * (o.Vs + 34)./(1 - exp(-(o.Vs+34)/10));
            betan = 0.125 * exp(-(o.Vs + 44)/25);
            
            dif = o.p.phin.*(alphan.*(1-o.n) - betan.*o.n);
        end
        
        function dif = dni(o) %persistent
            %potassium channel inhibitory neuron
            alphani = 0.05 * (o.Vi + 34)./(1 - exp(-(o.Vi+34)/10));
            betani = 0.625 * exp(-(o.Vi + 44)/80);
            
            dif = o.p.phini.*(alphani.*(1-o.ni) - betani.*o.ni);
        end
        
        function dif = dmca(o)
            %calcium activation @dendrite. Mainen 1996; Reuveni 1993
            %alphamca = 0.055*(o.Vd + 27)./(1 - exp(-(27+o.Vd)/3.8));
            tha = -27+o.p.Vdif;
            thb = -75+o.p.Vdif;
            alphamca = 0.209*o.efun(-(o.Vd - tha)/3.8);
            betamca = 0.94*exp(-(o.Vd - thb)/17);
            
            dif = o.p.phid*(alphamca.*(1-o.mca) - betamca.*o.mca);
        end
                
        function dif = dhca(o)
            %calcium inactivation @dendrite. Mainen 1996; Reuveni 1993
            tha = -13+o.p.Vdif;
            thb = -15+o.p.Vdif;
            alphahca = 4.57 * 1e-4 .* exp(-(o.Vd - tha)/50);
            betahca = 0.0065 ./ (1 + exp(-(o.Vd + thb)/28));
            dif = o.p.phid*(alphahca.*(1-o.hca) - betahca.*o.hca);
        end
        
        function dif = dxi(o)
            %KCa activation @dendrite Mainan 1996
            alphaxi = 0.01*o.Ca; %according to https://senselab.med.yale.edu/ModelDB/showmodel?model=2488&file=/cells/kca.mod#tabs-2
            betaxi = 0.02;
            xiinf = alphaxi./(alphaxi+betaxi);
            xitau = 1./o.p.phid./(alphaxi+betaxi);
            %dif = o.p.phid*(alphaxi.*(1-o.xi) - betaxi.*o.xi);
            dif = o.p.phid*(xiinf - o.xi)./xitau;
        end
        
        function dif = dq(o)
            %KM activation @ dendrite Mainen 1996
            %alphaq = 1e-4*(o.Vd+30)./(1 - exp(-(o.Vd+30)/9));
            %betaq = -1e-4*(o.Vd+30)./(1 - exp((o.Vd+30)/9)); %https://senselab.med.yale.edu/ModelDB/showmodel?model=2488&file=/cells/cad.mod#tabs-2
            %1.1e-4 rather than 1.1^(-4) produced similar qinf to Gutfreund 1995 (A1.3)
            qa=9; %mV inf slope
            tha = -30 + o.p.Vdif; %Vhalf for inf
            Ra   = 0.001;%	(/ms)		: max act rate  (slow)
            Rb   = 0.001;%	(/ms)		: max deact rate  (slow)
            alphaq = Ra * qa * o.efun(-(o.Vd - tha)/qa);
            betaq = Rb * qa * o.efun((o.Vd - tha)/qa);
            dif = o.p.phid*(alphaq .* (1-o.q) - betaq .* o.q);
        end
        
        function dif = dmd(o)
            %sodium activation @dendrite Mainen 1996
            vshift = 10;
            Vhalf = -35-vshift + o.p.Vdif;%https://senselab.med.yale.edu/ModelDB/showmodel?model=2488&file=/cells/na.mod#tabs-2
            %Vhalf = -25;%Mainen 1996 paper
            
            qa   = 9;%	(mV)		: act slope
            Ra   = 0.182;%	(/ms)		: open (v)
            Rb   = 0.124;%	(/ms)		: close (v)
            alpham = Ra * qa * o.efun(-(o.Vd - Vhalf)/qa);
            betam = Rb * qa * o.efun((o.Vd - Vhalf)/qa);
            taumd = 1./o.p.phid./(alpham+betam);
            minf = alpham./(alpham+betam);
            %dif = o.p.phid*(alpham .* (1-o.md) - betam .* o.md);
            dif = o.p.phid*(minf - o.md)./taumd;

        end
        
        function dif = dhd(o)
            %sodium channel excitatory dendrite
            %Mainen 1996; 1995:
            %                 alphah = 0.024 * (o.Vd + 40) ./ ( 1 - exp(-(o.Vd + 40)/5));
            %                 betah = -0.0091 * (o.Vd + 65) ./ ( 1 - exp((o.Vd + 65)/5));
            % hinf = 1./(1+exp(o.Vd+55)/6.2);
            
            vshift = 10;%-10??
            
            Rg   = 0.0091;%	(/ms)		: inact (v)
            Rd   = 0.024;%	(/ms)		: inact recov (v)
            thi1  = -50-vshift+o.p.Vdif;%	(mV)		: v 1/2 for inact
            thi2  = -75-vshift+o.p.Vdif;%	(mV)		: v 1/2 for inact
            qi   = 5;%	(mV)	        : inact tau slope
            thinf  = -55-vshift+o.p.Vdif;%-65;	(mV)		: inact inf slope	
            qinf  = 6.2;%	(mV)		: inact inf slope
            
            alphah = Rd * qi * o.efun(-(o.Vd - thi1)/qi);
            betah = Rg * qi * o.efun((o.Vd - thi2)/qi);
            tauh = 1./o.p.phid./(alphah+betah);
            hinf = 1./(1+exp((o.Vd-thinf)/qinf)); %o.Vd should be vm but NG

            dif = o.p.phid*(hinf - o.hd)./tauh;
            %dif = o.p.phid*(alphah .* (1-o.hd) - betah .* o.hd); %NG

        end
        
        %% synaptic currents
        function current = Isyns(o) %[1e-8 A]
            %synaptic current to excitatory soma
            current = 1e-4 * (o.p.gIEs .* o.ssGABA .* (o.Vs - o.p.VsynGABA) ...
                + o.p.gEEsAMPA .* o.ssAMPA .* (o.Vs - o.p.VsynAMPA) ...
                + o.p.gEEsNMDA .* o.ssNMDA .* (o.Vs - o.p.VsynNMDA));
        end
        
        function current = Isynd(o) %[1e-8 A]
            %synpatic current to excitatory dendrite
            current = 1e-4 * (o.p.gIEd .* o.sdGABA .* (o.Vd - o.p.VsynGABA) ...
                + o.p.gEEdAMPA .* o.sdAMPA .* (o.Vd - o.p.VsynAMPA) ...
                + o.p.gEEdNMDA .* o.sdNMDA .* (o.Vd - o.p.VsynNMDA));
        end
        
        function current = Isyni(o) %[1e-8 A]
            %synaptic current to inhibitory neuron
            current = 1e-4 * (o.p.gEIAMPA .* o.siAMPA .* (o.Vi - o.p.VsynAMPA) ...
                + o.p.gEINMDA .* o.siNMDA .* (o.Vi - o.p.VsynNMDA) ...
                + o.p.gII .* o.siGABA .* (o.Vi - o.p.VsynGABA));
        end
        
        %% neurotransmitter binding probability
        function dif = dssGABA(o)
            %                 dif = o.p.alphaGABA .* (o.p.WIE * o.ff(o.Vi)) .* (1-o.ssGABA) - o.ssGABA./o.p.tauGABA;% Tegner 2002
            dif = o.p.alphaGABA .* (o.p.WIEs * o.ff(o.ViL)) - o.ssGABA./o.p.tauGABA;
        end
        function dif = dsdGABA(o)
            %                 dif = o.p.alphaGABA .* (o.p.WII * o.ff(o.Vi)) .* (1 - o.siGABA) - o.siGABA./o.p.tauGABA;
            dif = o.p.alphaGABA .* (o.p.WIEd * o.ff(o.ViL))  - o.sdGABA./o.p.tauGABA;
        end
        function dif = dsdAMPA(o)
            %                  dif= o.p.alphaAMPA .* (o.p.WEE * o.ff(o.Vs)) .* (1-o.sdAMPA) - o.sdAMPA./o.p.tauAMPA;% Tegner 2002
            dif= o.p.alphaAMPA .* (o.p.WEEd * o.ff(o.VsL)) - o.sdAMPA./o.p.tauAMPA;%Pinsky 1994; Compte 2003
        end
        
        function dif = dsdNMDA(o)
            %                  dif = o.p.alphasNMDA .* o.xdNMDA .* (1 - o.sdNMDA) - o.sdNMDA./o.p.tausNMDA; %Tegner 2002; added xiNMDA to Compte 2003
            %                dif = o.p.alphasNMDA .* (1 - o.sdNMDA) - o.sdNMDA./o.p.tausNMDA; %written in Compte 2003
            dif = o.p.alphasNMDA .* (1 - o.sdNMDA) .* (o.p.WEEd * o.ff(o.VsL)) - o.sdNMDA./o.p.tausNMDA; %added network input to Compte 2003
        end
        
        function dif = dssAMPA(o)
            %                  dif= o.p.alphaAMPA .* (o.p.WEE * o.ff(o.VsL)) .* (1-o.sdAMPA) - o.sdAMPA./o.p.tauAMPA;% Tegner 2002
            dif= o.p.alphaAMPA .* (o.p.WEEs * o.ff(o.VsL)) - o.ssAMPA./o.p.tauAMPA;%Pinsky 1994; Compte 2003
        end
        
        function dif = dssNMDA(o)
            %                  dif = o.p.alphasNMDA .* o.xdNMDA .* (1 - o.sdNMDA) - o.sdNMDA./o.p.tausNMDA; %Tegner 2002; added xiNMDA to Compte 2003
            %                dif = o.p.alphasNMDA .* (1 - o.sdNMDA) - o.sdNMDA./o.p.tausNMDA; %written in Compte 2003
            dif = o.p.alphasNMDA .* (1 - o.ssNMDA) .* (o.p.WEEs * o.ff(o.VsL)) - o.ssNMDA./o.p.tausNMDA; %added network input to Compte 2003
        end
        
        function dif = dsiAMPA(o)
            %                 dif = o.p.alphaAMPA .* (o.p.WEI * o.ff(o.Vs)) .* (1 - o.siAMPA) - o.siAMPA./o.p.tauAMPA; %Tegner 2002
            dif = o.p.alphaAMPA .* (o.p.WEI * o.ff(o.VsL)) - o.siAMPA./o.p.tauAMPA; %Pinsky 1994; Compte 2003
        end
        
        function dif = dsiNMDA(o)
            %                  dif = o.p.alphasNMDA .* o.xiNMDA .* (1 - o.siNMDA) - o.siNMDA./o.p.tausNMDA; %Tegner 2002; added xiNMDA to Compte 2003
            %                dif = o.p.alphasNMDA .*  (1 - o.siNMDA) - o.siNMDA./o.p.tausNMDA; %written in Compte 2003 > produce constant siNMDA
            dif = o.p.alphasNMDA .*  (1 - o.siNMDA) .* (o.p.WEI * o.ff(o.VsL)) - o.siNMDA./o.p.tausNMDA; %written in Compte 2003
        end
        
        function dif = dsiGABA(o)
            %                 dif = o.p.alphaGABA .* (o.p.WII * o.ff(o.Vi)) .* (1 - o.siGABA) - o.siGABA./o.p.tauGABA;
            dif = o.p.alphaGABA .* (o.p.WII * o.ff(o.ViL))  - o.siGABA./o.p.tauGABA;
        end
        
        %% current balance
        function dif = dVs(o)%, extCurrent) %[mV]
%             if nargin<2
%                 extCurrent = 0;%0.25; %[nA] %Compte 2003
%             end
            %soma of excitatory neuron
            dif = (-o.p.As.*(o.ILs+o.INa+o.IK+o.IA+o.IKS+o.IKNa) ...
                - o.Isyns - 0.1*o.p.gsd .* (o.Vs-o.Vd) ...
                + 0.1*o.p.VsExtCurrent)./o.p.Cmd./o.p.As; %+ 0.025
        end
        
        function dif = dVd(o) %, extCurrent) %[mV]
%             if nargin<2
%                 extCurrent = 0.4; %[nA] Mainen 96 fig2b
%             end
            %dendrite of excitatory neuron
            dif = (-o.p.Ad.*(o.ILd+o.INad+o.ICa+o.IKCa+o.IKm) ...
                - 1.*o.Isynd - 0.1*o.p.gsd .*(o.Vd-o.Vs) ...
                + 0.1*o.p.VdExtCurrent)./o.p.Cm./o.p.Ad;
        end
        
        function dif = dVi(o)%, extCurrent) %[mV]
            %inhibitory neuron
%             if nargin<2
%                 extCurrent = 0; %[nA]
%             end
            dif = (-o.p.Ai.*(o.ILi+o.INai+o.IKi) - o.Isyni ...
                + 0.1*o.p.ViExtCurrent)./o.p.Cm./o.p.Ai;
        end
        
        function dif = dCa(o)
            %intracellular calcium concentration (dendrite)
            %dif = -1e-2 * o.p.alphaca .* o.p.Ad * o.ICa - (o.Ca - o.p.Cainf)./o.p.tauca; %[mM/ms] % Compte 2003
            
            F = 96485.3321233100184; %Faraday constant [C/m]
            depth = .1;	%(um) depth of shell https://senselab.med.yale.edu/ModelDB/showmodel?model=2488&file=/cells/cad.mod#tabs-2
            drive_channel = -1e5*o.ICa/2/F/depth;%*1e-2
            drive_channel(drive_channel<0) = 0; %cannot pump inward
            dif = drive_channel - (o.Ca - o.p.Cainf)/o.p.tauca;
        end
        
        function dif = dNa(o)
            %                  %intracellular sodium concentration (soma = dendrite?)
            %                    dif = - 10 * o.p.alphana .* (o.p.As .* o.INa + o.p.Ad .* o.INaP) - ...
            %                      o.p.Rpump .* (o.Na.^3 ./(o.Na.^3 + 15^3) - o.p.Naeq.^3 ./ (o.p.Naeq.^3 + 15^3)); %[mM/ms]
            
            % DS without dendritic contribution
            dif = - 10 * o.p.alphana .* (o.p.As .* o.INa) - ...
                o.p.Rpump .* (o.Na.^3 ./(o.Na.^3 + 15^3) - o.p.Naeq.^3 ./ (o.p.Naeq.^3 + 15^3)); %[mM/ms]
        end
        
        function allDif = dVar(o)
            allDif(1:o.p.Ne,1) = o.dVs; %membrane potential of excitatory neuron, soma
            allDif(o.p.Ne+1:2*o.p.Ne,1) = o.dVd; %membrane potential of excitatory neuron, dendrite
            allDif(2*o.p.Ne+1:3*o.p.Ne,1) = o.dCa;%intracellular calcium concentration, dendrite
            allDif(3*o.p.Ne+1:4*o.p.Ne,1) = o.dNa;%intracellular sodium concentration, soma
            allDif(4*o.p.Ne+1:5*o.p.Ne,1) = o.dssGABA;%fraction of open GABA receptor channels on soma of excitatory neuron
            allDif(5*o.p.Ne+1:6*o.p.Ne,1) = o.dsdGABA;%fraction of open GABA receptor channels on dendrite of excitatory neuron
            allDif(6*o.p.Ne+1:7*o.p.Ne,1) = o.dsdAMPA;%fraction of open AMPA receptor channels on dendrite of excitatory neuron
            allDif(7*o.p.Ne+1:8*o.p.Ne,1) = o.dsdNMDA;%fraction of open NMDA receptor channels on dendrite of excitatory neuron
            allDif(8*o.p.Ne+1:9*o.p.Ne,:) = o.dssAMPA;%fraction of open AMPA receptor channels on dendrite of excitatory neuron
            allDif(9*o.p.Ne+1:10*o.p.Ne,:) = o.dssNMDA;%fraction of open NMDA receptor channels on dendrite of excitatory neuron
            allDif(10*o.p.Ne+1:11*o.p.Ne,:) = o.dh;%sodium channel for soma of excitatory neuron
            allDif(11*o.p.Ne+1:12*o.p.Ne,:) = o.dhd;%sodium channel for dendrite of excitatory neuron
            allDif(12*o.p.Ne+1:13*o.p.Ne,:) = o.dn;%delayed rectifier channel for soma of excitatory neuron
            allDif(13*o.p.Ne+1:14*o.p.Ne,:) = o.dha;%fast A-type sodium channel for soma of excitatory neuron
            allDif(14*o.p.Ne+1:15*o.p.Ne,:) = o.dmks;%non-inactivating sodium channel for soma of excitatory neuron
            allDif(15*o.p.Ne+1:16*o.p.Ne,:) = o.dmca;%calcium activation of excitatory dendrite Mainen 96
            allDif(16*o.p.Ne+1:17*o.p.Ne,:) = o.dhca;%calcium inactivation of excitatory dendrite Mainen 96
            allDif(17*o.p.Ne+1:18*o.p.Ne,:) = o.dxi;%KCa activation of excitatory dendrite Mainen 96
            allDif(18*o.p.Ne+1:19*o.p.Ne,:) = o.dq;%KM activation of excitatory dendrite Mainen 96
            allDif(19*o.p.Ne+1:20*o.p.Ne,:) = o.dmd;%Nad activation of excitatory dendrite Mainen 96
            
            allDif(o.p.Netot+1:o.p.Netot+o.p.Ni,1) = o.dVi; %membrane potential of inhibitoroy neuron
            allDif(o.p.Netot+o.p.Ni+1:o.p.Netot+2*o.p.Ni,1) = o.dsiAMPA; %fraction of open AMPA receptor channels on inhibitory neuron
            allDif(o.p.Netot+2*o.p.Ni+1:o.p.Netot+3*o.p.Ni,1) = o.dsiNMDA; %fraction of open NMDA receptor channels on inhibitory neuron
            allDif(o.p.Netot+3*o.p.Ni+1:o.p.Netot+4*o.p.Ni,1) = o.dsiGABA; %fraction of open GABA receptor channels on inhibitory neuron
            allDif(o.p.Netot+4*o.p.Ni+1:o.p.Netot+5*o.p.Ni,1) = o.dhi; %sodium channel for inhibitory neuron
            allDif(o.p.Netot+5*o.p.Ni+1:o.p.Netot+6*o.p.Ni,1) = o.dni; %potassium channel for inhibitory neuron
        end
        
        function conductance = condSyn_exc_d(o) %correct?
            %total excitatory synaptic conductance of dendrite[nS]
            conductance = o.p.gEEdAMPA .* o.sdAMPA + o.p.gEEdNMDA .* o.sdNMDA;
        end
        
        function conductance = condSyn_exc_s(o) %correct?
            %total excitatory synaptic conductance of dendrite[nS]
            conductance = o.p.gEEsAMPA .* o.ssAMPA + o.p.gEEsNMDA .* o.ssNMDA;
        end
        function conductance = condSyn_exc_i(o) %correct?
            %total excitatory synaptic conductance of inhibitory neuron[nS]
            conductance = o.p.gEIAMPA .* o.siAMPA + o.p.gEINMDA .* o.siNMDA;
        end
        function conductance = condSyn_inh_s(o) %correct?
            %total inhibitory synaptic conductance of soma[nS]
            conductance = o.p.gIEs .* o.ssGABA;
        end
        function conductance = condSyn_inh_d(o) %correct?
            %total inhibitory synaptic conductance of dendrite[nS]
            conductance = o.p.gIEd .* o.sdGABA;
        end
       function conductance = condSyn_inh_i(o) %correct?
            %total inhibitory synaptic conductance of inhibitory neuron[nS]
            conductance = o.p.gII .* o.siGABA;
        end
    end %methods
    
    methods(Static)
        function spike = ff(Vpre)
            %works as a spike detector
            Vth = 0; %[mV] changed from +20 on 12/1/23
            spike = 1./(1 + exp(-(Vpre - Vth)/2));
        end
        
        function out = efun(z)
            %for mca
            %https://senselab.med.yale.edu/modeldb/ShowModel?model=2488&file=/cells/ca.mod#tabs-2
            if (abs(z) < 1e-6)
                out = 1 - z./2;
            else
                out = z./(exp(z) - 1);
            end
        end
    end %methods(static)
end %classdef