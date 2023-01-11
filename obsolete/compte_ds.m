classdef compte_ds
    
    properties
        Vs
        Vd
        Ca
        Na
        ssGABA
        sdAMPA
        sdNMDA
        ssAMPA %added
        ssNMDA %added
        h
        hd %added
        n
        ha
        mks
        Vi
        siAMPA
        xiNMDA
        siNMDA
        siGABA
        hi
        ni
        p
    end
        methods
            
            function o = compte_ds(params,Var,varargin)
                o.p = params;
                o.Vs = Var(1:o.p.Ne,:); %membrane potential of excitatory neuron, soma mV
                o.Vd = Var(o.p.Ne+1:2*o.p.Ne,:); %membrane potential of excitatory neuron, dendrite mV
                o.Ca = Var(2*o.p.Ne+1:3*o.p.Ne,:);%intracellular calcium concentration, dendrite [mM]
                o.Na = Var(3*o.p.Ne+1:4*o.p.Ne,:);%intracellular sodium concentration, soma [mM]
                o.ssGABA = Var(4*o.p.Ne+1:5*o.p.Ne,:);%fraction of open GABA receptor channels on soma of excitatory neuron
                o.sdAMPA = Var(5*o.p.Ne+1:6*o.p.Ne,:);%fraction of open AMPA receptor channels on dendrite of excitatory neuron
                o.sdNMDA = Var(6*o.p.Ne+1:7*o.p.Ne,:);%fraction of open NMDA receptor channels on dendrite of excitatory neuron
                o.ssAMPA = Var(7*o.p.Ne+1:8*o.p.Ne,:);%fraction of open AMPA receptor channels on dendrite of excitatory neuron
                o.ssNMDA = Var(8*o.p.Ne+1:9*o.p.Ne,:);%fraction of open NMDA receptor channels on dendrite of excitatory neuron
                o.h = Var(9*o.p.Ne+1:10*o.p.Ne,:);%sodium channel for soma of excitatory neuron
                o.hd = Var(10*o.p.Ne+1:11*o.p.Ne,:);%sodium channel for dendrite of excitatory neuron
                o.n = Var(11*o.p.Ne+1:12*o.p.Ne,:);%delayed rectifier channel for soma of excitatory neuron
                o.ha = Var(12*o.p.Ne+1:13*o.p.Ne,:);%fast A-type sodium channel for soma of excitatory neuron
                o.mks = Var(13*o.p.Ne+1:14*o.p.Ne,:);%non-inactivating sodium channel for soma of excitatory neuron
                
                o.Vi = Var(o.Netot+1:o.Netot+o.p.Ni,:); %membrane potential of inhibitoroy neuron mV
                o.siAMPA  = Var(o.Netot+o.p.Ni+1:o.Netot+2*o.p.Ni,:); %fraction of open AMPA receptor channels on inhibitory neuron
                o.xiNMDA = Var(o.Netot+2*o.p.Ni+1:o.Netot+3*o.p.Ni,:); %second messenger for NMDA receptor channels on inhibory neuron
                o.siNMDA = Var(o.Netot+3*o.p.Ni+1:o.Netot+4*o.p.Ni,:); %fraction of open NMDA receptor channels on inhibitory neuron
                o.siGABA = Var(o.Netot+4*o.p.Ni+1:o.Netot+5*o.p.Ni,:); %fraction of open GABA receptor channels on inhibitory neuron
                o.hi = Var(o.Netot+5*o.p.Ni+1:o.Netot+6*o.p.Ni,:); %sodium channel for inhibitory neuron
                o.ni = Var(o.Netot+6*o.p.Ni+1:o.Netot+7*o.p.Ni,:); %sodium channel for inhibitory neuron
            end
            
            function nn = Netot(o)
                nn = 14*o.p.Ne;
            end
            
            function nn = Nitot(o)
                nn = 7*o.p.Ni;
            end
            
            function currentDensity = INa(o)
                %sodium current (soma)
                alpham = 0.1*(o.Vs+33)./(1 - exp(-(o.Vs+33)/10));
                betam = 4 * exp(-(o.Vs +53.7)/12);
                minf = alpham ./ (alpham + betam);
                currentDensity = o.p.gNa .* minf.^3 .* o.h .* (o.Vs - o.p.VNa); %[mS/cm^2]*[mV] = [1e-2 * A/m^2]
            end
            
            function currentDensity = INad(o)
                %sodium current (dendrite)
                alpham = 0.1*(o.Vd+33)./(1 - exp(-(o.Vd+33)/10));
                betam = 4 * exp(-(o.Vd +53.7)/12);
                minf = alpham ./ (alpham + betam);
                currentDensity = o.p.gNa .* minf.^3 .* o.hd .* (o.Vd - o.p.VNa); %[mS/cm^2]*[mV] = [1e-2 * A/m^2]
            end
            
            function currentDensity = INai(o)
                %sodium current (inhibitory neuron)
                alphami = 0.5*(o.Vi+35)./(1 - exp(-(o.Vi+35)/10));
                betami = 20 * exp(-(o.Vi + 60)/18);
                minfi = alphami ./ (alphami + betami);
                currentDensity = o.p.gNai .* minfi.^3 .* o.hi .* (o.Vi - o.p.VNai); %[1e-2 * A/m^2]
            end
            
            function currentDensity = IK(o)
                %potassium current (soma)
                currentDensity = o.p.gK .* o.n.^4 .* (o.Vs - o.p.VK); %[1e-2 * A/m^2]
            end
            
            function currentDensity = IKi(o)
                %potassium current (inhibitory neuron)
                  currentDensity = o.p.gKi .* o.ni.^4 .* (o.Vi - o.p.VKi); %[1e-2 * A/m^2]
            end
        
            function currentDensity = ILs(o)
                %leak current (soma)
                currentDensity = o.p.gLs .* (o.Vs - o.p.VL); %[1e-2 * A/m^2]
            end
            
            function currentDensity = ILd(o)
                %leak current (dendrite)
                currentDensity = o.p.gLd .* (o.Vd - o.p.VL); %[1e-2 * A/m^2]
            end
            function currentDensity = ILi(o)
                %leak current (inhibitory neuron)
                currentDensity = o.p.gLi .* (o.Vi - o.p.VLi); %[1e-2 * A/m^2]
            end
            
            function currentDensity = IA(o)
                %fast A-type potassium current (soma)
                mainf = 1./(1+exp(-(o.Vs+50)/20));
                
                %should use main
                %should use ha NOT hainf
                currentDensity = o.p.gA .* mainf.^3 .* o.ha .* (o.Vs - o.p.VK);
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
                mcainf = 1./(1 + exp(-(o.Vd + 20)/9));
                currentDensity = o.p.gCa .* mcainf.^2 .* (o.Vd - o.p.VCa); %[1e-2 * A/m^2]
                conductance = 1e4 * o.p.Ad * o.p.gCa .* mcainf.^2;
            end
            
            function [currentDensity, conductance] = IKCa(o)
                % calcium dependent potassium current (dendrite)
                currentDensity = o.p.gKCa .* o.Ca./(o.Ca + o.p.KD) .* (o.Vd - o.p.VK); %[1e-2 * A/m^2]
                conductance = 1e4 * o.p.Ad * o.p.gKCa .* o.Ca./(o.Ca + o.p.KD);
            end
            
            function [currentDensity, conductance] = INaP(o)
                % non-inactivating (persistent) sodium current  (dendrite)
                mnapinf = 1./(1+exp(-(o.Vd+55.7)/7.7));
                currentDensity = o.p.gNaP .* mnapinf.^3 .* (o.Vd - o.p.VNa); %[1e-2 * A/m^2]
                conductance = 1e4 * o.p.Ad * o.p.gNaP .* mnapinf.^3; %[nS]
            end
            
            function [currentDensity, conductance] = IAR(o)
                % inward rectifier (activated by hyperpolarization) (dendrite)
                % non-inactivating current
                harinf = 1./(1+exp((o.Vd+75)/4));
                currentDensity = o.p.gAR.*harinf.*(o.Vd - o.p.VK); %[1e-2 * A/m^2]
                conductance = 1e4 * o.p.Ad * o.p.gAR.*harinf; %[nS]
            end
            
            
            %% synaptic currents
            function current = Isyns(o) %[1e-8 A]
                %synaptic current to excitatory soma 
                current = 1e-4 * (o.p.gIE .* o.ssGABA .* (o.Vs - o.p.VsynGABA) ...
                    + o.p.gEEsAMPA .* o.ssAMPA .* (o.Vs - o.p.VsynAMPA) ...
                    + o.p.gEEsNMDA .* o.ssNMDA .* (o.Vs - o.p.VsynNMDA));
            end
            
            function current = Isynd(o) %[1e-8 A]
                %synpatic current to excitatory dendrite 
                current = 1e-4 * (o.p.gEEdAMPA .* o.sdAMPA .* (o.Vd - o.p.VsynAMPA) ...
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
                dif = o.p.alphaGABA .* (o.p.WIE * o.ff(o.Vi)) - o.ssGABA./o.p.tauGABA;
            end
            
            function dif = dsdAMPA(o)
%                  dif= o.p.alphaAMPA .* (o.p.WEE * o.ff(o.Vs)) .* (1-o.sdAMPA) - o.sdAMPA./o.p.tauAMPA;% Tegner 2002
                 dif= o.p.alphaAMPA .* (o.p.WEEd * o.ff(o.Vs)) - o.sdAMPA./o.p.tauAMPA;%Pinsky 1994; Compte 2003
            end
                       
            function dif = dsdNMDA(o)
                %                  dif = o.p.alphasNMDA .* o.xdNMDA .* (1 - o.sdNMDA) - o.sdNMDA./o.p.tausNMDA; %Tegner 2002; added xiNMDA to Compte 2003
                 %                dif = o.p.alphasNMDA .* (1 - o.sdNMDA) - o.sdNMDA./o.p.tausNMDA; %written in Compte 2003
                 dif = o.p.alphasNMDA .* (1 - o.sdNMDA) .* (o.p.WEEd * o.ff(o.Vs)) - o.sdNMDA./o.p.tausNMDA; %added network input to Compte 2003
            end
            
            function dif = dssAMPA(o)
%                  dif= o.p.alphaAMPA .* (o.p.WEE * o.ff(o.Vs)) .* (1-o.sdAMPA) - o.sdAMPA./o.p.tauAMPA;% Tegner 2002
                 dif= o.p.alphaAMPA .* (o.p.WEEs * o.ff(o.Vs)) - o.ssAMPA./o.p.tauAMPA;%Pinsky 1994; Compte 2003
            end
                       
            function dif = dssNMDA(o)
                %                  dif = o.p.alphasNMDA .* o.xdNMDA .* (1 - o.sdNMDA) - o.sdNMDA./o.p.tausNMDA; %Tegner 2002; added xiNMDA to Compte 2003
                 %                dif = o.p.alphasNMDA .* (1 - o.sdNMDA) - o.sdNMDA./o.p.tausNMDA; %written in Compte 2003
                 dif = o.p.alphasNMDA .* (1 - o.ssNMDA) .* (o.p.WEEs * o.ff(o.Vs)) - o.ssNMDA./o.p.tausNMDA; %added network input to Compte 2003
            end
            
            function dif = dsiAMPA(o)
%                 dif = o.p.alphaAMPA .* (o.p.WEI * o.ff(o.Vs)) .* (1 - o.siAMPA) - o.siAMPA./o.p.tauAMPA; %Tegner 2002
                dif = o.p.alphaAMPA .* (o.p.WEI * o.ff(o.Vs)) - o.siAMPA./o.p.tauAMPA; %Pinsky 1994; Compte 2003
            end
            
            function dif = dxiNMDA(o)
                dif = o.p.alphaxNMDA .* (o.p.WEI * o.ff(o.Vs)) - o.xiNMDA./o.p.tauxNMDA;%wang 1999 eq(4); Compte 2003
            end
            
            function dif = dsiNMDA(o)
%                  dif = o.p.alphasNMDA .* o.xiNMDA .* (1 - o.siNMDA) - o.siNMDA./o.p.tausNMDA; %Tegner 2002; added xiNMDA to Compte 2003
%                dif = o.p.alphasNMDA .*  (1 - o.siNMDA) - o.siNMDA./o.p.tausNMDA; %written in Compte 2003 > produce constant siNMDA
                dif = o.p.alphasNMDA .*  (1 - o.siNMDA) .* (o.p.WEI * o.ff(o.Vs)) - o.siNMDA./o.p.tausNMDA; %written in Compte 2003
            end
            
            function dif = dsiGABA(o)
%                 dif = o.p.alphaGABA .* (o.p.WII * o.ff(o.Vi)) .* (1 - o.siGABA) - o.siGABA./o.p.tauGABA;
                dif = o.p.alphaGABA .* (o.p.WII * o.ff(o.Vi))  - o.siGABA./o.p.tauGABA;
            end
            
            %% channel activation variables
            function dif = dh(o)
                %sodium channel excitatory soma
                alphah = 0.07 * exp(-(o.Vs + 50)/10);
                betah = 1./(1 + exp(-(o.Vs + 20)/10));
                dif = o.p.phih.*(alphah.*(1-o.h) - betah.*o.h);
            end
            
             function dif = dhd(o) %ADDED
                %sodium channel excitatory dendrite
                alphah = 0.07 * exp(-(o.Vd + 50)/10);
                betah = 1./(1 + exp(-(o.Vd + 20)/10));
                dif = o.p.phih.*(alphah.*(1-o.hd) - betah.*o.hd);
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
           
             %% current balance
             function dif = dVs(o, extCurrent) %[mV]
                 if nargin<2
                     extCurrent = 0; %[nA]
                 end
                 %soma of excitatory neuron
                 dif = (-o.p.As.*(o.ILs+o.INa+o.IK+o.IA+o.IKS+o.IKNa) ...
                     - o.Isyns - 0.1*o.p.gsd.*(o.Vs-o.Vd) + 0.1*extCurrent)./o.p.Cm./o.p.As; %+ 0.025
             end
             
             function dif = dVd(o, extCurrent) %[mV]
                  if nargin<2
                     extCurrent = 0; %[nA]
                 end
                 %dendrite of excitatory neuron
                 dif = (-o.p.Ad.*(o.ILd+o.ICa+o.IKCa+o.INad+o.IAR) ...
                     - o.Isynd - 0.1*o.p.gsd.*(o.Vd-o.Vs) + 0.1*extCurrent)./o.p.Cm./o.p.Ad;
                 %replaced INaP w INad
             end
             %for default gL:
             %ILd=8
             %ICa = -25
             %IKCa = 16
             %>Vd stays at 55mV
                
             function dif = dVi(o, extCurrent) %[mV]
                 %inhibitory neuron
                 if nargin<2
                     extCurrent = 0; %[nA]
                 end
                 dif = (-o.p.Ai.*(o.ILi+o.INai+o.IKi) - o.Isyni + 0.1*extCurrent)./o.p.Cm./o.p.Ai;
             end
        
             function dif = dCa(o)
                 %intracellular calcium concentration (dendrite)
                 dif = -1e-2 * o.p.alphaca .* o.p.Ad * o.ICa - o.Ca./o.p.tauca; %[mM/ms]
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
                allDif(5*o.p.Ne+1:6*o.p.Ne,1) = o.dsdAMPA;%fraction of open AMPA receptor channels on dendrite of excitatory neuron
                allDif(6*o.p.Ne+1:7*o.p.Ne,1) = o.dsdNMDA;%fraction of open NMDA receptor channels on dendrite of excitatory neuron
                allDif(7*o.p.Ne+1:8*o.p.Ne,:) = o.dssAMPA;%fraction of open AMPA receptor channels on dendrite of excitatory neuron
                allDif(8*o.p.Ne+1:9*o.p.Ne,:) = o.dssNMDA;%fraction of open NMDA receptor channels on dendrite of excitatory neuron
                allDif(9*o.p.Ne+1:10*o.p.Ne,:) = o.dh;%sodium channel for soma of excitatory neuron
                allDif(10*o.p.Ne+1:11*o.p.Ne,:) = o.dhd;%sodium channel for dendrite of excitatory neuron
                allDif(11*o.p.Ne+1:12*o.p.Ne,:) = o.dn;%delayed rectifier channel for soma of excitatory neuron
                allDif(12*o.p.Ne+1:13*o.p.Ne,:) = o.dha;%fast A-type sodium channel for soma of excitatory neuron
                allDif(13*o.p.Ne+1:14*o.p.Ne,:) = o.dmks;%non-inactivating sodium channel for soma of excitatory neuron
                
                allDif(o.Netot+1:o.Netot+o.p.Ni,1) = o.dVi; %membrane potential of inhibitoroy neuron
                allDif(o.Netot+o.p.Ni+1:o.Netot+2*o.p.Ni,1) = o.dsiAMPA; %fraction of open AMPA receptor channels on inhibitory neuron
                allDif(o.Netot+2*o.p.Ni+1:o.Netot+3*o.p.Ni,1) = o.dxiNMDA; %second messenger for NMDA receptor channels on inhibory neuron
                allDif(o.Netot+3*o.p.Ni+1:o.Netot+4*o.p.Ni,1) = o.dsiNMDA; %fraction of open NMDA receptor channels on inhibitory neuron
                allDif(o.Netot+4*o.p.Ni+1:o.Netot+5*o.p.Ni,1) = o.dsiGABA; %fraction of open GABA receptor channels on inhibitory neuron
                allDif(o.Netot+5*o.p.Ni+1:o.Netot+6*o.p.Ni,1) = o.dhi; %sodium channel for inhibitory neuron
                allDif(o.Netot+6*o.p.Ni+1:o.Netot+7*o.p.Ni,1) = o.dni; %potassium channel for inhibitory neuron
            end
            
            function conductance = condSyn_exc_d(o) %correct?
                %total excitatory synaptic conductance of dendrite[nS]
                conductance = o.p.gEEdAMPA .* o.sdAMPA + o.p.gEEdNMDA .* o.sdNMDA;
            end
            
            function conductance = condSyn_exc_s(o) %correct?
                %total excitatory synaptic conductance of dendrite[nS]
                conductance = o.p.gEEsAMPA .* o.ssAMPA + o.p.gEEsNMDA .* o.ssNMDA;
            end
            
            function conductance = condSyn_inh(o) %correct?
                %total inhibitory synaptic conductance of soma[nS]
                conductance = o.p.gIE .* o.ssGABA;
            end
        end %methods
        
        methods(Static)
            function spike = ff(Vpre)
                %works as a spike detector
                spike = 1./(1 + exp(-(Vpre - 20)/2));
            end          
        end
end %classdef