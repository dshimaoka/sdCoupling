      
p.Ne=1;%: #excitatory neurons
p.Ni=1;%: #inhibitory neurons

p.Cm = 1;%: membrane capacitance 1 uF/cm^2
p.As = 0.015;%: soma of excitatory neuron surface area  mm^2
p.Ad = 0.035;%: dendrite of excitatory neuron surface area  mm^2
p.Ai = 0.02;%: inhibitory neuron surface area mm^2

p.gsd = 1.75;%: coupling between soma-dendrite +-0.1 uS

p.gNa = 50; %: conductance  mS/cm^2
p.gNai = 35;%: inhibitory neuron mS/cm^2
p.gK = 10.5;% mS/cm^2
p.gKi = 9;%: inhibitory neuron  mS/cm^2
p.gLs = 0.0667;% +- 0.0067 mS/cm^2 gaussian distributed in population
p.gLd = 0.3;% mS/cm^2 ADDED
p.gLi = 0.1025; %+- 0.0025 mS/cm^2;
p.gA = 1;% mS/cm^2;
p.gKS = 0.576;% mS/cm^2
p.gNaP = 0.0686;% mS/cm^2
p.gAR = 0.0257;% mS/cm^2
p.gCa = 0.43;% mS/cm^2
p.gKCa = 0.57;% mS/cm^2
p.gKNa = 1.33; %mS/cm^2

%synaptic condactance
p.gIE = 4.15;% nS %inhibitory to excitatory soma
p.gEEdAMPA = 5.4;% nS %excitatory soma to excitatory dendrite
p.gEEdNMDA = 0.9;% nS
p.gEEsAMPA = 5.4;% nS %excitatory soma to excitatory soma (added)
p.gEEsNMDA = 0.9;% nS added
p.gEIAMPA = 2.25;% nS
p.gEINMDA = 0.5;% nS
p.gII = 0.165;% nS inhibitory to inhibitory (GABAergic?)

%channel opening time constant
p.alphaGABA = 1; 
p.alphaAMPA = 3.48;
p.alphaxNMDA = 3.48;
p.alphasNMDA = 0.5;

%channel closing time constant
p.tauGABA = 10;% ms
p.tauAMPA = 2;% ms
p.tauxNMDA = 2;% ms
p.tausNMDA = 100;% ms

%temperature factor
p.phih = 4;
p.phin = 4;
p.phihi = 1;
p.phini = 1;
p.phiha = 1; %Golomb Amitai 97 A.14
p.phimks = 1;

%channel time constant
p.tauha = 15;% ms

%intracellular Ca2+ concentration constant
p.alphaca = 0.005;% uM/(nA*ms)
p.tauca = 150;% ms
p.KD = 1e-3*30; %mM

%intracellular Na+ concentration constant
p.alphana = 0.01; %mA/(nA*ms)
p.Rpump = 0.018;% mM/ms
p.Naeq = 9.5; %mM

%resting potential
p.VNa = 55;% mV
p.VCa = 120;% mV %140 in Pinsky 1994
p.VK = -100;% mV
p.VL = -60.95+0.3*randn(p.Ne,1);

p.VNai = 55; %mV
p.VKi = -90; %mV
p.VLi = -63.8+0.15*randn(p.Ni,1);

p.VsynAMPA = 0; %mV
p.VsynNMDA = 0; %mV
p.VsynGABA = -70; %mV

% firing voltage threshold (not used in simulation)
p.Vth = 20; %[mV]


%connectivity
p.sigma_e = 1e-3*250; %mm
p.sigma_i = 1e-3*125; %mm;
p.pos_e = linspace(0,5,p.Ne);
p.pos_i = linspace(0,5,p.Ni);

p.Cees = 20; %# excitatory to excitatory soma connection. 
p.Ceed = 20; %# excitatory to excitatory dendrite connection. 
p.Cie = 20;%20; %# inhibitory to excitatory connection. 
p.Cei = 5; %# excitatory to inhibitory connection. 
p.Cii = 5; %# inhibitory to inhibitory connection. 

p.WIE = 0;%cell2mat(arrayfun(@(x)getConnectivity(p.pos_e, x, p.sigma_i, p.Cie), p.pos_i, 'UniformOutput', false));
p.WEEs = 0;%cell2mat(arrayfun(@(x)getConnectivity(p.pos_e, x, p.sigma_e, p.Cees), p.pos_e, 'UniformOutput', false));
p.WEEd = 0;%cell2mat(arrayfun(@(x)getConnectivity(p.pos_e, x, p.sigma_e, p.Ceed), p.pos_e, 'UniformOutput', false));
p.WEI = 0;%cell2mat(arrayfun(@(x)getConnectivity(p.pos_i, x, p.sigma_e, p.Cei), p.pos_e, 'UniformOutput', false));
p.WII = 0;%cell2mat(arrayfun(@(x)getConnectivity(p.pos_i, x, p.sigma_i, p.Cii), p.pos_i, 'UniformOutput', false));
