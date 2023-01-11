%% default parameters for compte_ds_mainen.m      

if doSingle
    p.Ne = 1;
    p.Ni = 1;
else
    p.Ne=1024;%: #excitatory neurons
    p.Ni=256;%: #inhibitory neurons
end

p.Cm = 1;%: membrane capacitance 1 uF/cm^2
p.Cmd = 0.75;%:Mainen 1996
p.As = 0.015;%: soma of excitatory neuron surface area  mm^2
p.Ad = 100*1e-6 * 165;%0.035;%: dendrite of excitatory neuron surface area  mm^2
%< from https://senselab.med.yale.edu/modeldb/ShowModel?model=2488&file=/cells/demofig2.hoc#tabs-2
p.Ai = 0.02;%: inhibitory neuron surface area mm^2

%soma-dendrite geometry
p.gsd = 1.75;%: coupling between soma-dendrite +-0.1 uS
%rho = p.Ad/p.As ;%[0-1] area(dendrite) / area(soma). Pinsky94; Mainen96

%temperature factor
p.phih = 4;
p.phin = 4;
p.phihi = 1;
p.phini = 1;
p.phiha = 1; %Golomb Amitai 97 A.14
p.phimks = 1;
p.phid = 3.2094; %for all dendrite-related components.
%computed as q10^((celsius - temp)/10), where q10=2.3, celsius=37, temp=23


p.gNa = 50; %: conductance  mS/cm^2
p.gNai = 35;%: inhibitory neuron mS/cm^2
p.gK = 10.5;% mS/cm^2%
p.gKi = 9;%: inhibitory neuron  mS/cm^2
p.gLs = 0.0667+0.01;%067*randn(p.Ne,1);%mS/cm^2
p.gLi = 0.1025+0.0025*randn(p.Ni,1);%mS/cm^2;
p.gA = 1;% mS/cm^2;
p.gKS = 0.576;% mS/cm^2
p.gKNa = 1.33; %mS/cm^2
fac = 0.1;%
p.gLd = fac/3e4;%Mainen computed from membrane resistivity (ohm-cm^2)
p.gNad = fac * 15; %dendrite mS/cm^2 Mainen 96
p.gCa = fac * 0.3;% mS/cm^2 Mainen 96
p.gKCa = fac * 3;%0.57;% mS/cm^2 Mainen 96
p.gKm = fac * 0.1;% Mainen 96

%synaptic condactance
p.gIEs = 4.15;% nS %inhibitory to excitatory soma
p.gIEd = 4.15;% nS %inhibitory to excitatory dendrite (added)
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

%channel time constant
p.tauha = 15;% ms

%intracellular Ca2+ concentration constant
p.alphaca = 0.005;% uM/(nA*ms)
p.tauca = 200;%ms. Mainen 1996.  150;% Compte 2003
%p.Cainf = 0.1*1e-3; %mM Mainen 1996
p.Cainf = 100e-6; %mM https://senselab.med.yale.edu/ModelDB/showmodel?model=2488&file=/cells/cad.mod#tabs-2
%p.KD = 1e-3*30; %mM

%intracellular Na+ concentration constant
p.alphana = 0.01; %mA/(nA*ms)
p.Rpump = 0.018;% mM/ms
p.Naeq = 9.5; %mM

%soma reversal potential
p.VNa = 55;% mV
p.VK = -100;% mV
p.VL = -60.95+0.3*randn(p.Ne,1);

p.Vdif = 25;%without external input, Vs=-65, Vd = -90

%dendrite reversal potential
p.VCa = 140+p.Vdif;%neuroDB. in Pinsky 1994 %120 Compte 2003
p.VLd = -70+p.Vdif; %mainen96 TOBE FIXED
p.VKd = -90+p.Vdif; %mainen96 TOBE FIXED
p.VNad = 50+p.Vdif; %mainen96 TOBE FIXED


%inhibitory neuron reversal potential
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

cFac = 5;
p.Cees = cFac*20; %# excitatory to excitatory soma connection. 
p.Ceed = cFac*20; %# excitatory to excitatory dendrite connection. 
p.Cies = cFac*20; %# inhibitory to excitatory soma connection. 
p.Cied = cFac*20; %# inhibitory to excitatory dendrite connection. (added)
p.Cei = cFac*5; %# excitatory to inhibitory connection. 
p.Cii = cFac*5; %# inhibitory to inhibitory connection. 

if doSingle
    p.WIEs = 0;
    p.WIEd = 0; %(added)
    p.WEEs = 0;
    p.WEEd = 0;
    p.WEI = 0;
else
    p.WIEs = cell2mat(arrayfun(@(x)getConnectivity(p.pos_e, x, p.sigma_i, p.Cies), p.pos_i, 'UniformOutput', false));
    p.WIEd = cell2mat(arrayfun(@(x)getConnectivity(p.pos_e, x, p.sigma_i, p.Cied), p.pos_i, 'UniformOutput', false));%(added)
    p.WEEs = cell2mat(arrayfun(@(x)getConnectivity(p.pos_e, x, p.sigma_e, p.Cees), p.pos_e, 'UniformOutput', false));
    p.WEEd = cell2mat(arrayfun(@(x)getConnectivity(p.pos_e, x, p.sigma_e, p.Ceed), p.pos_e, 'UniformOutput', false));
    p.WEI = cell2mat(arrayfun(@(x)getConnectivity(p.pos_i, x, p.sigma_e, p.Cei), p.pos_e, 'UniformOutput', false));
    p.WII = cell2mat(arrayfun(@(x)getConnectivity(p.pos_i, x, p.sigma_i, p.Cii), p.pos_i, 'UniformOutput', false));
end

%% external current
p.VsExtCurrent = zeros(p.Ne,1);
p.VdExtCurrent = zeros(p.Ne,1);
p.ViExtCurrent = zeros(p.Ni,1);

