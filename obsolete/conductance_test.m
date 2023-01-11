cd('Z:\Shared\Daisuke\sandbox\Compte2003');
run('param_ds_single.m');

q10 = 2.3;
temp = 23;
celsius = 37;
tadj = q10^((celsius - temp)/10);

%% soma
o.Vs = -100:150;
%h
alphah = 0.07 * exp(-(o.Vs + 50)/10);
betah = 1./(1 + exp(-(o.Vs + 20)/10));
yyaxis left
plot(o.Vs, alphah./(alphah+betah));
hold on
yyaxis right
plot(o.Vs, 1./(alphah+betah));
legend('inf','tau');

%A current
hainf = 1./(1+exp((o.Vs+80)/6));
yyaxis left
plot(o.Vs,hainf);

%Na activation
alpham = 0.1*(o.Vs+33)./(1 - exp(-(o.Vs+33)/10));
betam = 4 * exp(-(o.Vs +53.7)/12);
plot(o.Vs, alpham./(alpham+betam));
hold on
yyaxis right
plot(o.Vs, 1./(alpham+betam));
legend('inf','tau');


%% dendrite
o.Vd = -100:150;

%sodium inactivation @dendrite Mainen 1996
alphah = 0.024 * (o.Vd + 40) ./ ( 1 - exp(-(o.Vd + 40)/5));
betah = -0.0091 * (o.Vd + 65) ./ ( 1 - exp((o.Vd + 65)/5));
tauh = 1./(alphah+betah);
hinf = 1./(1+exp(o.Vd+55)/6.2);
yyaxis left
plot(o.Vd,hinf);
hold on
yyaxis right
plot(o.Vd,tauh);

%sodium activation @dendrite Mainen 1996
alpham = 0.1*(o.Vd+25)./(1 - exp(-(o.Vd+25)/9));
betam = -0.124 * (o.Vd+25) ./ (1 - exp((o.Vd +25)/9));
plot(o.Vs, alpham./(alpham+betam));
hold on
yyaxis right
plot(o.Vs, 1./(alpham+betam));

%calcium activation @dendrite Mainen 1996 ???
alphamca = 0.055*(o.Vd + 27)./(1 - exp(-(27+o.Vd)/3.8));
betamca = 0.94*exp(-(o.Vd + 75)/17);
yyaxis left
plot(o.Vd, alphamca./(alphamca+betamca));
hold on
yyaxis right
plot(o.Vd, 1./(alphamca+betamca));
             
%calcium inactivation @dendrite. Mainen 1996 ???
%inactivation happens at lower potential than activation
alphahca = 4.57 * 1e-4 .* exp(-(o.Vd+13)/50);
betahca = 0.0065 ./ (1 + exp(-(o.Vd + 15)/28));
yyaxis left
plot(o.Vd, alphahca./(alphahca+betahca));
hold on
yyaxis right
plot(o.Vd, 1./(alphahca+betahca));

%slow non-inactivating K Mainen96
 alphaq = 1e-4*(o.Vd+30)./(1 - exp(-(o.Vd+30)/9));
 betaq = -1.1e-4*(o.Vd+30)./(1 - exp((o.Vd+30)/9));%1.1e-4?
 yyaxis left
 plot(o.Vd, alphaq./(alphaq+betaq));
hold on
yyaxis right
plot(o.Vd, 1./(alphaq+betaq));

 %slow non-inactivating K
 T=34;
 tau = 1000./(3.3*(exp((o.Vd+35)/40)+exp(exp(-(o.Vd+35)/20)))) ./ 3.^((T-22)/10); 
 inf = 1./(1+exp(-(o.Vd+35)/10));
 yyaxis left
 plot(o.Vd, inf);
hold on
yyaxis right
plot(o.Vd, tau);

%KCa activation (xi) Mainen96
o.Ca = p.KD; %mM
alphaxi = 0.01*1e-3*o.Ca; %Mainen assumes uM whereas this class assumes mM
betaxi = 0.02;
tau = 1/(alphaxi+betaxi);

%KCa opening probability Compte 2003
compteKca = o.Ca./(o.Ca + p.KD);


          