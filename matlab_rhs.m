NNa1=kmrgd(1);
NK1=kmrgd(2);
NCl1=kmrgd(3);
Wi1=kmrgd(4);
NNa2=kmrgd(5);
NK2=kmrgd(6);
NCl2=kmrgd(7);
Wi2=kmrgd(8);
O2e12=kmrgd(9);
r1=kmrgd(10);
r2=kmrgd(11);
NNa3=kmrgd(12);
NK3=kmrgd(13);
NCl3=kmrgd(14);
Wi3=kmrgd(15);
NNa4=kmrgd(16);
NK4=kmrgd(17);
NCl4=kmrgd(18);
Wi4=kmrgd(19);
O2e34=kmrgd(20);
r3=kmrgd(21);
r4=kmrgd(22);


C                 = 20;                          % [pF], membrane capacitance
F                 = 96485.333;                    % [C/mol], Faraday's constant
R                 = 8.3144598*1000;               % [(CmV)/(mol K)], Universal gas constant
T                 = 310;                          % [K], Absolute temperature
rcell             = 7.815926417967719;            % [mu m], radius of spherical cell
Ari0      =  4*pi*(rcell)^2/1000;             % [mu m^2], Membrane surface of neuroning                       % [1000 mu m^2]
Wi0      = (4/3)*pi*(rcell)^3/1000;         % [mu m^3], Volume of spherical cell   ng                       % [1000 mu m^3]
NAi  = 300.01357719312625;
% Part 1B. ATP-related parameters and ATP-pump
PumpStrength           = 60;                                 % [pA], Pump current scaling (including strenght of pump) / NaK pump rate

% Part 1C. Transmembrane receptors
% Part 1E. KCL cotransporter
UKCl      = 1e-6; % [fmol/ms/mV] the cotransporter strength
PWi         = 2*1e-14;   %  Water permeability
% Part 1e. Ion related parameters
NaCe0      = 152;   % [mMol/l = mM], Initial concentration of Sodium (o = extracellular, e=excitatory)
NaCi0    = 13;    % [mMol/l = mM], ICS sodium concentration of excitatory population
KCe0       = 3;    % [mMol/l = mM], ECS potassium concentration
KCi0     = 145 ;  % [mMol/l = mM], Initial ICS potassium concentration of excitatory population
ClCe0      = 135;
ClCi0    = 8.0;   % [mMol/l = mM], ICS chloride concentration of excitatory population
Vi0  = -65.5;
% Part 1G. Leak permeabilities (with voltage-gated currents)
PNaL         = 1.2850711297046679e-6;% [1000 mum^3/ms], Leak Na+ permeability (for -65 mV)
PKL          = 1.2527625265033467e-5;% [1000 mum^3/ms], Leak K+ permeability (for -65mV)
PClL         = 2.812973433749513e-6;% [1000 mum^3/ms], Leak Cl- permeability (for -65mV)
% Part 1H. Permeabilities
PNaG        = 80*1e-5             ;                  % [1000 mum^3/ms] Maximal transient Na+ permeability
PKG         = 40*1e-5              ;                 % [1000 mum^3/ms] Maximal delayed rectifier K+ permeability
PClG        = 1.95*1e-5             ;                % [1000 mum^3/ms] Maximal voltage-gated Cl- permeability
PNaSyn      = 1;
PClSyn      = 1;
% Part 1I. Firing Rate related parameters
sigma       = 3;
% Part 1Ia. Coefficients of polynomial Ithreshold
O2e_fac  = 0.05; % 0.1875
syn_fac  = 1;
O2bath  = 2;
O2_baseline  = 1.75;
O2_alpha  = 1/6;
O2_lambda  = 1;
PvATP  = 0.24091128665044184;
O2_diff  = 0.00011747326098812187;
perc    = 0.5 ;
beta1   = 4  ;
beta2   = 4;
tfinal  = 50.0;
ratio   = 0.8;
saveat  = 0.01;
bandpass  = [5.0; 40.0];
O2e_th_NKA  = par_NKA_th;
O2e_th_vATP  = par_vATP_th;
syn_ion      = par_gsyn;
min_vATP  = 0.1;
I_Ext = 20;
CNa  = 2483.9999999999995;
CK  = 627.9999999999998;
CCl  = 2191.9999999999995;
We0  = 15.999999999999998;
Wtot  = 19.999999999999996;
NAe  = 416.1086175450109;


% First solve synapses
NaCi1 = NNa1/Wi1;
KCi1 = NK1/Wi1;
ClCi1 = NCl1/Wi1       ;
NaCi2 = NNa2/Wi2;
KCi2 = NK2/Wi2;
ClCi2 = NCl2/Wi2;
We12 = Wtot - Wi1 - Wi2;
NaCe12 = (CNa - NNa1 - NNa2)/We12;
KCe12 = (CK - NK1 - NK2)/We12;
ClCe12 = (CCl - NCl1 - NCl2)/We12;
V1 = (F/(C))*(NNa1+NK1-NCl1-NAi);
V2 = (F/(C))*(NNa2+NK2-NCl2-NAi);

NaCi3 = NNa3/Wi3;
KCi3 = NK3/Wi3;
ClCi3 = NCl3/Wi3;
NaCi4 = NNa4/Wi4;
KCi4 = NK4/Wi4;
ClCi4 = NCl4/Wi4;
We34 = Wtot - Wi3 - Wi4;
NaCe34 = (CNa - NNa3 - NNa4)/We34;
KCe34 = (CK - NK3 - NK4)/We34;
ClCe34 = (CCl - NCl3 - NCl4)/We34;
V3 = (F/(C))*(NNa3+NK3-NCl3-NAi);
V4 = (F/(C))*(NNa4+NK4-NCl4-NAi);

conn = [par_conn11 par_conn12 par_conn13 par_conn14; par_conn21 par_conn22 par_con23 par_conn24; par_conn31 par_conn32 par_conn33 par_conn34; par_conn41 par_conn42 par_conn43 par_conn44];
rsyn = [r1;r2;r3;r4];
syn_curr = [-(V1-R*T/F*log(NaCe12/NaCi1)) -(V1 - R*T/F*log(ClCi1/ClCe12));
    -(V2-R*T/F*log(NaCe12/NaCi2)) -(V2 - R*T/F*log(ClCi2/ClCe12));
    -(V3-R*T/F*log(NaCe34/NaCi3)) -(V3 - R*T/F*log(ClCi3/ClCe34));
    -(V4-R*T/F*log(NaCe34/NaCi4)) -(V4 - R*T/F*log(ClCi4/ClCe34))];
syn_curr = [syn_curr syn_curr];
syn_curr_full = rsyn' .* (conn .* syn_curr);
syn_curr_full = [sum(syn_curr_full(:,1:2:end),2) sum(syn_curr_full(:,2:2:end),2)];
syn_curr = (conn .* syn_curr)*rsyn;
% Pop1: Thalamic Excitatory
syn_th = 0.2 ;
syn_act = 1.25;
syn_deact = 0.3;
NaCi = NaCi1;
KCi = KCi1;
ClCi = ClCi1;
NaCe = NaCe12;
KCe = KCe12;
ClCe = ClCe12;
We = We12;
O2e = O2e12;
Wi = Wi1;
V = V1;
INaL = PNaL*GHK(1,NaCi,NaCe,V);
IKL = PKL*GHK(1,KCi,KCe,V);
IClL = PClL*GHK(-1,ClCi,ClCe,V);
alpha_n = 0.01 * (V+34.0)/( 1.0 - exp(-0.1 * (V+34.0)) ); %[no units];
beta_n  = 0.125 * exp(-(V+44.0)/80.0);
alpha_m = 0.1 * (V+30.0)/( 1.0 - exp(-0.1 * (V+30.0)) );
beta_m  = 4.0 * exp(-(V+55.0)/18.0);
alpha_h = 0.07 * exp(-(V+44.0)/20.0);
beta_h  = 1.0/( 1.0 + exp(-0.1 * (V+14.0)) );
m_inf = alpha_m/(alpha_m + beta_m);
n_inf = alpha_n/(alpha_n + beta_n);
h_inf = alpha_h/(alpha_h + beta_h);
INaG = m_inf^3*h_inf*PNaG*GHK(1,NaCi,NaCe,V);
IKG = n_inf^4*PKG*GHK(1,KCi,KCe,V);
IClG = 1/(1+exp(-(V+10)/10))*PClG*GHK(-1,ClCi,ClCe,V);
Ipump = 1/(1+exp((O2e_th_NKA-O2e)/O2e_fac))*PumpStrength*NKA(NaCi,KCe,NaCe,V,O2e);
Ipump1 = Ipump;
IvATP1 =  min_vATP + (1-min_vATP)/(1+exp((O2e_th_vATP-O2e12)/O2e_fac));
JKCl = UKCl*KCl(KCi,KCe,ClCi,ClCe);
SCi = NaCi + KCi + ClCi + NAi/Wi;
SCe = NaCe + KCe + ClCe + NAe/We;
ENa = R*T/F*log(NaCe/NaCi);
EK = R*T/F*log(KCe/KCi);
I = I_Ext - Ipump + syn_curr(1) ;
FR = firing_rate(I,EK, ENa);
rsyn_act = syn_act*(FR)/(FR + syn_th);
rsyn_rhs1 = rsyn_act*(1-rsyn(1)) - syn_deact*rsyn(1);
rhs1 =  [-1/(F)*(INaG + 3*Ipump + INaL) + 1/F*syn_ion*PNaSyn*syn_curr_full(1,1);
    -1/(F)*(IKG - 2*Ipump + IKL) - JKCl;
    1/(F)*(IClG + IClL) - JKCl - 1/F*syn_ion*PClSyn*syn_curr_full(1,2);
    PWi*R*T*(SCi-SCe)];
%Pop2 Thalamic Inhibitory
syn_th = 0.5;
syn_act = 0.5;
syn_deact = 0.003;
NaCi = NaCi2;
KCi = KCi2;
ClCi = ClCi2;
NaCe = NaCe12;
KCe = KCe12;
ClCe = ClCe12;
We = We12;
O2e = O2e12;
Wi = Wi2;
V = V2;
INaL = PNaL*GHK(1,NaCi,NaCe,V);
IKL = PKL*GHK(1,KCi,KCe,V);
IClL = PClL*GHK(-1,ClCi,ClCe,V);
alpha_n = 0.01 * (V+34.0)/( 1.0 - exp(-0.1 * (V+34.0)) ); %[no units]
beta_n  = 0.125 * exp(-(V+44.0)/80.0);
alpha_m = 0.1 * (V+30.0)/( 1.0 - exp(-0.1 * (V+30.0)) );
beta_m  = 4.0 * exp(-(V+55.0)/18.0);
alpha_h = 0.07 * exp(-(V+44.0)/20.0);
beta_h  = 1.0/( 1.0 + exp(-0.1 * (V+14.0)) );
m_inf = alpha_m/(alpha_m + beta_m);
n_inf = alpha_n/(alpha_n + beta_n);
h_inf = alpha_h/(alpha_h + beta_h);
INaG = m_inf^3*h_inf*PNaG*GHK(1,NaCi,NaCe,V);
IKG = n_inf^4*PKG*GHK(1,KCi,KCe,V);
IClG = 1/(1+exp(-(V+10)/10))*PClG*GHK(-1,ClCi,ClCe,V);
Ipump = 1/(1+exp((O2e_th_NKA-O2e)/O2e_fac))*PumpStrength*NKA(NaCi,KCe,NaCe,V,O2e);
Ipump2 = Ipump;
IvATP2 = min_vATP + (1-min_vATP)/(1+exp((O2e_th_vATP-O2e12)/O2e_fac));
JKCl = UKCl*KCl(KCi,KCe,ClCi,ClCe);
SCi = NaCi + KCi + ClCi + NAi/Wi;
SCe = NaCe + KCe + ClCe + NAe/We;
ENa = R*T/F*log(NaCe/NaCi);
EK = R*T/F*log(KCe/KCi);
I = - Ipump + syn_curr(2);
FR = firing_rate(I,EK, ENa);
rsyn_act = syn_act*(FR)/(FR + syn_th);
rsyn_rhs2 = rsyn_act*(1-rsyn(2)) - syn_deact*rsyn(2);
rhs2 =  [-1/(F)*(INaG + 3*Ipump + INaL) + 1/F*syn_ion*PNaSyn*syn_curr_full(2,1);
    -1/(F)*(IKG - 2*Ipump + IKL) - JKCl;
    1/(F)*(IClG + IClL) - JKCl - 1/F*syn_ion*PClSyn*syn_curr_full(2,2);
    PWi*R*T*(SCi-SCe)];
O2_rhs_12 = -O2_alpha*O2_lambda*(1/F)*(Ipump1/Wi1 + Ipump2/Wi2 + PvATP*IvATP1/Wi1 + PvATP*IvATP2/Wi2) + O2_diff*(O2bath - O2e12);

%Pop3 Cortical Excitatory
syn_th = 0.2;
syn_act = 12.5;
syn_deact = 3;
NaCi = NaCi3;
KCi = KCi3;
ClCi = ClCi3;
NaCe = NaCe34;
KCe = KCe34;
ClCe = ClCe34;
We = We34;
O2e = O2e34;
Wi = Wi3;
V = V3;
INaL = PNaL*GHK(1,NaCi,NaCe,V);
IKL = PKL*GHK(1,KCi,KCe,V);
IClL = PClL*GHK(-1,ClCi,ClCe,V);
alpha_n = 0.01 * (V+34.0)/( 1.0 - exp(-0.1 * (V+34.0)) ); %[no units];
beta_n  = 0.345 * exp(-(V+44.0)/80.0);
alpha_m = 0.1 * (V+30.0)/( 1.0 - exp(-0.1 * (V+30.0)) );
beta_m  = 4.0 * exp(-(V+55.0)/18.0);
alpha_h = 0.07 * exp(-(V+44.0)/20.0);
beta_h  = 1.0/( 1.0 + exp(-0.1 * (V+14.0)) );
m_inf = alpha_m/(alpha_m + beta_m);
n_inf = alpha_n/(alpha_n + beta_n);
h_inf = alpha_h/(alpha_h + beta_h);
INaG = m_inf^3*h_inf*PNaG*GHK(1,NaCi,NaCe,V);
IKG = n_inf^4*PKG*GHK(1,KCi,KCe,V);
IClG = 1/(1+exp(-(V+10)/10))*PClG*GHK(-1,ClCi,ClCe,V);
Ipump = 1/(1+exp((O2e_th_NKA-O2e)/O2e_fac))*PumpStrength*NKA(NaCi,KCe,NaCe,V,O2e);
Ipump3 = Ipump;
IvATP3 = min_vATP + (1-min_vATP)/(1+exp((O2e_th_vATP-O2e34)/O2e_fac));
JKCl = UKCl*KCl(KCi,KCe,ClCi,ClCe);
SCi = NaCi + KCi + ClCi + NAi/Wi;
SCe = NaCe + KCe + ClCe + NAe/We;
ENa = R*T/F*log(NaCe/NaCi);
EK = R*T/F*log(KCe/KCi);
I = - Ipump + syn_curr(3) ;
FR = firing_rate(I,EK, ENa);
rsyn_act = syn_act*(FR)/(FR + syn_th);
block = min_vATP + (1-min_vATP)/(1+exp((O2e_th_vATP-O2e)/O2e_fac));
rsyn_rhs3 = block*(rsyn_act*(1-rsyn(3)) - syn_deact*rsyn(3));
rhs3 =  [-1/(F)*(INaG + 3*Ipump + INaL) + 1/F*syn_ion*PNaSyn*syn_curr_full(3,1);
    -1/(F)*(IKG - 2*Ipump + IKL) - JKCl;
    1/(F)*(IClG + IClL) - JKCl - 1/F*syn_ion*PClSyn*syn_curr_full(3,2);
    PWi*R*T*(SCi-SCe)];
%Pop4 Cortical Inhibitory
syn_th = 0.5;
syn_act = 5;
syn_deact = 0.03;
NaCi = NaCi4;
KCi = KCi4;
ClCi = ClCi4;
NaCe = NaCe34;
KCe = KCe34;
ClCe = ClCe34;
We = We34;
O2e = O2e34;
Wi = Wi4;
V = V4;
INaL = PNaL*GHK(1,NaCi,NaCe,V);
IKL = PKL*GHK(1,KCi,KCe,V);
IClL = PClL*GHK(-1,ClCi,ClCe,V);
alpha_n = 0.01 * (V+34.0)/( 1.0 - exp(-0.1 * (V+34.0)) ); %[no units];
beta_n  = 0.345 * exp(-(V+44.0)/80.0);
alpha_m = 0.1 * (V+30.0)/( 1.0 - exp(-0.1 * (V+30.0)) );
beta_m  = 4.0 * exp(-(V+55.0)/18.0);
alpha_h = 0.07 * exp(-(V+44.0)/20.0);
beta_h  = 1.0/( 1.0 + exp(-0.1 * (V+14.0)) );
m_inf = alpha_m/(alpha_m + beta_m);
n_inf = alpha_n/(alpha_n + beta_n);
h_inf = alpha_h/(alpha_h + beta_h);
INaG = m_inf^3*h_inf*PNaG*GHK(1,NaCi,NaCe,V);
IKG = n_inf^4*PKG*GHK(1,KCi,KCe,V);
IClG = 1/(1+exp(-(V+10)/10))*PClG*GHK(-1,ClCi,ClCe,V);
Ipump = 1/(1+exp((O2e_th_NKA-O2e)/O2e_fac))*PumpStrength*NKA(NaCi,KCe,NaCe,V,O2e);
Ipump4 = Ipump;
IvATP4 = min_vATP + (1-min_vATP)/(1+exp((O2e_th_vATP-O2e34)/O2e_fac));
JKCl = UKCl*KCl(KCi,KCe,ClCi,ClCe);
SCi = NaCi + KCi + ClCi + NAi/Wi;
SCe = NaCe + KCe + ClCe + NAe/We;
ENa = R*T/F*log(NaCe/NaCi);
EK = R*T/F*log(KCe/KCi);
I = - Ipump + syn_curr(4);
FR = firing_rate(I,EK, ENa);
rsyn_act = syn_act*(FR)/(FR + syn_th);
block = min_vATP + (1-min_vATP)/(1+exp((O2e_th_vATP-O2e)/O2e_fac));
rsyn_rhs4 = block*(rsyn_act*(1-rsyn(4)) - syn_deact*rsyn(4));
rhs4 =   [-1/(F)*(INaG + 3*Ipump + INaL) + 1/F*syn_ion*PNaSyn*syn_curr_full(4,1);
    -1/(F)*(IKG - 2*Ipump + IKL) - JKCl;
    1/(F)*(IClG + IClL) - JKCl - 1/F*syn_ion*PClSyn*syn_curr_full(4,2);
    PWi*R*T*(SCi-SCe)];
O2_rhs_34 = -O2_alpha*O2_lambda*(1/F)*(Ipump3/Wi3 + Ipump4/Wi4 + PvATP*IvATP3/Wi3 + PvATP*IvATP4/Wi4) + O2_diff*(par_ox_bath*O2bath - O2e34);

dydt = [rhs1; rhs2; O2_rhs_12; rsyn_rhs1; rsyn_rhs2; rhs3; rhs4; O2_rhs_34; rsyn_rhs3; rsyn_rhs4];

    function Ipump = NKA(NaCi, KCe, NaCe, V, O2e)
        FF                 = 96485.333;                    % [C/mol], Faraday's constant
        RR                 = 8.3144598*1000;               % [(CmV)/(mol K)], Universal gas constant
        TT                 = 310;                          % [K], Absolute temperature
        sigmapump = 1/7*(exp(NaCe/67.3)-1) ;
        nka_na            = 13;
        nka_k            =0.2 ;
        fpump = 1/(1+0.1245*exp(-0.1*FF/RR/TT*V) +0.0365*sigmapump*exp(-FF/RR/TT*V));
        Ipump = fpump*(NaCi^(1.5)/(NaCi^(1.5)+nka_na^(1.5)))*(KCe/(KCe+nka_k));
    end

    function curr = GHK(z,XCi, XCe, V)
        FF                 = 96485.333;                    % [C/mol], Faraday's constant
        RR                 = 8.3144598*1000;               % [(CmV)/(mol K)], Universal gas constant
        TT                 = 310;                          % [K], Absolute temperature
        curr = (FF^2/RR/TT)*V*((XCi-XCe*exp(-z*FF*V/RR/TT))/(1-exp(-z*FF*V/RR/TT)));
    end

    function kclcurr = KCl(KCi,KCe,ClCi,ClCe)
        FF                 = 96485.333;                    % [C/mol], Faraday's constant
        RR                 = 8.3144598*1000;               % [(CmV)/(mol K)], Universal gas constant
        TT                 = 310;                          % [K], Absolute temperature
        kclcurr = (RR*TT/FF)*log((KCi*ClCi)/(KCe*ClCe));
    end

    function FR = firing_rate(I,Ek,Ena)
        p00               = -213.8 ;
        p10              = -8.404 ;
        p01              = -3.007;
        p20              = -0.1111;
        p11              = -0.07862;
        p02              = -0.001213;
        p30              = -0.0005002;
        p21              = -0.0005221 ;
        p12              = -1.483e-06 ;
        p03              = 3.779e-05 ;
        % Part 1Ib. Coefficients of polynomial 'kappa'
        q00        = 0.287 ;
        q10        = 0.01094;
        q01        = 0.001027;
        q20        = 0.0001692;
        q11        = -0.00000842;
        q02        = -0.00004839;
        q30        = 0.0000008167;
        q21        = -0.0000003156;
        q12        = -0.0000004421;
        q03        = 0.0000002764;
        % Part 1Ic. Coefficients of polynomial Ithreshold2
        r00         = -185.9;
        r10         = -4.877;
        r01         = -3.674;
        r20         = -0.002785;
        r11         = -0.1044;
        r02         = 0.01561;
        r30         = -1.17e-05;
        r21         = 1.897e-05;
        r12         = 0.0006927;
        r03         = 0.0001089;
        Ith     = p00 + p10.*Ek + p01.*Ena + p20.*(Ek).^2 + p11.*Ek.*Ena + p02.*(Ena).^2 + p30.*(Ek).^3 + p21.*(Ek).^2.*Ena + p12.*Ek.*(Ena).^2 + p03.*(Ena).^3;  % 1st threshold
        kappa   = q00 + q10.*Ek + q01.*Ena + q20.*(Ek).^2 + q11.*Ek.*Ena + q02.*(Ena).^2 + q30.*(Ek).^3 + q21.*(Ek).^2.*Ena + q12.*Ek.*(Ena).^2 + q03.*(Ena).^3;
        Ith2    = r00 + r10.*Ek + r01.*Ena + r20.*(Ek).^2 + r11.*Ek.*Ena + r02.*(Ena).^2 + r30.*(Ek).^3 + r21.*(Ek).^2.*Ena + r12.*Ek.*(Ena).^2 + r03.*(Ena).^3;  % 2nd threshold
        
        fun = @(x)((1/(sigma*sqrt(2*pi)))*exp(-(I-x).^2/(2*sigma^2))).*(kappa.*sqrt(max(0,x-Ith))).*(1-(1+sign(x-Ith2))/2);
        % FR = integral(fun,-Inf,Inf,'AbsTol',1e-6,'RelTol',1e-3);
        FR = (kappa.*sqrt(max(0,I-Ith))).*(1-(1+sign(I-Ith2))/2);
        
    end
