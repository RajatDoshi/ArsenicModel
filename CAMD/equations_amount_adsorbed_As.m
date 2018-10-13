function[N1,Nb,Xb,Nstar,Xstar,SPC,Gamma_bulk,Gamma_star,mtotal,V,T,R,x1s,phistar]=equations_amount_adsorbed_As(m0,n1b,K2)

%---- Initial inforamtion (produce water after pretreatment) 
V=1000;         % Volume of solution <ml>  100
T=300;       % Temperature <K> assuming ambient temperature
R=8.314;        %  ideal gas constant <J/ mol K>

% from the report%
F1=[5 6.3 85.4]; %flowrate <mg/L> 
molecule=K2;

% -- Initial phase ----
% properties per group contribution 

[MW0,Rk,Qk,Vp,SPC,N1,x1,MWra]=input_data(molecule,F1,V);
MW0;
%note: moles fraction are taken from excel file: mass balance adsorption process 
n11=N1(1);     % Moles of As in the initial feed
n21=N1(2);     % Moles of Na in the initial feed
n31=N1(3);     % Moles of HCO3 in the initial feed
n41=N1(4);     % Moles of H2O in the initial feed
n1=sum(N1);    % Total moles in the feed 

% id of groups and number of groups based on molecule structure
[Nid,Nindex]=molecule_structure(molecule);

% -- Bulk phase -----
n2b=n21;            % Moles of H20 in the bulk phase
n3b=n31;            % Moles of Na in the bulk phase
n4b=n41;            % Moles of HCO3 in the bulk phase
Nb=[n1b n2b n3b n4b];
nb=sum(Nb);         % Total moles in the bulk phase
Xb=Nb./nb; % moles fraction in the bulk phase

% -- Adsorbate-solid-solution ASS phase ---- 
n0star=m0/MW0;      % Moles of solid in ASS phase
n1star=n11-n1b;     % Moles of As in ASS phase
Nstar=[n1star n0star]; %[As Adsorbent ]
nstar=sum(Nstar);   % Total moles in the ASS phase
Xstar=Nstar./nstar; % moles fraction in the ASS phase

% -- surface phase ---- 
x1s=Xstar(1)/(1-Xstar(2));   
ns=nstar*(1-Xstar(2));
Xs=[x1s];

% ---x01star and x02star represent the activity coef, of binary ASS, that is
% ---the pure adsorbed liquid and solid: therefore normalize: 
x01star=Xstar(1)/(Xstar(1)+Xstar(2));
n01star=x01star*(n1star+n0star);
% Natural Log of pure component adsorption activity coefficient
LNgamma0_star=log(1+(1/(SPC(1)*MW0)));   % As

%--- Activities coefficient bulk phase  -------
[LNgamma_bulk,Gamma_bulk]=UNIFAC_bulk_phase(Xb,T);

%--- Activities coefficient ASS phase  -------
[LNgammaGE_star,Gamma_star]=UNIFAC_ASS(T,LNgamma0_star,Xstar,Rk,Qk,Nid,Nindex);

[LNgamma_s,Gamma_s]=UNIFAC_SURFACE(Xs,T);

GEstar0(1)=n01star*LNgamma0_star(1);

GEstar=nstar*(Xstar(1)*LNgammaGE_star(1)+Xstar(2)*LNgammaGE_star(2));
GES=R*T*ns*(Xs(1)*LNgamma_s(1));
phistar(1)=(R*T/m0)*(GEstar-GES-GEstar0(1));

mstar(1)=n1star*MWra;

mtotal=m0+mstar(1);