function [MW0,Rk,Qk,Vp,SPC,n1,x1,MWra]=input_data(molecule,F1,V)
% -----------Notes--------------------------
% MW0: molecular weight of adsorbent
% Rk and Qk: volume and surfase van der waals
% Vp: pore volume cm3/g
% SPC: Surfase phase capacity <mol/g>
% F1: Initial feed <mg/l>
% V: volume of solution <ml>
% n1: Initial moles <mol>

%----initial moles 
%properties     MW(g/mol)   density(g/ml)
proper_initial=[74.92       5.78;   %As
                22.989      0.968;  %Na
                61.0168     2.2;    %HCO3
                18.02       1;];    %H20

MWra=proper_initial(1,1);
MW=zeros(1,length(molecule)); 
vp=zeros(1,length(molecule)); 
ini=length(F1)+1;     
vi=zeros(1,ini);
for jj=1:ini-1
    mass(jj)=F1(jj)*(1/1000)*(1/1000)*V;  % mass of each component <g>
    vi(jj)= mass(jj)*(1/proper_initial(jj,2)); % volume of each component <ml>
end
% mass and volume of water
vi(ini)=V-sum(vi);
mass(ini)=vi(ini)*proper_initial(ini,2);

for jj=1:ini
    n1(jj)=mass(jj)/proper_initial(jj,1);   % moles of 
end
 x1=n1/sum(n1);   % mole fraction of the feed 

%-----For Molecule------
%properties MW    Vpi
proper=[60.083    0.03;     %SiO2
        101.961   0.181;    %Al2O3
        61.979    0.004;    %Na2O
        40.304    0.848;    %MgO 
        159.687   0.08;     %Fe2O3
        94.187    0.483;    %K2O
        56.077    0.093;    %CaO
        22.990    0];       %Na

    
   %    Rk     Qk  for all group plus Arsenic  
RkQk=[ 3.610 5.380;    %As
       1.389 1.267;    %SiO2
       1.109 1.103;    %Al2O3
       1.688 1.103;    %Na2O
       1.013 1.041;    %MgO 
       1.221 1.167;    %Fe2O3
       2.716 1.844;    %K2O
       1.752 1.438;   %CaO
       3.000 3.000];   %Na

for ii=1:length(molecule)
    MW(ii)= proper(ii,1)*molecule(ii);
    vp(ii)= proper(ii,2)*molecule(ii);
end

MW0= sum(MW);
Vp = 1/(sum(vp)+3.356); 
Rk=RkQk(:,1);
Qk=RkQk(:,2);

MV(1)=12.96;    % Molar volume of Arsenic <cm3/mol>
SPC(1)=Vp/MV(1); % Surfase phase capacity Arsenic <mol/g>
end