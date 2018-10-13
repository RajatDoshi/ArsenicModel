%Adsorption of Ra on montmorillonite  paper: L.L Ames et al. 1983 
%last modify: March 24/14
function [LNgamma_bulk,Gamma_bulk]=UNIFAC_bulk_phase_all_adsorbents(Xb,T,at)

%components in the bulk phase: % Ra (1)   NaCl(2)  Water (3)
%  components: group  #  group  #  group  id     #
%   Radium:      Ra   1  
%   NaCl:        Na+  1  Cl-    1   
%   Water:       H20  1

%Input-----------------------
v_comp_b = [1 0 0 0 0 0 0; 0 1 0 0 0 0 0; 0 0 1 0 0 0 0; 0 0 0 1 0 0 0; 0 1 0 0 1 0 0; 0 0 0 0 0 1 0; 0 0 0 0 0 0 1;];

%X=[x1 x2 x3];
%T=293.15;
% Van der Waals volume and surface parameters 
Rk=[3.61 0.92 0.163 0.152 0.631 1.389 2.856];
Qk=[5.28 1.4 .297 .285 0.7370 1.267 2.015];

Nc=length(Xb);  % number of components 
Nb=length(Qk);  % number of groups in the bulk phase

group_b = [1 0 0 0 0 0 0; 0 1 0 0 0 0 0; 0 0 1 0 0 0 0; 0 0 0 1 0 0 0; 0 1 0 0 1 0 0; 0 0 0 0 0 1 0; 0 0 0 0 0 0 1;];

z=10;
r1=2.05;
%-----------------------------
%------- COMBINATORIAL------------------
for i=1:Nc      % per component
    for k=1:Nb  % per gruop 
        rr(i,k)=v_comp_b(i,k)*Rk(k);
        qq(i,k)=v_comp_b(i,k)*Qk(k);
    end
    r(i)=sum(rr(i,:));
    q(i)=sum(qq(i,:));
    l(i)=(z/2)*(r(i)-q(i))-(r(i)-1);
end

%Molecular surface area fraction and Volume fraction:
for i=1:Nc
    thetaCi(i)=q(i)*Xb(i); % surface area 
    phiCi(i)=r(i)*Xb(i);   % volume fraction
end
thetaCT=sum(thetaCi);
thetaC=thetaCi./thetaCT;  %theta(1)=q(1)*X(1) / (q(1)* X(1) + q(2)*X(2));

phiCT=sum(phiCi);
phiC=phiCi./phiCT;        %phi(1)=r(1)*x1 / (r(1)* x1 + r(2)*x2);
 
% combinatorial activity coefficient for component i:

for i=1:Nc
        LL(i)=Xb(i)*l(i);
end
sumL=sum(LL);

for i=1:Nc
    LNgammaC_b(i)=log(phiC(i)/Xb(i))+(z/2)*q(i)*log(thetaC(i)/phiC(i))+l(i)-(phiC(i)/Xb(i))*sumL;
end

%------ RESIDUAL--------------------------------------------
% group_b interaction parameters for the residual coefficient
ab12=at(1); ab13=at(2);  ab14=at(3); ab15=at(4); ab16=at(5); ab17=at(6); ab21=at(7); ab26=at(8); ab31=at(9); ab36=at(10); ab37=at(11); ab41=at(12); ab46=at(13); ab47=at(14); ab51=at(15); ab52=at(16); ab53=at(17); ab54=at(18); ab56=at(19); ab57=at(20); ab61=at(21); ab63=at(22); ab64=at(23); ab65=at(24); ab67=at(25); ab71=at(26); ab73=at(27); ab74=at(28); ab75=at(29); ab76=at(30);
a=[0 ab12 ab13 ab14 ab15 ab16 ab17; ab21 0 -897.2 22.38 -982.5 ab26 17869.84; ab31 -838.2 0 0 4166.3 ab36 ab37; ab41 -165 0 0 6242.2 ab46 ab47; ab51 ab52 ab53 ab54 0 ab56 ab57; ab61 -185624.64 ab63 ab64 ab65 0 ab67; ab71 -1156.9 ab73 ab74 ab75 ab76 0;];

for n=1:Nb
    for m=1:Nb
        Y(n,m)=exp(-a(n,m)/T);
    end
end

% for the pure component i - group_b fraction parameters:
for i=1:Nc
    XT=0;
    for k=1:Nb
        Xi1(i,k)=v_comp_b(i,k)*Xb(i);
        XTi1(i)=XT+Xi1(i,k);
        XT=XTi1(i);
    end
    for kk=1:Nb
        xp(i,kk)=Xi1(i,kk)/XTi1(i);  %x11
    end          
end

% group_b surface area fraction:
for i=1:Nc
    TT=0;
    for k=1:Nb
        thetai1(i,k)=Qk(k)*xp(i,k);
        thetaTi1(i)=TT+thetai1(i,k);
        TT=thetaTi1(i);
    end
    for kk=1:Nb
        theta(i,kk)=thetai1(i,kk)/thetaTi1(i);  %x11
    end          
end

% gamma for pure component i
for i=1:Nc 
    for k=1:Nb 
        pp1=theta(i,1)*Y(k,1)/(theta(i,1)*Y(1,1)+theta(i,2)*Y(2,1)+theta(i,3)*Y(3,1)+theta(i,4)*Y(4,1)+theta(i,5)*Y(5,1)+theta(i,6)*Y(6,1)+theta(i,7)*Y(7,1));
        pp2=theta(i,2)*Y(k,2)/(theta(i,1)*Y(1,2)+theta(i,2)*Y(2,2)+theta(i,3)*Y(3,2)+theta(i,4)*Y(4,2)+theta(i,5)*Y(5,2)+theta(i,6)*Y(6,2)+theta(i,7)*Y(7,2));
        pp3=theta(i,3)*Y(k,3)/(theta(i,1)*Y(1,3)+theta(i,2)*Y(2,3)+theta(i,3)*Y(3,3)+theta(i,4)*Y(4,3)+theta(i,5)*Y(5,3)+theta(i,6)*Y(6,3)+theta(i,7)*Y(7,3));
        pp4=theta(i,4)*Y(k,4)/(theta(i,1)*Y(1,4)+theta(i,2)*Y(2,4)+theta(i,3)*Y(3,4)+theta(i,4)*Y(4,4)+theta(i,5)*Y(5,4)+theta(i,6)*Y(6,4)+theta(i,7)*Y(7,4));
        pp5=theta(i,5)*Y(k,5)/(theta(i,1)*Y(1,5)+theta(i,2)*Y(2,5)+theta(i,3)*Y(3,5)+theta(i,4)*Y(4,5)+theta(i,5)*Y(5,5)+theta(i,6)*Y(6,5)+theta(i,7)*Y(7,5));
        pp6=theta(i,6)*Y(k,6)/(theta(i,1)*Y(1,6)+theta(i,2)*Y(2,6)+theta(i,3)*Y(3,6)+theta(i,4)*Y(4,6)+theta(i,5)*Y(5,6)+theta(i,6)*Y(6,6)+theta(i,7)*Y(7,6));
        pp7=theta(i,7)*Y(k,7)/(theta(i,1)*Y(1,7)+theta(i,2)*Y(2,7)+theta(i,3)*Y(3,7)+theta(i,4)*Y(4,7)+theta(i,5)*Y(5,7)+theta(i,6)*Y(6,7)+theta(i,7)*Y(7,7));
        LNgammaRP(i,k)=Qk(k)*(1-log(theta(i,1)*Y(1,k)+theta(i,2)*Y(2,k)+theta(i,3)*Y(3,k)+theta(i,4)*Y(4,k)+theta(i,5)*Y(5,k)+theta(i,6)*Y(6,k)+theta(i,7)*Y(7,k))-(pp1+pp2+pp3+pp4+pp5+pp6+pp7));
    end
end

%%group_b fraction parameters: for all components
for i=1:Nc  % components
    for k=1:Nb  % group_bs
        p1Xi(i,k)=v_comp_b(i,k)*Xb(i);  
    end
end
SumNb=sum(p1Xi);
SumNc=sum(SumNb);

for k=1:Nb  % for each group_b
    Xi(k)= SumNb(k)/SumNc;
end

% theta calculation for all components 
    for k=1:Nb  %3
        p1theta(k)=Qk(k)*Xi(k);
    end
    thetaT=sum(p1theta);

    for k=1:Nb  %3
        theta(k)=p1theta(k)/thetaT;  
    end

%LNgamma 
for i=1:Nc
    for k=1:Nb
        p1=theta(1)*Y(k,1)/(theta(1)*Y(1,1)+theta(2)*Y(2,1)+theta(3)*Y(3,1)+theta(4)*Y(4,1)+theta(5)*Y(5,1)+theta(6)*Y(6,1)+theta(7)*Y(7,1));
        p2=theta(2)*Y(k,2)/(theta(1)*Y(1,2)+theta(2)*Y(2,2)+theta(3)*Y(3,2)+theta(4)*Y(4,2)+theta(5)*Y(5,2)+theta(6)*Y(6,2)+theta(7)*Y(7,2));
        p3=theta(3)*Y(k,3)/(theta(1)*Y(1,3)+theta(2)*Y(2,3)+theta(3)*Y(3,3)+theta(4)*Y(4,3)+theta(5)*Y(5,3)+theta(6)*Y(6,3)+theta(7)*Y(7,3));
        p4=theta(4)*Y(k,4)/(theta(1)*Y(1,4)+theta(2)*Y(2,4)+theta(3)*Y(3,4)+theta(4)*Y(4,4)+theta(5)*Y(5,4)+theta(6)*Y(6,4)+theta(7)*Y(7,4));
        p5=theta(5)*Y(k,5)/(theta(1)*Y(1,5)+theta(2)*Y(2,5)+theta(3)*Y(3,5)+theta(4)*Y(4,5)+theta(5)*Y(5,5)+theta(6)*Y(6,5)+theta(7)*Y(7,5));
        p6=theta(6)*Y(k,6)/(theta(1)*Y(1,6)+theta(2)*Y(2,6)+theta(3)*Y(3,6)+theta(4)*Y(4,6)+theta(5)*Y(5,6)+theta(6)*Y(6,6)+theta(7)*Y(7,6));
        p7=theta(7)*Y(k,7)/(theta(1)*Y(1,7)+theta(2)*Y(2,7)+theta(3)*Y(3,7)+theta(4)*Y(4,7)+theta(5)*Y(5,7)+theta(6)*Y(6,7)+theta(7)*Y(7,7));
        LNgammaR(i,k)=Qk(k)*(1-log(theta(1)*Y(1,k)+theta(2)*Y(2,k)+theta(3)*Y(3,k)+theta(4)*Y(4,k)+theta(5)*Y(5,k)+theta(6)*Y(6,k)+theta(7)*Y(7,k))-(p1+p2+p3+p4+p5+p6+p7));
    end
end

% caluclation activity coefficient Residual:
for i=1:Nc
    for k=1:Nb
        LNgammaR_bi(i,k)=group_b(i,k)*v_comp_b(i,k)*(LNgammaR(i,k)-LNgammaRP(i,k));
    end
    LNgammaR_b(i)=sum(LNgammaR_bi(i,:));
    LNgamma_bulk(i)=LNgammaC_b(i)+LNgammaR_b(i);
    Gamma_bulk(i)=exp(LNgamma_bulk(i));
end