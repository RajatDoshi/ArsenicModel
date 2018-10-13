function [LNgamma_bulk,Gamma_bulk]=UNIFAC_bulk_phase(Xb,T)

%Input-----------------------
v_comp_b=[1 0 0 0; %Arsenic
          0 1 0 0; %Na
          0 0 1 0; %HCO3
          0 0 0 1;]; %H20
         

% Van der Waals volume and surface parameters 
Rk=[3.61 3 0.631 0.92];
Qk=[5.28 3 0.737 1.4];

Nc=length(Xb);  % number of components 
Nb=length(Qk);  % number of groups in the bulk phase
group_b=[1 0 0 0 ;0 1 0 1; 0 0 1 0;0 0 0 1 ]; % 1 = yes 0= no yes: the component has the first group.

% ------COMBINATORIAL-----------------------
for i=1:Nc      % per component
    for k=1:Nb  % per gruop 
        rr(i,k)=v_comp_b(i,k)*Rk(k);
        qq(i,k)=v_comp_b(i,k)*Qk(k);
    end
    ri(i)=sum(rr(i,:));
    qi(i)=sum(qq(i,:));
end
for ii=1:Nc
Sri_Xstari(ii)=ri(ii)*Xb(ii);
Sqi_Xstari(ii)=qi(ii)*Xb(ii);
end
Sri_Xstar=sum(Sri_Xstari);
Sqi_Xstar=sum(Sqi_Xstari);

for i=1:Nc
    V(i)=ri(i)/Sri_Xstar;
    F(i)=qi(i)/Sqi_Xstar;
    LNgammaC_b(i)=1-V(i)+log(V(i))-5*qi(i)*(1-(V(i)/F(i))+log(V(i)/F(i)));
end

%------ RESIDUAL--------------------------------------------
%activity coefficient parameters bulk phase
%AS 
    a=[0    6857.10    2715.10      -7914.10;
   1080.30  0          6342.2       -165;   
   2086.80  10798.0    0            -24755.0;
  758.38   22.38      -982.5       0;];
  
for n=1:Nb
    for m=1:Nb
        Y(n,m)=exp(-a(n,m)/T);
    end
end

% for the pure component i - group_b fraction parameters:
%---group fraction parameters for pure component:
for i=1:Nc
    for k=1:Nb
        X0i(i,k)=v_comp_b(i,k)/sum(v_comp_b(i,:));
    end
end

% group_b surface area fraction pure component:
for i=1:Nc
    for k=1:Nb
        theta0i(i,k)=Qk(k)*X0i(i,k);
    end
    for kk=1:Nb
        theta0(i,kk)=theta0i(i,kk)/sum(theta0i(i,:));  %x11
    end          
end

% gamma for pure component i
for i=1:Nc 
    for k=1:Nb 
        pp1=theta0(i,1)*Y(k,1)/(theta0(i,1)*Y(1,1)+theta0(i,2)*Y(2,1)+theta0(i,3)*Y(3,1)+theta0(i,4)*Y(4,1));
        pp2=theta0(i,2)*Y(k,2)/(theta0(i,1)*Y(1,2)+theta0(i,2)*Y(2,2)+theta0(i,3)*Y(3,2)+theta0(i,4)*Y(4,2));
        pp3=theta0(i,3)*Y(k,3)/(theta0(i,1)*Y(1,3)+theta0(i,2)*Y(2,3)+theta0(i,3)*Y(3,3)+theta0(i,4)*Y(4,3));
        pp4=theta0(i,4)*Y(k,4)/(theta0(i,1)*Y(1,4)+theta0(i,2)*Y(2,4)+theta0(i,3)*Y(3,4)+theta0(i,4)*Y(4,4));
        LNgammaRP(i,k)=group_b(i,k)*Qk(k)*(1-log(theta0(i,1)*Y(1,k)+theta0(i,2)*Y(2,k)+theta0(i,3)*Y(3,k)+theta0(i,4)*Y(4,k))-(pp1+pp2+pp3+pp4));
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

for k=1:Nb
    p1=theta(1)*Y(k,1)/(theta(1)*Y(1,1)+theta(2)*Y(2,1)+theta(3)*Y(3,1)+theta(4)*Y(4,1));
    p2=theta(2)*Y(k,2)/(theta(1)*Y(1,2)+theta(2)*Y(2,2)+theta(3)*Y(3,2)+theta(4)*Y(4,2));
    p3=theta(3)*Y(k,3)/(theta(1)*Y(1,3)+theta(2)*Y(2,3)+theta(3)*Y(3,3)+theta(4)*Y(4,3));
    p4=theta(4)*Y(k,4)/(theta(1)*Y(1,4)+theta(2)*Y(2,4)+theta(3)*Y(3,4)+theta(4)*Y(4,4));
    LNgammaR(k)=Qk(k)*(1-log(theta(1)*Y(1,k)+theta(2)*Y(2,k)+theta(3)*Y(3,k)+theta(4)*Y(4,k))-(p1+p2+p3+p4));
end


% caluclation activity coefficient Residual:
for i=1:Nc
    for k=1:Nb
        LNgammaR_bi(i,k)=group_b(i,k)*v_comp_b(i,k)*(LNgammaR(k)-LNgammaRP(i,k));
    end
    LNgammaR_b(i)=sum(LNgammaR_bi(i,:));
    LNgamma_bulk(i)=LNgammaC_b(i)+LNgammaR_b(i);
    Gamma_bulk(i)=exp(LNgamma_bulk(i));
end
