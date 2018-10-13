function [LNgammaGE_star,gamma_star]=unifac_adsorbate_solid_solution_zeolite(T,LNgamma01_star,x0_star,at)

%groups: Ra SiO2 Al2O3 Na2O3 MgO Fe3O2 
Rk=[3.61 1.389 1.688 1.109 1.752 2.716 1.013 1.221];
Qk=[5.28 1.267 1.407 1.103 1.438 1.844 1.041 1.167];

%Input-----------------------
v_comp_star=[1 0 0 0 0 0 0 0; 0 1.1951 0.1187 0.0196 0.0112 0.0032 0.0617 0.0251;];
group_star=[1 0 0 0 0 0 0 0; 0 1 1 1 1 1 1 1;];
x1_star=(1-x0_star); 
x1_s=1;
x1_s=x1_star/(1-x0_star);           
X_star=[x1_star x0_star];           

Nc=length(X_star);                  % number of components in the mixture (adsorbed solid solution)
                                    % 1:Ra 2:solid which is represented by Zero in the paper
Nb=length(Qk);                      % number of groups in the adsorbed-solid solution

%------COMBINATORIAL-----------------------
for i=1:Nc      % per component
    for k=1:Nb  % per gruop 
        rr(i,k)=v_comp_star(i,k)*Rk(k);
        qq(i,k)=v_comp_star(i,k)*Qk(k);
    end
    ri(i)=sum(rr(i,:));
    qi(i)=sum(qq(i,:));
end
Sri_Xstar=ri(1)*X_star(1)+ ri(2)*X_star(2);
Sqi_Xstar=qi(1)*X_star(1)+ qi(2)*X_star(2);
 
% calculation from the paper: 
for i=1:Nc
    V(i)=ri(i)/Sri_Xstar;
    F(i)=qi(i)/Sqi_Xstar;
    LNGammaGE_star_C(i)=1-V(i)+log(V(i))-5*qi(i)*(1-(V(i)/F(i))+log(V(i)/F(i)));
end
 
%---------RESIDUAL----------------------%
% group interaction parameters for the residual coefficient
aa12=at(31); aa13=at(32); aa14=at(33); aa15=at(34); aa16=at(35); aa17=at(36); aa18=at(37); aa21=at(38); aa31=at(39); aa41=at(40); aa51=at(41); aa61=at(42); aa71=at(43); aa81=at(44);
aa=[0 aa12 aa13 aa14 aa15 aa16 aa17 aa18; aa21 0 0 0 0 0 0 0; aa31 0 0 0 0 0 0 0; aa41 0 0 0 0 0 0 0; aa51 0 0 0 0 0 0 0; aa61 0 0 0 0 0 0 0; aa71 0 0 0 0 0 0 0; aa81 0 0 0 0 0 0 0];

for n=1:Nb
    for m=1:Nb
        Y(n,m)=exp(-aa(n,m)/T);
    end
end

%---group fraction parameters for pure component:
for i=1:Nc
    for k=1:Nb
        X0i(i,k)=v_comp_star(i,k)/sum(v_comp_star(i,:));
    end
end
%---group surface area fraction for pure component:
for i=1:Nc
    for k=1:Nb
        theta0ii(i,k)=Qk(k)*X0i(i,k);
    end
    for kk=1:Nb
        theta0i(i,kk)=theta0ii(i,kk)/sum(theta0ii(i,:));  
    end          
end
%---gamma for pure components, 
for i=1:Nc 
    for k=1:Nb 
        p01=theta0i(i,1)*Y(k,1)/(theta0i(i,1)*Y(1,1)+theta0i(i,2)*Y(2,1)+theta0i(i,3)*Y(3,1)+theta0i(i,4)*Y(4,1)+theta0i(i,5)*Y(5,1)+theta0i(i,6)*Y(6,1)+theta0i(i,7)*Y(7,1)+theta0i(i,8)*Y(8,1));
        p02=theta0i(i,2)*Y(k,2)/(theta0i(i,1)*Y(1,2)+theta0i(i,2)*Y(2,2)+theta0i(i,3)*Y(3,2)+theta0i(i,4)*Y(4,2)+theta0i(i,5)*Y(5,2)+theta0i(i,6)*Y(6,2)+theta0i(i,7)*Y(7,2)+theta0i(i,8)*Y(8,2));
        p03=theta0i(i,3)*Y(k,3)/(theta0i(i,1)*Y(1,3)+theta0i(i,2)*Y(2,3)+theta0i(i,3)*Y(3,3)+theta0i(i,4)*Y(4,3)+theta0i(i,5)*Y(5,3)+theta0i(i,6)*Y(6,3)+theta0i(i,7)*Y(7,3)+theta0i(i,8)*Y(8,3));
        p04=theta0i(i,4)*Y(k,4)/(theta0i(i,1)*Y(1,4)+theta0i(i,2)*Y(2,4)+theta0i(i,3)*Y(3,4)+theta0i(i,4)*Y(4,4)+theta0i(i,5)*Y(5,4)+theta0i(i,6)*Y(6,4)+theta0i(i,7)*Y(7,4)+theta0i(i,8)*Y(8,4));
        p05=theta0i(i,5)*Y(k,5)/(theta0i(i,1)*Y(1,5)+theta0i(i,2)*Y(2,5)+theta0i(i,3)*Y(3,5)+theta0i(i,4)*Y(4,5)+theta0i(i,5)*Y(5,5)+theta0i(i,6)*Y(6,5)+theta0i(i,7)*Y(7,5)+theta0i(i,8)*Y(8,5));
        p06=theta0i(i,6)*Y(k,6)/(theta0i(i,1)*Y(1,6)+theta0i(i,2)*Y(2,6)+theta0i(i,3)*Y(3,6)+theta0i(i,4)*Y(4,6)+theta0i(i,5)*Y(5,6)+theta0i(i,6)*Y(6,6)+theta0i(i,7)*Y(7,6)+theta0i(i,8)*Y(8,6));
        p07=theta0i(i,7)*Y(k,7)/(theta0i(i,1)*Y(1,7)+theta0i(i,2)*Y(2,7)+theta0i(i,3)*Y(3,7)+theta0i(i,4)*Y(4,7)+theta0i(i,5)*Y(5,7)+theta0i(i,6)*Y(6,7)+theta0i(i,7)*Y(7,7)+theta0i(i,8)*Y(8,7));
        p08=theta0i(i,8)*Y(k,8)/(theta0i(i,1)*Y(1,8)+theta0i(i,2)*Y(2,8)+theta0i(i,3)*Y(3,8)+theta0i(i,4)*Y(4,8)+theta0i(i,5)*Y(5,8)+theta0i(i,6)*Y(6,8)+theta0i(i,7)*Y(7,8)+theta0i(i,8)*Y(8,8));
        LNgammaR0_star(i,k)=group_star(i,k)*Qk(k)*(1-log(theta0i(i,1)*Y(1,k)+theta0i(i,2)*Y(2,k)+theta0i(i,3)*Y(3,k)+theta0i(i,4)*Y(4,k)+theta0i(i,5)*Y(5,k)+theta0i(i,6)*Y(6,k)+theta0i(i,7)*Y(7,k)+theta0i(i,8)*Y(8,k))-(p01+p02+p03+p04+p05+p06+p07+p08));
    end
end

%---group fraction parameters: for all components
for i=1:Nc  
    for k=1:Nb 
        p1Xi(i,k)=v_comp_star(i,k)*X_star(i);  
    end
end
SumNb1=sum(p1Xi(:,1));  % As
SumNb2=sum(p1Xi(:,2));  % SiO2
SumNb3=sum(p1Xi(:,3));  % Na2O
SumNb4=sum(p1Xi(:,4));  % Al2O3
SumNb5=sum(p1Xi(:,5));  % CaO
SumNb6=sum(p1Xi(:,6));  % K2O
SumNb7=sum(p1Xi(:,7));  % MgO
SumNb8=sum(p1Xi(:,8));  % Fe2O3
SumNb=[SumNb1 SumNb2 SumNb3 SumNb4 SumNb5 SumNb6 SumNb7 SumNb8];
SumNc=sum(SumNb);

% SumNb=sum(p1Xi);
% SumNc=sum(SumNb);
% 
for k=1:Nb  % for each group
    Xi(k)= SumNb(k)/SumNc;
end
%---theta calculation for all components 
for k=1:Nb  
   thetaii(k)=Qk(k)*Xi(k);
end
thetaiT=sum(thetaii);

 for k=1:Nb  %for groups: 6
     thetai(k)=thetaii(k)/thetaiT; 
 end

for k=1:Nb
    p1=thetai(1)*Y(k,1)/(thetai(1)*Y(1,1)+thetai(2)*Y(2,1)+thetai(3)*Y(3,1)+thetai(4)*Y(4,1)+thetai(5)*Y(5,1)+thetai(6)*Y(6,1)+thetai(7)*Y(7,1)+thetai(8)*Y(8,1));
    p2=thetai(2)*Y(k,2)/(thetai(1)*Y(1,2)+thetai(2)*Y(2,2)+thetai(3)*Y(3,2)+thetai(4)*Y(4,2)+thetai(5)*Y(5,2)+thetai(6)*Y(6,2)+thetai(7)*Y(7,2)+thetai(8)*Y(8,2));
    p3=thetai(3)*Y(k,3)/(thetai(1)*Y(1,3)+thetai(2)*Y(2,3)+thetai(3)*Y(3,3)+thetai(4)*Y(4,3)+thetai(5)*Y(5,3)+thetai(6)*Y(6,3)+thetai(7)*Y(7,3)+thetai(8)*Y(8,3));
    p4=thetai(4)*Y(k,4)/(thetai(1)*Y(1,4)+thetai(2)*Y(2,4)+thetai(3)*Y(3,4)+thetai(4)*Y(4,4)+thetai(5)*Y(5,4)+thetai(6)*Y(6,4)+thetai(7)*Y(7,4)+thetai(8)*Y(8,4));
    p5=thetai(5)*Y(k,5)/(thetai(1)*Y(1,5)+thetai(2)*Y(2,5)+thetai(3)*Y(3,5)+thetai(4)*Y(4,5)+thetai(5)*Y(5,5)+thetai(6)*Y(6,5)+thetai(7)*Y(7,5)+thetai(8)*Y(8,5));
    p6=thetai(6)*Y(k,6)/(thetai(1)*Y(1,6)+thetai(2)*Y(2,6)+thetai(3)*Y(3,6)+thetai(4)*Y(4,6)+thetai(5)*Y(5,6)+thetai(6)*Y(6,6)+thetai(7)*Y(7,6)+thetai(8)*Y(8,6));
    p7=thetai(7)*Y(k,7)/(thetai(1)*Y(1,7)+thetai(2)*Y(2,7)+thetai(3)*Y(3,7)+thetai(4)*Y(4,7)+thetai(5)*Y(5,7)+thetai(6)*Y(6,7)+thetai(7)*Y(7,7)+thetai(8)*Y(8,7));
    p8=thetai(8)*Y(k,8)/(thetai(1)*Y(1,8)+thetai(2)*Y(2,8)+thetai(3)*Y(3,8)+thetai(4)*Y(4,8)+thetai(5)*Y(5,8)+thetai(6)*Y(6,8)+thetai(7)*Y(7,8)+thetai(8)*Y(8,8));    
    LNgammaR_star(k)=Qk(k)*(1-log(thetai(1)*Y(1,k)+thetai(2)*Y(2,k)+thetai(3)*Y(3,k)+thetai(4)*Y(4,k)+thetai(5)*Y(5,k)+thetai(6)*Y(6,k)+thetai(7)*Y(7,k)+thetai(8)*Y(8,k))-(p1+p2+p3+p4+p5+p6+p7+p8));
end

% caluclation activity coefficient Residual:
for i=1:Nc
    for k=1:Nb
        LNgammaGE_star_R_i(i,k)=group_star(i,k)*v_comp_star(i,k)*(LNgammaR_star(k)-LNgammaR0_star(i,k));
    end
    LNGammaGE_star_R(i)=sum(LNgammaGE_star_R_i(i,:));
end
%LNgammaGE_star(1)=LNgamma01_star+LnGammaGE_star_C(1)+LNgammaGE_star_R(1);  % adsorbate Ra
LNgammaGE_star(1)=LNGammaGE_star_C(1)+ LNGammaGE_star_R(1);       
LNgammaGE_star(2)= LNGammaGE_star_C(2)+ LNGammaGE_star_R(2);            % adsorbent  
gamma_star(1)=exp(LNgammaGE_star(1));
gamma_star(2)=exp(LNgammaGE_star(2));