%Adsorption of As for CAMD 
%This program calculates the Gamma in the Adsorbate - Solid - Solution
function [LNgammaGE_star,gamma_star]=UNIFAC_ASS(T,LNgamma0_star,Xstar,Rk,Qk,Nid,Nindex)

%components in the adsorbate-solid-solution: % As (1) 
%  components:  group #  
%  Radium       As    1
%  Adsorbent    CAMD  0: [ SiO2 Al2O3 Na2O MgO Fe2O3 K20 CaO] % any of this groups 

%------------------Input-----------------------
Nc=length(Xstar);       % Number of components in the mixture (adsorbed solid solution)
Nb=1+length(Nid);       % Number of groups in the adsorbed-solid solution
TNb=length(Qk);        
         
v_comp_star=zeros(Nc,TNb);  
group_star=zeros(Nc,TNb);

v_comp_star(1,1)=1; 
v_comp_star(2,:)=[0 Nindex]; % Total number of groups including As taking into account all possible groups

% for matrix group_star 1= yes 0= no yes: the component has the first group.
for h=1:Nc
    for hh=1:TNb
        if v_comp_star(h,hh) ~=0
        group_star(h,hh)=1;
        else 
        group_star(h,hh)=0;
        end
    end
end


%------COMBINATORIAL-----------------------
% for xi*:
for i=1:Nc      % per component
    for k=1:TNb  % per gruop 
        rr(i,k)=v_comp_star(i,k)*Rk(k);
        qq(i,k)=v_comp_star(i,k)*Qk(k);
    end
    ri(i)=sum(rr(i,:));
    qi(i)=sum(qq(i,:));
end
for jj=1:Nc
 ri_Xstar(jj)=ri(jj)*Xstar(jj);
 qi_Xstar(jj)=qi(jj)*Xstar(jj);
end

for i=1:Nc
    V(i)=ri(i)/sum(ri_Xstar);
    F(i)=qi(i)/sum(qi_Xstar);
    LNGammaGE_star_C(i)=1-V(i)+log(V(i))-5*qi(i)*(1-(V(i)/F(i))+log(V(i)/F(i)));
end

%---------RESIDUAL----------------------%
% group interaction parameters for the residual coefficient
Gyes=sum(group_star);  %to determine if the gruop is included here or no yes=1 no=0  
Anm=interaction_parameter_ASS(Gyes);  % the column or rowa filled of zeros means that that group is not included


for n=1:TNb %Nb
    for m=1:TNb%Nb
        Y(n,m)=exp(-Anm(n,m)/T);
    end
end

%---group fraction parameters for pure component:
for i=1:Nc
    for k=1:TNb%Nb
        X0i(i,k)=v_comp_star(i,k)/sum(v_comp_star(i,:));
    end
end

%---group surface area fraction for pure component:
for i=1:Nc
    for k=1:TNb%Nb
        theta0ii(i,k)=Qk(k)*X0i(i,k);
    end
    for kk=1:TNb%Nb
        theta0i(i,kk)=theta0ii(i,kk)/sum(theta0ii(i,:));  
    end          
 end
%---gamma for pure components, 
for i=1:Nc 
    for k=1:TNb%Nb 
        p01=theta0i(i,1)*Y(k,1)/(theta0i(i,1)*Y(1,1)+theta0i(i,2)*Y(2,1)+theta0i(i,3)*Y(3,1)+theta0i(i,4)*Y(4,1)+theta0i(i,5)*Y(5,1)+theta0i(i,6)*Y(6,1)+theta0i(i,7)*Y(7,1)+theta0i(i,8)*Y(8,1)+theta0i(i,9)*Y(9,1));
        p02=theta0i(i,2)*Y(k,2)/(theta0i(i,1)*Y(1,2)+theta0i(i,2)*Y(2,2)+theta0i(i,3)*Y(3,2)+theta0i(i,4)*Y(4,2)+theta0i(i,5)*Y(5,2)+theta0i(i,6)*Y(6,2)+theta0i(i,7)*Y(7,2)+theta0i(i,8)*Y(8,2)+theta0i(i,9)*Y(9,2));
        p03=theta0i(i,3)*Y(k,3)/(theta0i(i,1)*Y(1,3)+theta0i(i,2)*Y(2,3)+theta0i(i,3)*Y(3,3)+theta0i(i,4)*Y(4,3)+theta0i(i,5)*Y(5,3)+theta0i(i,6)*Y(6,3)+theta0i(i,7)*Y(7,3)+theta0i(i,8)*Y(8,3)+theta0i(i,9)*Y(9,3));
        p04=theta0i(i,4)*Y(k,4)/(theta0i(i,1)*Y(1,4)+theta0i(i,2)*Y(2,4)+theta0i(i,3)*Y(3,4)+theta0i(i,4)*Y(4,4)+theta0i(i,5)*Y(5,4)+theta0i(i,6)*Y(6,4)+theta0i(i,7)*Y(7,4)+theta0i(i,8)*Y(8,4)+theta0i(i,9)*Y(9,4));
        p05=theta0i(i,5)*Y(k,5)/(theta0i(i,1)*Y(1,5)+theta0i(i,2)*Y(2,5)+theta0i(i,3)*Y(3,5)+theta0i(i,4)*Y(4,5)+theta0i(i,5)*Y(5,5)+theta0i(i,6)*Y(6,5)+theta0i(i,7)*Y(7,5)+theta0i(i,8)*Y(8,5)+theta0i(i,9)*Y(9,5));
        p06=theta0i(i,6)*Y(k,6)/(theta0i(i,1)*Y(1,6)+theta0i(i,2)*Y(2,6)+theta0i(i,3)*Y(3,6)+theta0i(i,4)*Y(4,6)+theta0i(i,5)*Y(5,6)+theta0i(i,6)*Y(6,6)+theta0i(i,7)*Y(7,6)+theta0i(i,8)*Y(8,6)+theta0i(i,9)*Y(9,6));
        p07=theta0i(i,7)*Y(k,7)/(theta0i(i,1)*Y(1,7)+theta0i(i,2)*Y(2,7)+theta0i(i,3)*Y(3,7)+theta0i(i,4)*Y(4,7)+theta0i(i,5)*Y(5,7)+theta0i(i,6)*Y(6,7)+theta0i(i,7)*Y(7,7)+theta0i(i,8)*Y(8,7)+theta0i(i,9)*Y(9,7));
        p08=theta0i(i,8)*Y(k,8)/(theta0i(i,1)*Y(1,8)+theta0i(i,2)*Y(2,8)+theta0i(i,3)*Y(3,8)+theta0i(i,4)*Y(4,8)+theta0i(i,5)*Y(5,8)+theta0i(i,6)*Y(6,8)+theta0i(i,7)*Y(7,8)+theta0i(i,8)*Y(8,8)+theta0i(i,9)*Y(9,8));
        p09=theta0i(i,9)*Y(k,9)/(theta0i(i,1)*Y(1,9)+theta0i(i,2)*Y(2,9)+theta0i(i,3)*Y(3,9)+theta0i(i,4)*Y(4,9)+theta0i(i,5)*Y(5,9)+theta0i(i,6)*Y(6,9)+theta0i(i,7)*Y(7,9)+theta0i(i,9)*Y(8,9)+theta0i(i,9)*Y(9,9));
        LNgammaR0_star(i,k)=group_star(i,k)*Qk(k)*(1-log(theta0i(i,1)*Y(1,k)+theta0i(i,2)*Y(2,k)+theta0i(i,3)*Y(3,k)+theta0i(i,4)*Y(4,k)+theta0i(i,5)*Y(5,k)+theta0i(i,6)*Y(6,k)+theta0i(i,7)*Y(7,k)+theta0i(i,8)*Y(8,k)+theta0i(i,9)*Y(9,k))-(p01+p02+p03+p04+p05+p06+p07+p08+p09));
    end
end
   
%---group fraction parameters: for all components
for i=1:Nc  
    for k=1:TNb %Nb 
        p1Xi(i,k)=v_comp_star(i,k)*Xstar(i);  
    end
end

SumNb=sum(p1Xi);
SumNc=sum(SumNb);

for k=1:TNb %Nb  % for each group
    Xi(k)= SumNb(k)/SumNc;
end
%---theta calculation for all components 
for k=1:TNb %Nb  
   thetaii(k)=Qk(k)*Xi(k);
end
thetaiT=sum(thetaii);

 for k=1:TNb % Nb
     
     thetai(k)=thetaii(k)/thetaiT; 
 end

for k=1:TNb % Nb
    p1=thetai(1)*Y(k,1)/(thetai(1)*Y(1,1)+thetai(2)*Y(2,1)+thetai(3)*Y(3,1)+thetai(4)*Y(4,1)+thetai(5)*Y(5,1)+thetai(6)*Y(6,1)+thetai(7)*Y(7,1)+thetai(8)*Y(8,1)+thetai(9)*Y(9,1));
    p2=thetai(2)*Y(k,2)/(thetai(1)*Y(1,2)+thetai(2)*Y(2,2)+thetai(3)*Y(3,2)+thetai(4)*Y(4,2)+thetai(5)*Y(5,2)+thetai(6)*Y(6,2)+thetai(7)*Y(7,2)+thetai(8)*Y(8,2)+thetai(9)*Y(9,2));
    p3=thetai(3)*Y(k,3)/(thetai(1)*Y(1,3)+thetai(2)*Y(2,3)+thetai(3)*Y(3,3)+thetai(4)*Y(4,3)+thetai(5)*Y(5,3)+thetai(6)*Y(6,3)+thetai(7)*Y(7,3)+thetai(8)*Y(8,3)+thetai(9)*Y(9,3));
    p4=thetai(4)*Y(k,4)/(thetai(1)*Y(1,4)+thetai(2)*Y(2,4)+thetai(3)*Y(3,4)+thetai(4)*Y(4,4)+thetai(5)*Y(5,4)+thetai(6)*Y(6,4)+thetai(7)*Y(7,4)+thetai(8)*Y(8,4)+thetai(9)*Y(9,4));
    p5=thetai(5)*Y(k,5)/(thetai(1)*Y(1,5)+thetai(2)*Y(2,5)+thetai(3)*Y(3,5)+thetai(4)*Y(4,5)+thetai(5)*Y(5,5)+thetai(6)*Y(6,5)+thetai(7)*Y(7,5)+thetai(8)*Y(8,5)+thetai(9)*Y(9,5));
    p6=thetai(6)*Y(k,6)/(thetai(1)*Y(1,6)+thetai(2)*Y(2,6)+thetai(3)*Y(3,6)+thetai(4)*Y(4,6)+thetai(5)*Y(5,6)+thetai(6)*Y(6,6)+thetai(7)*Y(7,6)+thetai(8)*Y(8,6)+thetai(9)*Y(9,6));
    p7=thetai(7)*Y(k,7)/(thetai(1)*Y(1,7)+thetai(2)*Y(2,7)+thetai(3)*Y(3,7)+thetai(4)*Y(4,7)+thetai(5)*Y(5,7)+thetai(6)*Y(6,7)+thetai(7)*Y(7,7)+thetai(8)*Y(8,7)+thetai(9)*Y(9,7));
    p8=thetai(8)*Y(k,8)/(thetai(1)*Y(1,8)+thetai(2)*Y(2,8)+thetai(3)*Y(3,8)+thetai(4)*Y(4,8)+thetai(5)*Y(5,8)+thetai(6)*Y(6,8)+thetai(7)*Y(7,8)+thetai(8)*Y(8,8)+thetai(9)*Y(9,8));
    p9=thetai(9)*Y(k,9)/(thetai(1)*Y(1,9)+thetai(2)*Y(2,9)+thetai(3)*Y(3,9)+thetai(4)*Y(4,9)+thetai(5)*Y(5,9)+thetai(6)*Y(6,9)+thetai(7)*Y(7,9)+thetai(8)*Y(8,9)+thetai(9)*Y(9,9));
    LNgammaR_star(k)=Qk(k)*(1-log(thetai(1)*Y(1,k)+thetai(2)*Y(2,k)+thetai(3)*Y(3,k)+thetai(4)*Y(4,k)+thetai(5)*Y(5,k)+thetai(6)*Y(6,k)+thetai(7)*Y(7,k)+thetai(8)*Y(8,k)+thetai(9)*Y(9,k))-(p1+p2+p3+p4+p5+p6+p7+p8+p9));
end

% caluclation activity coefficient Residual:
for i=1:Nc
    for k=1:TNb %Nb
        LNgammaGE_star_R_i(i,k)=group_star(i,k)*v_comp_star(i,k)*(LNgammaR_star(k)-LNgammaR0_star(i,k));
    end
    LNGammaGE_star_R(i)=sum(LNgammaGE_star_R_i(i,:));
end
LNgammaGE_star(1)=LNGammaGE_star_C(1)+LNGammaGE_star_R(1); % Adsorbate As
LNgammaGE_star(2)=LNGammaGE_star_C(2)+LNGammaGE_star_R(2); % Adsorbent
 
gamma_star(1)=exp(LNgammaGE_star(1));
gamma_star(2)=exp(LNgammaGE_star(2));


