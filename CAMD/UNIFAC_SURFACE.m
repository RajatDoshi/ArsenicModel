function [LNgamma_s,Gamma_s]=UNIFAC_SURFACE(Xs,T)

%components in the surface phase : % As (1) 
%  components:  group #  
%  Radium       As    1


%------------------Input-----------------------
Nc=length(Xs);       % Number of components in the surface phase
Nb=length(Xs);       % Number of groups in the surface phase 
         
v_comp_s=zeros(Nc,Nb);  
group_s=zeros(Nc,Nb);

v_comp_s(1,1)=1;  
group_s=v_comp_s;  %for matrix group_star 1= yes 0= no yes: the component has the first group.
Rk=[3.61];
Qk=[5.28];

%------COMBINATORIAL-----------------------
% for xis:
for i=1:Nc      % per component
    for k=1:Nb  % per gruop 
        rr(i,k)=v_comp_s(i,k)*Rk(k);
        qq(i,k)=v_comp_s(i,k)*Qk(k);
    end
    ri(i)=sum(rr(i,:));
    qi(i)=sum(qq(i,:));
end
for jj=1:Nc
 ri_Xs(jj)=ri(jj)*Xs(jj);
 qi_Xs(jj)=qi(jj)*Xs(jj);
end

for i=1:Nc
    V(i)=ri(i)/sum(ri_Xs);
    F(i)=qi(i)/sum(qi_Xs);
    LNGamma_s_C(i)=1-V(i)+log(V(i))-5*qi(i)*(1-(V(i)/F(i))+log(V(i)/F(i)));
end

%---------RESIDUAL----------------------%
% group interaction parameters for the residual coefficient
Anm=[0];

for n=1:Nb
    for m=1:Nb
        Y(n,m)=exp(-Anm(n,m)/T);
    end
end

%---group fraction parameters for pure component:
for i=1:Nc
    for k=1:Nb
        X0i(i,k)=v_comp_s(i,k)/sum(v_comp_s(i,:));
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
        p01=theta0i(i,1)*Y(k,1)/(theta0i(i,1)*Y(1,1));
        LNgammaR0_s(i,k)=group_s(i,k)*Qk(k)*(1-log(theta0i(i,1)*Y(1,k))-(p01));
    end
end

%---group fraction parameters: for all components
for i=1:Nc  
    for k=1:Nb  
        p1Xi(i,k)=v_comp_s(i,k)*Xs(i);  
    end
end

SumNb=sum(p1Xi);
SumNc=sum(SumNb);

for k=1:Nb   % for each group
    Xi(k)= SumNb(k)/SumNc;
end
%---theta calculation for all components 
for k=1:Nb   
   thetaii(k)=Qk(k)*Xi(k);
end
thetaiT=sum(thetaii);

 for k=1:Nb 
     thetai(k)=thetaii(k)/thetaiT; 
 end

for k=1:Nb % Nb
    p1=thetai(1)*Y(k,1)/(thetai(1)*Y(1,1));
    LNgammaR_s(k)=Qk(k)*(1-log(thetai(1)*Y(1,k))-(p1));
end

% caluclation activity coefficient Residual:
for i=1:Nc
    for k=1:Nb
        LNgamma_s_R_i(i,k)=group_s(i,k)*v_comp_s(i,k)*(LNgammaR_s(k)-LNgammaR0_s(i,k));
    end
    LNGamma_s_R(i)=sum(LNgamma_s_R_i(i,:));
end
LNgamma_s(1)=LNGamma_s_C(1)+LNGamma_s_R(1);  % Adsorbate As
Gamma_s(1)=exp(LNgamma_s(1));

