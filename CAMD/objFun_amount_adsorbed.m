% Objective Function and Constraints 

function [ObjFunAd res nA1 pAdsorbed]=objFun_amount_adsorbed(TSOL)

NDIM = length(TSOL);    % 2 continuous variables and 7 discrete
m0 = TSOL(1,1);
n1b = TSOL(1,2);
%n1b = .0334;
K2 = TSOL(1,3:NDIM);

[N1,Nb,Xb,Nstar,Xstar,SPC,Gamma_bulk,Gamma_star,mtotal,V,T,R,x1s,phistar]=equations_amount_adsorbed_As(m0,n1b,K2);
n11=N1(1);    


nA1=n11;
nAb=n1b;

F1star=Nstar(1)/mtotal;    

% Nonlinear inequality constraints
c(1)= 0-m0;              % 0<m0>500
c(2)= m0-500;             
c(3)= -F1star;           % 0<F1star< 1.106
%c(4)= F1star-1.106; 
c(4)= F1star-0.5;
%--------boundaries----------
%--bulk phase:
c(5)= Xb(1)-1;   c(6)= Xb(2)-1; c(7)= Xb(3)-1; c(8)= Xb(4)-1;
c(9)= -Xb(1);    c(10)= -Xb(2); c(11)= -Xb(3); c(12)= -Xb(4);   
%-- ASS:
c(13)= Xstar(1)-1; c(14)= Xstar(2)-1;  c(15)= -Xstar(1); c(16)= -Xstar(2);  
c(17)= -n1b; c(18)=n1b-n11;  %%% maximun amount of N1b means no adsorption

% To check the valence of the gruops so the molecule is no charged
[SumVal]=sumValance(K2);
c(19)=SumVal;


% constraints based on the thermodynamic equation:
tol=1e-8;
c(20)=(Xb(1)*Gamma_bulk(1))-(Xstar(1)*Gamma_star(1)*exp(-phistar(1)/(R*T*SPC)));

ff=Xb(1)*Gamma_bulk(1);
qq=Xstar(1)*Gamma_star(1)*exp(-phistar(1)/(R*T*SPC));

for i=1:length(c-2)
    if c(i)<= 0 
        resi(i)=min(0,-c(i));
    else
        resi(i)=c(i);
    end
end

if c(19) ~= 0 
    resi(19)=abs(c(19));
else
    resi(19)=min(0,abs(c(19)));
end


if abs(ff-qq) <= tol 
    resi(20)=min(0,abs(c(20)));
else
    resi(20)=abs(c(20));
end

res=sum(resi);
% objective function 
ObjFunAd=(nA1-nAb)/m0;

pAdsorbed=((nA1-nAb)/nA1)*100;