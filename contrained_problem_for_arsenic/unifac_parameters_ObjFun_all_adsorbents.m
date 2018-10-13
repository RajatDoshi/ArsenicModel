%---- this program containts the Objective Function of the adsorption
%---problem of Ra into montmorillonite
% Last Modification: April 5/14

function [SumObjFun res]=unifac_parameters_ObjFun_all_adsorbents(at)
T=308.15;            % Temperature  (K)
R=8.314;             % Ideal constant (J/molK)
m0=18.079944;        % Amount of solid for Sepiolite and Beidellite     
m01=20.9100601;      % Amount of solid for Zeolite
for ad=1:3
    [x1_star,Xb,qexp,Gamma_ms1,M0]=input_data_adsorption(ad);
    numex=length(x1_star);          % Amount of data for optimization problems 
 % Gamma for the pure component
    LNgamma01_star=log(1+(1/(Gamma_ms1*M0)));
    for l=1:numex
        x0_star(l)=1-x1_star(l); 
        %-------- Activity coefficient from bulk phase  
        [LNgamma_bulk,Gamma_bulk]=UNIFAC_bulk_phase_all_adsorbents(Xb(l,:),T,at);
        %--------Activity coefficient from Adsorbate-solid solution
   %--Sepiolite:
        if ad==1
            [LNgammaGE_star,gamma_star]=unifac_adsorbate_solid_solution_sepiolite(T,LNgamma01_star,x0_star(l),at);
        end
   %--Beidellite:
        if ad==2
            [LNgammaGE_star,gamma_star]=unifac_adsorbate_solid_solution_beidellite(T,LNgamma01_star,x0_star(l),at);
        end
    %--Zeolite:
        if ad==3
            [LNgammaGE_star,gamma_star]=unifac_adsorbate_solid_solution_zeolite(T,LNgamma01_star,x0_star(l),at);
        end
        % % Calculation of phi*-phi*_oi
        if ad==1 || ad==2
            phi_star_phi_star01(l)=-R*T*Gamma_ms1*log((Gamma_bulk(1)*Xb(l,1))/(gamma_star(1)*x1_star(l)));
            n_star(l)=(phi_star_phi_star01(l)*m0)/(R*T*(x0_star(l)*LNgammaGE_star(2)+x1_star(l)*LNgammaGE_star(1)-x1_star(l)*LNgamma01_star));
            qcal(l,ad)=(n_star(l)*x1_star(l))/m0;
        end 
        if ad==3
            phi_star_phi_star01(l)=-R*T*Gamma_ms1*log((Gamma_bulk(1)*Xb(l,1))/(gamma_star(1)*x1_star(l)));
            n_star(l)=(phi_star_phi_star01(l)*m01)/(R*T*(x0_star(l)*LNgammaGE_star(2)+x1_star(l)*LNgammaGE_star(1)-x1_star(l)*LNgamma01_star));
            qcal(l,ad)=(n_star(l)*x1_star(l))/m01;
        end
        %-----OBJECTIVE FUNCTION 2:
        ObjFuni(l,ad)=abs((qcal(l,ad)-qexp(l))/qexp(l));

        %---- TO CHECK THAT THE SUM OF THE CHEMICAL POTENTIAL IS ZERO
        if ad==1 || ad==2
            GE_star(l)=(x0_star(l)*n_star(l)*LNgammaGE_star(2)+n_star(l)*x1_star(l)*LNgammaGE_star(1));      % Modified Gibbs excess energy
            GE_s=0;                                                                                          % Since x1_s=1 the activity coefficient is 1 and ln(1)=0
            GE_star0i(l)=x1_star(l)*n_star(l)*LNgamma01_star;
            phi_star_phi_star012(l)=(R*T/m0)*(GE_star(l)-GE_s-GE_star0i(l));
            ObjFuni2(l,ad)=(Xb(l,1)*Gamma_bulk(1)-(x1_star(l)*gamma_star(1)*exp(-(phi_star_phi_star012(l))/(R*T*Gamma_ms1))));
        end
        if ad==3
            GE_star(l)=(x0_star(l)*n_star(l)*LNgammaGE_star(2)+n_star(l)*x1_star(l)*LNgammaGE_star(1));      % Modified Gibbs excess energy
            GE_s=0;                                                                                          % Since x1_s=1 the activity coefficient is 1 and ln(1)=0
            GE_star0i(l)=x1_star(l)*n_star(l)*LNgamma01_star;
            phi_star_phi_star012(l)=(R*T/m01)*(GE_star(l)-GE_s-GE_star0i(l));
            ObjFuni2(l,ad)=(Xb(l,1)*Gamma_bulk(1)-(x1_star(l)*gamma_star(1)*exp(-(phi_star_phi_star012(l))/(R*T*Gamma_ms1))));  
        end
        %------------lnGamma Summary-------
        % gammasmax(l,1)=LNgamma_bulk(1);    % Ra
        %  LNgammasmax(l,2)=LNgammaGE_star(1) ; % Ra
        %  LNgammasmax(l,3)=LNgammaGE_star(2);  % solid
    end
    %-----OBJECTIVE FUNCTION 1:     
    %ObjFun=sum(ObjFuni);
    %-----OBJECTIVE FUNCTION 2:
    ObjFun(ad)=sum(ObjFuni(:,ad))/numex;
    SumObjFun=sum(ObjFun);
end 
LB=-50000;    % Lower bound of the interaction parameters
UB= 80000;     % Upper bound of the interaction parameters  
for j=1:length(at)
    ineqcona(j)=at(j)- LB;    % inequality a: at(i)- LB > 0
    ineqconb(j)=-at(j)+ UB;   % inequality b: -at(i)+ UP > 0
    mina(j)=min(0,ineqcona(j));
    minb(j)=min(0,ineqconb(j));
end 

res=-(mina(1)+minb(1)+mina(2)+minb(2)+mina(3)+minb(3)+mina(4)+minb(4)+mina(5)+minb(5)+mina(6)+minb(6)+mina(7)+minb(7)+mina(8)+minb(8)+mina(9)+minb(9)+mina(10)+minb(10)+mina(11)+minb(11)+mina(12)+minb(12)+mina(13)+minb(13)+mina(14)+minb(14)+mina(15)+minb(15)+mina(16)+minb(16)+mina(17)+minb(17)+mina(18)+minb(18)+mina(19)+minb(19)+mina(20)+minb(20)+mina(21)+minb(21)+mina(22)+minb(22)+mina(23)+minb(23)+mina(24)+minb(24)+mina(25)+minb(25)+mina(26)+minb(26)+mina(27)+minb(27)+mina(28)+minb(28)+mina(29)+minb(29)+mina(30)+minb(30)+mina(31)+minb(31)+mina(32)+minb(32)+mina(33)+minb(33)+mina(34)+minb(34)+mina(35)+minb(35)+mina(36)+minb(36)+mina(37)+minb(37)+mina(38)+minb(38)+mina(39)+minb(39)+mina(40)+minb(40)+mina(41)+minb(41)+mina(42)+minb(42)+mina(43)+minb(43)+mina(44)+minb(44));