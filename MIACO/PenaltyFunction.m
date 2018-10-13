function P = PenaltyFunction(TSOL);
 acc = 1E-6;   % tolerance
% acc = 0 ;
Oracle = 0.1;   % this vaslues greater than the objective funct
%--------------------------------------------------------------------------
%
%
%
%
%--------------------------------------------------------------------------
% [solute_distribution res]= UNIFAC(N);

% [solute_distribution selectivity  solvent_loss    res    T_bp  Mw_S ineqcons1 ineqcons2 ineqcons3 ineqcons4 eqcons1]= UNIFAC(N) ;

% [solute_distribution selectivity  solvent_loss    res    T_bp  Mw_S  ineqcons1  ineqcons2  ineqcons3  ineqcons4 eqcons1]= UNIFAC(N) ;

% P=solute_distribution + ineqcons1 + 100*ineqcons2 + 0.01*ineqcons3 + 0.01*ineqcons4 ;
% P = -800*solute_distribution - (min(0,0.8*ineqcons1) + min(0,90*ineqcons2) + min(0,0.001*ineqcons3) + min(0,0.001*ineqcons4)) ;
%
[ObjFunAd res]=objFun_amount_adsorbed(TSOL);

 f= -ObjFunAd ;


if (f <= Oracle & res <= acc)
    P =-abs(f - Oracle);
elseif (f > Oracle & res <= acc)
    P =abs(f - Oracle);
elseif (f <= Oracle & res > acc)
    P = res;
else
    alpha = 0;
    if (res < abs(f-Oracle)/3)
        alpha = (abs(f-Oracle)*((6*sqrt(3)-2)/(6*sqrt(3)))-res)/(abs(f-Oracle)-res);
        
    elseif (res >=abs(f-Oracle)/3 & res <=abs(f-Oracle))
        alpha = 1- 1/sqrt(abs(f-Oracle)/res);
        
    elseif (res > abs(f-Oracle))
        alpha = (1/2)*sqrt(abs(f-Oracle)/res);
    end
    P = alpha*abs(f-Oracle) + (1-alpha)*res;
%     if (f < Oracle & res <= acc)
%         Oracle = f
%     else
%         Oracle = Oracle
%     end
    
end
end


