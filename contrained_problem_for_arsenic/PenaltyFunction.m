function P = PenaltyFunction(at);
 acc = 1E-6;
Oracle = 1;

[ SumObjFun res ]=unifac_parameters_ObjFun_all_adsorbents(at);

% P=solute_distribution + ineqcons1 + 100*ineqcons2 + 0.01*ineqcons3 + 0.01*ineqcons4 ;
% P = -800*solute_distribution - (min(0,0.8*ineqcons1) + min(0,90*ineqcons2) + min(0,0.001*ineqcons3) + min(0,0.001*ineqcons4)) ;
%
f = SumObjFun;

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


