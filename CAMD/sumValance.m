function [SumVal]=sumValance(K2)
LK=length(K2);
valvector=[0 0 0 0 0 0 0 1];   %[SiO2 Al2O3 Na2O MgO Fe2O3 K20 CaO]
val=[0 0 0 0 0 0 0 0];        % pre allocating

LK=length(K2);
for i=1:LK
  val(i)=valvector(i)*K2(i);
end
SumVal=sum(val);

