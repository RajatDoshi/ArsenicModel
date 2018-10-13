function wt=CALWeight(K,q);
% K=6;q=0.001;
 wt = zeros(K,1);
for j = 1:K
%     sum=sum+y(i).^2;
    wt(j) =(1/(q*K*sqrt(2*pi)))*exp(-(j-1)^2/(2*(q*K)^2));
end