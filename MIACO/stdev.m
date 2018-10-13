function [mean sigma]=stdev(y,rho,index);
[K NDIM] = size(y);
sig = zeros(K,NDIM);
sigma = zeros(1,NDIM);
 for j=1:NDIM
     sum(j) = 0;
    for k=1:K
        sig(k,j) = abs((y(k,j)-y(index,j)));
        sum(j) = sum(j)+sig(k,j);
    end
    sigma = rho*(sum/(K-1));
    mean(j) = y(index,j);
 end

