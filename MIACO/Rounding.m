function RSA=Rounding(SA,lbound,ubound)
[M NP]=size(SA);
% ds=(smax-smin)/(NP-1);
ds=1;
npoints =(ubound-lbound)/ds + 1;
s(1) = lbound;
for i=2:npoints
    s(i)=s(i-1)+ds;
end
for i = 1:M
for j = 1:NP
    if (SA(i,j)<= (s(1) +ds/2))
        RSA(i,j) = s(1);
    end
    for k = 2:npoints
        if((s(k-1) +ds/2) < SA(i,j) & SA(i,j) <= (s(k) +ds/2))
            RSA(i,j)=s(k);
        elseif (SA(i,j) > ubound )
            RSA(i,j) = ubound;
        elseif (SA(i,j) < lbound)
            RSA(i,j) = lbound;
        end
    end
end
end
% dels
end

