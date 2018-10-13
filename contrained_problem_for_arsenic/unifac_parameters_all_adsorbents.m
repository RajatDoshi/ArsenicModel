clc
clear all

%options = optimset(options,'TolFun', 1e-12);
a=[-7.9141e+03	2.0733e+03	6.8571e+03	2.7151e+03	-1.8966e+03	1.9262e+04	758.3810	461.1567	2.8789e+03	-1.3674e+03	2.7005e+03	1.0803e+03	4.5071e+03	1.4821e+03	2.0868e+03	-2.4755e+04	1.9633e+03	1.0798e+04	-8.3304e+03	-846.1069	2.6364e+03	6.7798e+03	-6.6990e+03	377.6565	7.5746e+03	376.7528	2.2947e+03	-1.0017e+04	2.9205e+03	-7.8008e+03	-1.5008e+05	1.2221e+03	1.7682e+03	-1.4879e+03	84.9686	-3.1345e+03	6.1718e+03	1.8116e+03	4.0599e+03	2.0151e+03	1.0096e+04	2.0988e+03	823.0430	-1.1950e+04;];

%for i=1:1
%[X,FVAL,EXITFLAG,OUTPUT]=fminsearch('unifac_parameters_ObjFun_all_adsorbents',aa(i,:),optimset('MaxFunEvals',80000,'TolX',1e-30,'TolFun',1e-30,'MaxIter',100000));
[X,FVAL,EXITFLAG,OUTPUT]=fminsearch('unifac_parameters_ObjFun_all_adsorbents',a,optimset('MaxFunEvals',80000,'TolX',1e-30,'TolFun',1e-30,'MaxIter',100000));
%XX(i,:)=X;
%FVALL(i)=FVAL;
%end
X;








