%Constraint Functions and Their Gradients Evaluation M-File for Example 12.3
% File name: ConstAndGrad3_5_3
function [c, ceq] = confun_parameter_all_adsorbents(at)
ab1=at(1);   ab2=at(2); 
ab3=at(3);   ab4=at(4);
ab5=at(5);   ab6=at(6); 
ab7=at(7);   ab8=at(8); 
ab9=at(9);   ab10=at(10);
ab11=at(11); ab12=at(12); 
ab13=at(13); ab14=at(14);
ab15=at(15); ab16=at(16);
ab17=at(17); ab18=at(18); 
ab19=at(19); ab20=at(20); 
ab21=at(21); ab22=at(22); 
ab23=at(23); ab24=at(24);
ab25=at(25); ab26=at(26); 
ab27=at(27); ab28=at(28);
ab29=at(29); ab30=at(30);

aa1=at(31);  aa2=at(32); 
aa3=at(33);  aa4=at(34);
aa5=at(35);  aa6=at(36);  
aa7=at(37);  aa8=at(38);
aa9=at(39);  aa10=at(40);
aa11=at(41); aa12=at(42);
aa13=at(43); aa14=at(44);

% Nonlinear inequality constraints

c=[ab1-100000;-50000-ab1;ab2-100000;-50000-ab2;ab3-100000;-50000-ab3;ab4-100000;-50000-ab4;ab5-100000;-50000-ab5;ab6-100000;-50000-ab6;ab7-100000;-50000-ab7;ab8-100000;-50000-ab8;ab9-100000;-50000-ab9;ab10-100000;-50000-ab10;ab11-100000;-50000-ab11;ab11-100000;-50000-ab11;ab12-100000;-50000-ab12;ab13-100000;-50000-ab13;ab14-100000;-50000-ab14;ab15-100000;-50000-ab15;ab16-100000;-50000-ab16;ab17-100000;-50000-ab17;ab18-100000;-50000-ab18;ab19-100000;-50000-ab19;ab20-100000;-50000-ab20;ab21-100000;-50000-ab21;ab22-100000;-50000-ab22;ab23-100000;-50000-ab23;ab24-100000;-50000-ab24;ab25-100000;-50000-ab25;ab26-100000;-50000-ab26;ab4-100000;-50000-ab26;ab27-100000;-50000-ab27;ab28-100000;-50000-ab28;ab29-100000;-50000-ab29;ab30-100000;-50000-ab30;aa1-100000;-50000-aa1;aa2-100000;-50000-aa2;aa3-100000;-50000-aa3;aa4-100000;-50000-aa4;aa5-100000;-50000-aa5;aa6-100000;-50000-aa6;aa7-100000;-50000-aa7;aa8-100000;-50000-aa8;aa9-100000;-50000-aa9;aa10-100000;-50000-aa10;aa11-100000;-50000-aa11;aa12-100000;-50000-aa12;aa13-100000;-50000-aa13;aa14-100000;-50000-aa14;];

% Nonlinear equality constraints
ceq = [];

end