function [Sheet Range]=Resultexcel(prob,func,dim, sample);
%--------------------------------------------------------------------------
% This function is helps to identify the the sheet and the range in excele 
% file where the results will be transported upon solving the problem. 
% sheet  -  The sheet number where the results has to be transported
% Range  -  Range of the results in the excel file 
% prob - either continuous or discrete function cont = 1 & disc = 2
%        - Continuous variable      = 1
%        - Mixed variable           = 2
% func - the function optimized using the CACO 
%        - Parabolic function       = 1
%        - Ellipsoid  function      = 2
%        - Cigar function           = 3
%        - Diwakar et al. (example) = 4
% dim - the dimension of function optimized using the CACO 
% sample - the sampling method used for intializing the solution archive 
%        - Hammerslay Sampling(HSS) = 1
%        - Monte Carlo Sampling(MCS)= 2


%--------------------SHEET SELECTION---------------------------------------

%--------------------SHEET FOR CONTINUOUS FUNCTIONS------------------------

if (prob == 1 & dim == 5)
    Sheet = 1;
end

if (prob == 1 & dim == 10)
    Sheet = 2;
end

if (prob == 1 & dim == 15)
    Sheet = 3;
end

if (prob == 1 & dim == 20)
    Sheet = 4;
end
%--------------------SHEET FOR MIXED VARIABLE FUNCTIONS--------------------
if (prob == 2 & dim == 5)
    Sheet = 5;
end

if (prob == 2 & dim == 10)
    Sheet = 6;
end

if (prob == 2 & dim == 15)
    Sheet = 7;
end

if (prob == 2 & dim == 20)
    Sheet = 8;
end
%--------------------RANGE SELECTION-----DIM = 5---------------------------
if (func == 1 & dim == 5 & sample == 1)
    Range= 'B19';
end
if (func == 1 & dim == 5 & sample == 2)
    Range= 'M19';
end

if (func == 2 & dim == 5 & sample == 1)
    Range= 'B43';
end

if (func == 2 & dim == 5 & sample == 2)
    Range= 'M43';
end
if (func == 3 & dim == 5 & sample == 1)
    Range= 'B67';
end

if (func == 3 & dim == 5 & sample == 2)
    Range= 'M67';
end

if (func == 4 & dim == 5 & sample == 1)
    Range= 'B91';
end
if (func == 4 & dim == 5 & sample == 2)
    Range= 'M91';
end

if (func == 5 & dim == 5 & sample == 1)
    Range= 'B116';
end
if (func == 5 & dim == 5 & sample == 2)
    Range= 'M116';
end

%--------------------RANGE SELECTION-----DIM = 10--------------------------
if (func == 1 & dim == 10 & sample == 1)
    Range= 'B19';
end
if (func == 1 & dim == 10 & sample == 2)
    Range= 'Q19';
end

if (func == 2 & dim == 10 & sample == 1)
    Range= 'B43';
end

if (func == 2 & dim == 10 & sample == 2)
    Range= 'Q43';
end

if (func == 3 & dim == 10 & sample == 1)
    Range= 'B67';
end

if (func == 3 & dim == 10 & sample == 2)
    Range= 'Q67';
end

if (func == 4 & dim == 10 & sample == 1)
    Range= 'B91';
end
if (func == 4 & dim == 10 & sample == 2)
    Range= 'Q91';
end

if (func == 5 & dim == 10 & sample == 1)
    Range= 'B116';
end
if (func == 5 & dim == 10 & sample == 2)
    Range= 'Q116';
end
%--------------------RANGE SELECTION-----DIM = 15--------------------------
if (func == 1 & dim == 15 & sample == 1)
    Range= 'B19';
end
if (func == 1 & dim == 15 & sample == 2)
    Range= 'W19';
end

if (func == 2 & dim == 15 & sample == 1)
    Range= 'B43';
end

if (func == 2 & dim == 15 & sample == 2)
    Range= 'W43';
end

if (func == 3 & dim == 15 & sample == 1)
    Range= 'B67';
end

if (func == 3 & dim == 15 & sample == 2)
    Range= 'W67';
end

if (func == 4 & dim == 15 & sample == 1)
    Range= 'B91';
end
if (func == 4 & dim == 15 & sample == 2)
    Range= 'W91';
end

if (func == 5 & dim == 15 & sample == 1)
    Range= 'B116';
end
if (func == 5 & dim == 15 & sample == 2)
    Range= 'W116';
end
%--------------------RANGE SELECTION-----DIM = 20--------------------------

if (func == 1 & dim == 20 & sample == 1)
    Range= 'B19';
end
if (func == 1 & dim == 20 & sample == 2)
    Range= 'AA19';
end
if (func == 2 & dim == 20 & sample == 1)
    Range= 'B43';
end

if (func == 2 & dim == 20 & sample == 2)
    Range= 'AA43';
end

if (func == 3 & dim == 20 & sample == 1)
    Range= 'B67';
end

if (func == 3 & dim == 20 & sample == 2)
    Range= 'AA67';
end

if (func == 4 & dim == 20 & sample == 1)
    Range= 'B91';
end
if (func == 4 & dim == 20 & sample == 2)
    Range= 'AA91';
end

if (func == 5 & dim == 20 & sample == 1)
    Range= 'B116';
end
if (func == 5 & dim == 20 & sample == 2)
    Range= 'AA116';
end



end


