clear all
clc
NCONV   = 44;
NDISCV  = 0;
rho     = 0.85;
q       = 1E-5;
nAnt    = 15;
% q     = 0.5;
K       = 300;
NDIM = NCONV + NDISCV ;
running = 10;
CONSTOL = 3000;

Max_Iter = 50;            % Maximum iteration
TOL =1E-6;                   


% [Sheet Range]=Resultexcel(prob,func,dim, sample)
  Sample_HSS = 1 ; 
  Sample_MCS = 2;
  func       = 2 ;                         % input('Enter the function number PR_CV & MV = 1, EL_CV  & MV =2, CG_CV  & MV = 3, EX3_MV =4  = ');
  fdomain    = 1;                          % input('Enter the domain type continous = 1 or mixed variable = 2  = ');
  
%  [Sheet Range] = Resultexcel(fdomain,func,NDIM, Sample_HSS);


MCS = fopen('Samples/MCS44.txt');
[SAI, siz] = fscanf(MCS, '%g %g', [NDIM K]);
fclose(MCS);
SAI = SAI';

solution=zeros(running,NDIM+3);
solution2=zeros(running,11);
for m=1:running 
%   nAnt    = m*5;
 bb= CACO_Continuous_Integer_MI(NCONV,NDISCV,K,SAI,q,rho, nAnt,Max_Iter, CONSTOL, TOL);
  [M N] = size(bb); 
%   [MM NN] = size(dy); 
  ff(m) = M;
  best(m) = bb(M,2);
  solution(m,1:NDIM+3) = bb(M,1:NDIM +3);
%   solution2(m,1:NN)= dy(M,1:NN);
end

 ftr = [num2str(bb)];
%--------------------------------------------------------------------------
% [solute_distribution selectivity  solvent_loss    res    T_bp  Mw_S ineqcons1 ineqcons2 ineqcons3 ineqcons4 eqcons1]= UNIFAC(dd(iter,4:NDIM+3)) ;

%  xlswrite('EACO_RUN.xls', [solution],1,'B18')
%  xlswrite('CACO.xls', [solution])
 Averare_iter = sum(ff)/length(ff);
 Min_iter = min(ff);
 Max_iter = max(ff);
  
 Min_OF = min(best);
 Max_OF = max(best);
 Averare_OF = sum(best)/length(ff);
 
 % disp(['Iter ' num2str(iter) ': Best So-far Cost = ' num2str(dd(iter,2)) ': Iter Best = ' num2str(dd(iter,3)) ': Iter Var = ' num2str(dd(iter,4:NDIM+3))]);
    
 
 
%  disp(['-----------------------------------------------------------------------------------------------------------------------------------' ]);
%  disp(['Function         = Solvent Selection']);
%  disp(['sampling         = MCS ']);
%  disp(['Dimension        = 100 ']);
%  disp(['Archive Size     = 500 ']);
%  disp(['nAnts            = 10 ']);
%  disp(['CON Iter Counter = 40 ']);
%  disp(['TOL              = 1E-6 ']);
%  disp(['q                = 1E-5 ']);
%  disp(['lower bound      = -30']);
%  disp(['upper bound      = 30 ']);
%  
%  
disp(['---------------------------------------------------------------------------------------------------------------------------------------' ]);
disp(['it1 ', ' it2 ', ' it3 ', ' it4', ' it5 ',' it6 ','it7 ', ' it8 ',' it9 ',' it10 ']);
disp([num2str(ff)]);
disp(['Aveg_it ']);
disp([num2str(Averare_iter)]);

disp(['---------------------------------------------------------------------------------------------------------------------------------------' ]);
disp(['     OF1  ','     OF2   ','       OF3    ','     OF4    ','      OF5    ','     OF6    ','     OF7    ','     OF8    ','    OF9    ','    OF10   ', ]);
disp([num2str(best)]);

disp(['Aveg_OF ']);
disp([num2str(Averare_OF)]);
disp(['------------------------------------------------------------------------------------------------------------------------------------------' ]);