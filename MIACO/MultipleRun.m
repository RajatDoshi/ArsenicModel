clear all
clc
NCONV   = 2 ;               % number of continuous variable
NDISCV  = 8;                % number of discrete variables
rho     = 0.75;             % EVAPORATION FACTOR 
q       = 0.0001;           % ALGORITHM PARAMETER
nAnt    = 50;
K       = 500;
NDIM = NCONV + NDISCV ;
running = 10;       
CONSTOL = 150;               % CONSEQUETIVES ITERATIONS

Max_Iter = 14000;            % Maximum iterations
TOL = 1E-7;                   % Tolerernce, the objective function value difference

%%addpath('C:\Users\Student\Documents\POST-DOC_VRI\NATURAL GAS\Adsorption of As in clays\UNIFAC_CLAYS\EACO -new version desicion variables10\EACO -10 - Case B_all groups for radium\CAMD')
addpath('C:/rajat_doshi/Diwekar/part_two_code/EACO_for_Arsenic/CAMD')

% [Sheet Range]=Resultexcel(prob,func,dim, sample)
  Sample_HSS = 1 ; 
  Sample_MCS = 2;
  func       = 2 ;                         % input('Enter the function number PR_CV & MV = 1, EL_CV  & MV =2, CG_CV  & MV = 3, EX3_MV =4  = ');
  fdomain    = 1;                          % input('Enter the domain type continous = 1 or mixed variable = 2  = ');
  
%  [Sheet Range] = Resultexcel(fdomain,func,NDIM, Sample_HSS);
 
 
 % %-----------Initializing of the solution archive---------------------------
 
%  for i = NCONV+1:NDIM
%     lbound(i) = 0;
%     ubound(i) = 24;
% end
% for i = 1:NCONV
%     lbound(i) = 0;
%     ubound(i) = 24;
%                                                                          
% end
% for j= 1:K
%     for i= 1:NDIM
%         meanj(j,i)= lbound(i) + (2*j-1)*(ubound(i)-lbound(i))/(2*K);
%         sigmaj(j,i) = (ubound(i) - lbound(i))/(2*K);
%         SAI(j,i)=meanj(j,i)+sigmaj(j,i)*randn;
%     end
% end

% Monte carlo sampling;
% MCS = fopen('Sample/HSSAdsRa500_20');
% [SAI, siz] = fscanf(MCS, '%g %g', [NDIM K]);
% fclose(MCS);
% SAI = SAI';

%Hammersley Sequence Sampling
HSS = fopen('Sample/HSSAdsRa500_20.txt');
[SAI, siz] = fscanf(HSS, '%g %g', [NDIM K]);
fclose(HSS);
SAI = SAI';

solution=zeros(running,NDIM+3);
solution2=zeros(running,11);
for m=1:running 
    iseed = m ;
    rng(iseed);
%   nAnt    = m*5;
 [bb candidate]= MIACO_ORD(NCONV,NDISCV,K,SAI,q,rho, nAnt,Max_Iter, CONSTOL, TOL);
  [M N] = size(bb); 
%   [MM NN] = size(dy); 
  ff(m) = M;
  best(m) = bb(M,2);
  solution(m,1:NDIM+6) = bb(M,1:NDIM +6);
%   solution2(m,1:NN)= dy(M,1:NN);
end

 ftr = [num2str(bb)];           % all iterations are printed 
 fcandidate = num2str(candidate);
%   ftr = [num2str(bb) num2str(dy)]; % all iterations are printed 
%--------------------------------------------------------------------------
% [solute_distribution selectivity  solvent_loss    res    T_bp  Mw_S ineqcons1 ineqcons2 ineqcons3 ineqcons4 eqcons1]= UNIFAC(dd(iter,4:NDIM+3)) ;

 xlswrite('EACO_CAMD.xls', [solution],1,'A5')
%  xlswrite('CACO.xls', [solution])
 Averare_iter = sum(ff)/length(ff);
 Min_iter = min(ff);
 Max_iter = max(ff);
  
 Min_OF = min(best);
 Max_OF = max(best);
 Averare_OF = sum(best)/length(ff);
 
 % disp(['Iter ' num2str(iter) ': Best So-far Cost = ' num2str(dd(iter,2)) ': Iter Best = ' num2str(dd(iter,3)) ': Iter Var = ' num2str(dd(iter,4:NDIM+3))]);
    
 
 
%  disp(['-----------------------------------------------------------------------------------------------------------------------------------' ]);
%  disp(['Function         = objFun_amount_adsorbed_RaBa']);
%  disp(['sampling         = MCS ']);
%  disp(['Dimension        = 200 ']);
%  disp(['Archive Size     = 200 ']);
%  disp(['nAnts            = 15 ']);
%  disp(['CON Iter Counter = 40 ']);
%  disp(['TOL              = 1E-6 ']);
%  disp(['q                = 1E-5 ']);
%  disp(['lower bound      = -30']);
%  disp(['upper bound      = 30 ']);
 
 
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