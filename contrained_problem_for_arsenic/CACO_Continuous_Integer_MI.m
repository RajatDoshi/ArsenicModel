 function dd = CACO_Continuous_Integer_MI(NCONV,NDISCV,K,SAI,q,rho, nAnt,Max_Iter, CONSTOL, TOL);
%  format long
%--------------------------------------------------------------------------
% This program is for solving an mixed integer optimizatin problem (MINLP)
% using an ant colony optimization method. This is first draft established 01-14-2014.
%               Author: Berhane H. Gebreslassie
%               E-mail: berhane@vri-custom.org
%               Vishwamitra Research Institute
% Center of Uncertain Systems: Tools for Optimization and Management
%
%--------------------------------------------------------------------------
% clear all
% clc

%--------The optimization problem from the test function-------------------


objfunc = 'PenaltyFunction';                                % objective function
%--------------------------- Parameters------------------------------------
NDIM = NCONV + NDISCV;    % Total number of decision variables continuous plus discrete variables
Best_OF = repmat(9999999,Max_Iter,1);
TSOL = repmat(99999*ones(1,NDIM),Max_Iter,1);
%-----------lower & uper bound of the solution space-----------------------

smin = 0; smax = 24;
% smin = -30; smax = 30;
for i = NCONV+1:NDIM
    lbound(i) = 0;
    ubound(i) = 24;
end
for i = 1:NCONV
    lbound(i) = 0;
    ubound(i) = 24;
end

% %-----------Initializing of the solution archive---------------------------
% for j= 1:K
%     for i= 1:NDIM
%         meanj(j,i)= lbound(i) + (2*j-1)*(ubound(i)-lbound(i))/(2*K);
%         sigmaj(j,i) = (ubound(i) - lbound(i))/(2*K);
%         SAI(j,i)=meanj(j,i)+sigmaj(j,i)*randn;
%     end
% end
% %---------Solution Archives Initializing-----------------------------------
TF = zeros(K,NDIM+2);                % solution archive
nASF = zeros(nAnt,NDIM+2);           % New ants solution archive
AOASF = [TF;nASF]   ;                % Ant solutions and old solution archive
Wt = zeros(K,1);                     % weight of each solution
%  SAI=SAIA ;
%---------------discrete values--------------------------------------------
% SACON     - Solution archive for continuous variable
% SADISC    - Solution archive for integer variables
% SAFOROF=zeros(K,NDIM);
SACON(:,1:NCONV) = SAI(:,1:NCONV) ;

if NDISCV > 0
    SADISC(:,1:NDISCV) = Rounding(SAI(:,NCONV+1:NDIM),smin,smax);
else
    SADISC(:,1:NDISCV)= double.empty;
end

SAFOROF(:,1:NDIM) = [SACON SADISC];      % initial solution configration
%--------------VECTORS-----------------------------------------------------
%   OF     -  Objective function value
%   OFAS   -  Objective functtion value of ant solutions
%   OFOAS  -  Objective functtion value of old archive and ant solutions
%   WT     -  Weighting of the objective function value in solution archive
%--------------------------------------------------------------------------
OF=zeros(K,1); OFAS = zeros(nAnt,1); OFOAS=zeros(K+nAnt,1);
%--------------objective function evaluation ------------------------------
%---Rounding the integer variable solution is only for evaluating the OF---
for k=1:K
    OF(k)=feval(objfunc,SAFOROF(k,1:NDIM));
end

%--Sorting solution archives according to the objective function value ----
%--------and determining the weighting according to their ranking----------
TT= [SAI OF];                   % The solution archive stores only the in the continuous value
T=sortrows(TT,NDIM+1);
WT=CALWeight(K,q);
TF= [T WT];
WTsum=sum(WT);                  % Weighting sum
for k=1:K
    prob(k)=WT(k)/WTsum;
end
% -----------Ant Colony Optimization for Mixed Variable Search Domain ---------
% ------------- until termination criteria is satisfied -------------------
%  -------------- in this case the number of iterations -------------------
%--------------------------------------------------------------------------
ASADISC = zeros(nAnt,NDISCV);
ASACONV = zeros(nAnt,NCONV);
nASn = zeros(nAnt,NDIM);
SAFOROF_ANT = zeros(nAnt,NDIM);
counter = 0 ;
iter = 0;
check=0;

while (iter <= Max_Iter & counter <= CONSTOL)
    iter = iter + 1;

    for ia=1:nAnt
          %-------------Solution construction -----------------------------------
    %--------Probablistically choose the guiding solution------------------
        gaussianj=rand;
    for jnew=1:K
        if gaussianj < sum(prob(1:jnew))
            newGuassindex=jnew ;    % Index of the selected gaussiani function
            break
        end % if
    end %for
    %----------------------------------------------------------------------
    %--------retrive mean and and calculate the standared deviation -------
    %-------------of the guiding solution for the next iteration-----------
    [mean Sigma] = stdev(T(1:K,1:NDIM),rho,newGuassindex); % 
        for id=1:NDIM
            nASn(ia,id) = mean(id)+Sigma(id)*randn;     % Sampling the choosen guasian distribution
        end %id
        ASACON(ia,1:NCONV) = nASn(ia,1:NCONV) ;
        if NDISCV > 0
            ASADISC(ia,1:NDISCV) = Rounding(nASn(ia,NCONV+1:NDIM),smin,smax);
        else
            ASADISC(ia,1:NDISCV) = [];
        end
    
    SAFOROF_ANT(ia,1:NDIM) = [ASACON(ia,1:NCONV) ASADISC(ia,1:NDISCV)];
    OFAS(ia) = feval(objfunc,SAFOROF_ANT(ia,1:NDIM));
    
    end
    %---Rounding the discrete variable solution is only for evaluating the OF---
    %         ASADISC(ia,1:NDISCV) = round(nASn(ia,NCONV+1:NDIM));
%     ASADISC(:,1:NDISCV)=rounding(nASn(:,NCONV+1:NDIM),smin,smax);
%     SAFOROF_ANT(:,1:NDIM) = [ASACON(:,1:NCONV) ASADISC(:,1:NDISCV)];
%     OFAS(1:nAnt) = feval(objfunc,SAFOROF_ANT(1:nAnt,1:NDIM));
    %ia
    %----------------------------------------------------------------------
    % the solution archive update only considers the continuous variables
    % and vertually continuous variables for the decrtete decision variables as follow
    %----------------------------------------------------------------------
    nASS = [nASn OFAS];                      % solution of each ant from the current iteration with continuous
    AOASS = [T; nASS];                       % combine the current solution archive and the ant solutions
    AOASS_sort = sortrows(AOASS,NDIM+1);      % Sort the combined solutions according to their objective functions
    TNEW = AOASS_sort(1:K,1:NDIM+1);          % Keep the first K solutions according to thier rank and remove solutions equivalent to the number of ant solutions
    
    %---------Update Solution Archive-------------------------------------
    [T, TJUNK] = UPDATE(T,TNEW);
    [Iter_Best_OF,ind] = min(nASS(:,NDIM+1)) ;  
    if NDISCV > 0
        soll = T(1,NCONV+1:NDIM);
        solll = Rounding(soll,smin,smax);
    else
        solll= [];
    end        
    TSOL_NEW = [T(1,1:NCONV) solll];        
    %---------Objective function updates-------------------------------------  
    if T(1,NDIM+1) < Best_OF
        Best_OF = T(1,NDIM+1);
        TSOL = TSOL_NEW;
    end % if
    
    CPUTIME=cputime;
    dd(iter,:) = [iter Best_OF Iter_Best_OF TSOL];
%     [solute_distribution selectivity  solvent_loss    res    T_bp  Mw_S ineqcons1 ineqcons2 ineqcons3 ineqcons4 eqcons1] = UNIFAC(dd(iter,4:NDIM+3)) ;
%     
%     DY(iter,:) = [solute_distribution selectivity solvent_loss  res  T_bp  Mw_S ineqcons1  ineqcons2  ineqcons3  ineqcons4 eqcons1];
    disp(['Iter ' num2str(iter) ': Best So-far Cost = ' num2str(dd(iter,2)) ': Iter Best = ' num2str(dd(iter,3)) ': Iter Var = ' num2str(dd(iter,4:NDIM+3))]); 
%    disp(['Iter ' num2str(iter) ': Best So-far Cost = ' num2str(dd(iter,2)) ': Iter Best = ' num2str(dd(iter,3)) ': Iter Var = ' num2str(dd(iter,4:NDIM+3))]);
%     
%     
    %---------Termination criteria counter---------------------------------
    % S     - New solution ( T(1,NDIM+1) is the objective function value )
    % SOLD  - Old solution ( TJUNK(1,NDIM+1) is the objective function value of the old solution )
    %----------------------------------------------------------------------
    OBF_NEW = T(1,NDIM+1);
    OBF_OLD = TJUNK(1,NDIM+1);

    if (abs(OBF_OLD-OBF_NEW) <= TOL)
        counter = counter + 1;
        if counter == 1
            check = iter;
        end
        if counter >= 2
            if (abs(check-iter) > (counter-1) | abs(check-iter) < (counter-1))
                counter = 0;
            end
        end
    else
        counter = 0;
    end
       
%     Best_OF;
end %iter
%  xlswrite('EACO.xls', [dd(:,:) DY(:,:)],1,'A4');

% figure(1);
% plot([1:iter],dd(:,2),[1:iter],dd(:,3));
subplot(2,2,1);
plot([1:iter],dd(:,2),[1:iter],dd(:,3));
legend('Best PF','iter Best PF');
title('Sofar best and Iteration best');

% subplot(2,2,2);
% plot([1:iter],dd(:,2),[1:iter],DY(:,1));
% legend('Best PF','Best OBF');
% title('Penalty function and Objective function');
% 
% subplot(2,2,3);
% plot([1:iter],DY(:,1),[1:iter],DY(:,2),[1:iter],DY(:,3));
% legend('m','beta','SL');
% title('Proporties');
% 
% subplot(2,2,4);
% plot([1:iter],DY(:,4),[1:iter],DY(:,5),[1:iter],DY(:,6),[1:iter],DY(:,7),[1:iter],DY(:,8),[1:iter],DY(:,9));
% legend('res','ineqcns1','ineqcns2','ineqcns3','ineqcns4','eqcns1');
% title('residual and constraints');

 end