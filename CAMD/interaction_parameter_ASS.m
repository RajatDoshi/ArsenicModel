%Activity coefficient parameters adsorbent solid solution
function Anm=interaction_parameter_ASS(Gyes)
col=length(Gyes);
Anm=zeros(col,col);

          %As           %SiO2    %Al2O3      %Na2O       %MgO    % Fe2O3       %K2O      %CaO        %Na+   
anmASSi=[  0          -150080.00 1768.20     1222.10 	-3134.5    6171.8	  84.9686	  -1487.9   6857.1;    
           1811.60    0          0          0           0          0           0          0         0;
           2015.10    0          0          0           0          0           0          0         0;
           4059.90    0          0          0           0          0           0          0         0;
           823.043    0          0          0           0          0           0          0         0;
           -11950     0          0          0           0          0           0          0         0;
           2098.8     0          0          0           0          0           0          0         0;
           10096      0          0          0           0          0           0          0         0;
           1080.3     0          0          0           0          0           0          0         0;];
       
 for kk=1:col
	if Gyes(kk) ==1	
       Anm(:,kk)=anmASSi(:,kk);
       Anm(kk,:)=anmASSi(kk,:);
    end
     
    if Gyes(kk) ==0
       Anm(:,kk)=zeros(1,col);
       Anm(kk,:)=zeros(col,1);
    end
 end