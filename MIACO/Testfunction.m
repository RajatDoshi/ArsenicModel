function f=Testfunction(x)
funnum = 2;
if funnum == 1%F1 Parabolic Function
    sum=0;
    for i=1:length(x)
        sum=sum+x(i).^2;
    end
    f=sum;
elseif funnum == 2 %F Ellipsoid Function
    %     [K NDIM]=size(x);
    NDIM=length(x);
    sum=0;
    for i = 1:NDIM
        sum = sum +(5^((i-1)/(NDIM-1))*x(i))^2;
    end
    f=sum;
elseif funnum == 3 %F Cigar Function
    sum=0;
    for i = 2:length(x)
        sum = sum + x(i)^2;
    end
    f=10^4*sum + x(1)^2 ;
elseif funnum == 4 %F3 Prof. Diwaker (Example) Function
    suma=0;
    NDIM=length(x);
    NCONV = 7;
    for i=1:NCONV
        suma=suma+abs(x(:,i)-i/NCONV).^2;
    end
    sumb=0;
    for i=NCONV+1:NDIM
        sumb=sumb+x(:,i).^2;
    end
    prodc=1.0;
    for i=NCONV+1:NDIM
        prodc=prodc*cos(4*pi*x(:,i));
    end
    f=suma+(sumb-prodc);
    
elseif funnum == 5 %F Rosenbrock Function
    n=length(x);
    sum=0;
    for i = 1:n-1
        sum = sum + 100*(x(i)^2-x(i+1))^2 + (x(i)-1)^2;
    end
    f=sum;
    
    
elseif funnum == 6 %F2
    summ=0;
    f=summ;
    for i=1:round(x(:,1))
        summ=summ+(x(:,1)-3).^2+ (x(:,2)-3).^2 +(x(:,3)-3).^2 ;
    end
    f=summ;
    
elseif funnum == 7 %F5
    %     [K NDIM]=size(x);
    NDIM=length(x);
    sum=0;
    for i = 1:NDIM
        sum = sum +(100^((i-1)/(NDIM-1))*x(i))^2;
    end
    f=sum;
    
elseif funnum == 8 %F6
    f=abs(x)+sin(x);
elseif funnum == 9 %F6
    f=x(:,11)^2-6*(x(:,11)+x(:,1)+x(:,2)+x(:,3)+x(:,4)+x(:,5)+x(:,6)+x(:,7)+x(:,8)+x(9)+x(:,10))+x(:,1)^2+x(:,2)^2+x(:,3)^2+x(:,4)^2+x(:,5)^2+x(:,6)^2+x(:,7)^2+x(:,8)^2+x(:,9)^2+x(:,10)^2+99;
      
    
    %----------------------Mixed Variable functions
    
elseif funnum == 10 %F Ellipsoid
    %     [K NDIM]=size(x);
    n=length(x);
    suma = 0;
    sumb = 0;
    for i = 1 : n/2
        suma = suma +(5^((i-1)/(n-1))*x(i))^2;
    end
    for i = n/2+1 : n
        sumb = sumb +(5^((i-1)/(n-1))*x(i))^2;
    end
    f=suma +sumb;
elseif funnum == 11 %F Rosenbrock
    n=length(x);
    suma = 0;
    %     sumb = 0;
    for i = 1 : n-1
        suma = suma + 100*(x(i)^2-x(i+1))^2 + (x(i)-1)^2;
    end
    %     for i = n/2+1 : n-1
    %         sumb = sumb + 100*(x(i)^2-x(i+1))^2 + (x(i)-1)^2;
    %     end
    %     f=suma +sumb;
    f=suma;
elseif funnum == 12 %F Cigar
    n=length(x);
    suma = 0;
    sumb = 0;
    for i = 2:n/2
        suma = suma + x(i)^2;
    end
    sumc=10^4*suma + x(1)^2 ;
    for i = (n/2 + 1) : n
        sumb = sumb + x(i)^2;
    end
    sumd=10^4*sumb;
    f = sumc + sumd ;
    
elseif funnum == 13 %F Cigar
    f= (x(1)-8.5)^2 + (x(2)-12.5)^2 + (x(3)-30)^2 + (x(4)-15)^2;
end
