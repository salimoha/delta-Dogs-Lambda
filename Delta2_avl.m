% Delta_Dogs for Simulation Based on AVL with 2 Bazier
% AddToPath and seten
clear all
close all
clc


cd FilesToPath
addTopath;
cd ..

setenvAVL


global n Ain bin y0   
% n: Dimension, min:number of additional inequality constraints 

n=7; 
%n=2;
%n=2;
%fun=@(x) log(1+sum(1./x(1))+sum(x(1))-2+(x(2)-0.1)^2); 
%fun=@(x) sum(1./x(1))+sum(x(1))-2+(x(2)-0.1)^2; 
%fun=@(x) log(1+10*(x(1)+0.1)+1/(10*(x(1)+0.1))*(x(2)-0.3)^2-2);
fun=@ plygen_n;
%fun= @(x) sum((x-[0.3;-0.2;-0.1;0.8;0.1;0;0.3]).^2);
%fun= @rastriginn;
%fun=@(x) log(1+(1.5-x(1)+x(1)*x(2))^2+(2.25-x(1)+x(1)*x(2)^2)^2+(2.625-x(1)+x(1)*x(2)^3)^2);
%fun=@(x)(x(1)-1)^2+10*(x(2)-x(1)^2)^2;
%fun=@(x) (sum(x.^4)-16*sum(x.^2)+5*sum(x)+39.16616*2*n);
% Surface,z1,y2,z2,w1,ctip,wc, 

%lb=[0.2;0.5;0.5;-0.3;1/11;0.05;1/11]; 
%ub=[0.5;1.5;1.5;0.3;1/4;0.5;1/1.5];
%lb=[0.2;0.5;0.5;-0.3;4;0.05;1.5]; 
%ub=[0.5;1.5;1.5;0.3;11;0.5;11];
lb=zeros(n,1); ub=ones(n,1);
%lb=[-4.5 -4.5]'; ub=[4.5 4.5]';
%lb=-5*ones(n,1); ub=5*ones(n,1);
%lb=[0.1 ;0]; ub=[10 ;1];
%lb=[-2;-2]; ub=[2;2];

% Upper bound and lower bound for f(x)
%ymax=-1; y00= -38;
%ymax=0.05;
%y00=0.0273;
ymax=100;
%y00=log(1+100*1/37);
%ymax=7;
%ymin=4;
y00=0.0267;
ss0=lb; rk=ub-lb;
Ain=[eye(n); -eye(n)]; bin=[ones(n,1); zeros(n,1)];
% Calculate the position of initial points
  [xie ,xic]= Linear_constrained_initilazation(Ain,bin);
   O=xie*ones(n+1,1)/(n+1); xie=[O xie]; xiT=[xic xie];
   [xiT ,xie]=inter_add_fix(xiT,xie);
% Initial calculation ??????
yiT=(1:size(xiT,2))*0; 
for ii=n+2:size(xiT,2)
    Nm=20;
    xiT(:,ii)=round(xiT(:,ii)*Nm)/Nm;
    xmr=ss0+rk.*xiT(:,ii); 
    yiT(ii)=min(ymax,fun(xmr));
end
yiT(1:n+1)=inf;
 piT=yiT;


% Perform the iterations
iter_max=200; inter_method=1;
delta_min=0.001;
%y00=0;
for k=1:iter_max
%y0=(min(yiT)+y00)/2;
y0=y00; 
% Find the regression
%y0=y00;
    inter_par=interpolateparametarization(xiT(:,n+2:end),yiT(n+2:end),inter_method);
 %  error=0*yiT;
%for ii=n+2:size(xiT,2)
%      p(ii)=interpolate_val(xiT(:,ii),inter_par);
%      error(ii)=max(min(yiT(n+2:end))-p(ii),0);
%end

   
 %  inter_par=regressionparametarization(xiT(:,n+2:end),yiT(n+2:end),(yiT(n+2:end)-min(yiT))/2,1);
% Disceret search function
%y0=y00;
       [y,ind_min]=min(yiT); xmin=xiT(:,ind_min); 
      [xm ym]=inter_min(xmin,inter_par);
      Nm=100;
      xm=round(xm*Nm)/Nm;
      tri=delaunayn(xiT.');
      if (ym>y0)
    %    [t,ind]=min(yiT);
        xm = tringulation_search_bound1(inter_par,xiT);
      else
          %delta=0.2;
          %[tt,ind,x0]=mindis(xm,xiT(:,n+2:end));
          %xm=xiT(:,ind);
        %  tt=xm-xmin; tt=min(tt,delta); tt=max(tt,-delta);
        %  xm=tt+xmin;
          %if mindis(xm,xiT(n+2:end))<0.05
          %    xm1=xm;
          %end
      end
     % xm= lin_convex_bnd_project(xm,xiT,tri);
      xm=round(xm*Nm)/Nm;
      [xm]=new_add_trust_region(xm,xiT,lb,ub,1);
      %delta=0.1;
          %[tt,ind,x0]=mindis(xm,xiT(:,n+2:end));
          %xm=xiT(:,ind);
          %tt=xm-x0; tt=min(tt,delta); tt=max(tt,-delta);
          %xm=tt+x0;
     % xm=round(xm*20)/20;
        if mindis(xm,xiT(:,n+2:end))<delta_min
         break
        end
        ym=fun(ss0+rk.*xm);
     yiT=[yiT ym]; xiT=[xiT xm];
     piT=[piT interpolate_val(xm,inter_par)];
    figure(1)
    subplot(2,1,1)
   % plot(1:length(yiT(n+2:end)),yiT(n+2:end),'-', 1:length(yiT(n+2:end)),piT(n+2:end))
    plot(1:length(yiT(n+2:end)),yiT(n+2:end),'-')
   ylim([0 0.1])
    subplot(2,1,2)
    plot(xiT(:,n+2:end)')
    drawnow        
        
        
        
end
figure(2)
plot(1:length(yiT)-n-1,yiT(n+2:end),'-','linewidth',2)
axis([0 350 -1 1000])
grid on
set(gca,'fontsize',30)
axes('Position',[0.7,0.3,0.28,0.28])
plot(295:length(yiT),yiT(295:end),'-','linewidth',1.5)
set(gca,'fontsize',15)
grid on