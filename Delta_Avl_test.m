% Delta_Dogs for Simulation Based on AVL
% AddToPath and seten
clear all
close all
clc


cd FilesToPath
addTopath;
cd ..

setenvAVL


global n beta m Ain bin y0   
% n: Dimension, min:number of additional inequality constraints 

n=9;  
fun=@plygen;


% First 3 parameters: theta_0, theta_1, theta_2.
A=[1 1 1]; B=pi/2;
lb=[0.01 pi/3 0.01]'; ub=[pi/2 pi/2 pi/2]'; 
[A,B,s,r]= actual_scaling_lincon(A,B,ub,lb,[],[],3,0);

Ain=[A zeros(size(A,1),n-3)]; bin=B;
ss0=s; rk=r;
% Second 3 parameters: l_0, l_1, l_2.
A=[1 1 1]; B=2;
lb=[0.5 0.1 0.7]'; ub=[1.5 0.5 1.5]'; 
[A,B,s,r]= actual_scaling_lincon(A,B,ub,lb,[],[],3,0);

Ain=[Ain; [zeros(size(A,1),3) A zeros(size(A,1),n-6)]]; bin=[bin ;B];
ss0=[ss0 ;s]; rk=[rk ;r];

% Last four parameters:
%lb=[0.2 0 0 0]'; ub=[0.5 1 1 1]';
lb=[0.2 0 0]'; ub=[0.5 1 1]';
%[A,B,s,r]= actual_scaling_lincon([],[],ub,lb,[],[],,0);
[A,B,s,r]= actual_scaling_lincon([],[],ub,lb,[],[],length(lb),0);

Ain=[Ain; [zeros(size(A,1),6) A]]; bin=[bin ;B];
ss0=[ss0 ;s]; rk=[rk ;r];


% Calculate the position of initial points
  [xie ,xic]= Linear_constrained_initilazation(Ain,bin);
   O=xie*ones(n+1,1)/(n+1); xie=[O xie]; xiT=[xic xie];
   [xiT ,xie]=inter_add_fix(xiT,xie);
   
% Initial calculation ??????
yiT=(1:size(xiT,2))*0; 
for ii=n+2:size(xiT,2)
    xmr=ss0+rk.*xiT(:,ii); yiT(ii)=min(-10,fun(xmr));
end
%yiT(1:n+1)=Inf; sigmaT(1:n+1)=Inf; transtime(1:n+1)=Inf;

% Perform the iterations
iter_max=0; inter_method=1; y0=-35;
for k=1:200

% Find the regression
   inter_par=interpolateparametarization(xiT(:,n+2:end),yiT(n+2:end),inter_method);
   
% Disceret search function
       [y,ind_min]=min(yiT); xmin=xiT(:,ind_min); 
       [xm ym]=inter_min(xmin,inter_par); 
    if (ym>y0) 
        [t,ind]=min(yiT);
        tri=delaunayn(xiT.');
        [xm cse]= tringulation_search_bound(inter_par,xiT,tri,ind_min);
        xm= lin_convex_bnd_project(xm,xiT,tri);
    end
      % deltam=mindis(xm,xiT);
     if mindis(xm,xiT(n+2:end))<0.05
         break
     end
     ym=min(-10,fun(ss0+rk.*xm));
     yiT=[yiT ym]; xiT=[xiT xm];
   
   % plotting 
   figure(1)
    subplot(2,1,1)
    plot(1:length(yiT(n+2:end)),yiT(n+2:end),'-')
    ylim([-45 0])
    subplot(2,1,2)
    plot(xiT(:,n+2:end)')
    drawnow
    
  
end