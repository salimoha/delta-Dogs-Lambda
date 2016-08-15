% Experimenal based optimization

% Give the input for each variable as the upper and lower bound
% ss0(ii), xx1(ii) are the lower and upper bound for ii-the variable

clear all
close all
clc

global n y0 Ain bin m kappa

% Dimension of problem
n=2;

% test function
%fun=@rastriginn;

% Get the upper and lower bound
ss0=[-2 -2].';
ss1=[pi pi].';

% Estimated value for the lower bound.
y00=1; delta_tol=0.05; k_max=20; 

% Calculations
kappa=1.2; m=2*n; inter_method=1;
rk=ss1-ss0;


% Initialization
Ain=[eye(n);-eye(n)]; bin=[ones(n,1);zeros(n,1)];
[xie xic]= Linear_constrained_initilazation(zeros(n,1));
 ff=1;
O=xie*ones(n+1,1)/(n+1);
   xie=[O xie];
   xiT=[xic xie];
   [xiT xie]=inter_add_fix(xiT,xie);
   
   yiT(1:n+1)=inf*ones(n+1,1);
% 
for ii=n+2:size(xiT,2)
    %xiT(:,ii)=linreflection(xiT(:,ii),Np,delta_tol);
    xmr=ss0+rk.*xiT(:,ii);
   
    % Calculate the function value at xmr by User
 disp( sprintf( '%d function evaluation', ff ) );
 ff=ff+1;
       for sss=1:n
 disp( sprintf( 'Parameter %d: %e', sss,xmr(sss)) );
       end
 prompt = 'What is Function value?\n';
 ymr = input(prompt); 
 clc
 % end of User Calculation 
 
% ymr=fun(xmr);
   yiT=[yiT ymr];
end

for k=1:k_max
% Interior calculation to estimate next point for 
 inter_par= interpolateparametarization(xiT(:,n+2:end),yiT(n+2:end),inter_method);
 [y,ind]=min(yiT); xmin=xiT(:,ind);  
 tri=delaunayn(xiT.');
 [xm ym]=inter_min(xmin,inter_par);
 deltam=mindis(xm,xiT);
 y0=(y00+min(yiT))/2; 
 [xm cse]= tringulation_search_bound(inter_par,xiT,tri);
   xm= lin_convex_bnd_project(xm,xiT,tri);

   delta=mindis(xm,xiT);
   if delta<delta_tol
       break
   end
   xiT=[xiT xm];
 for ii=1:n
     if xm(ii)<0.05 xm(ii)=0; end
     if xm(ii)>0.95 xm(ii)=1; end
 end
 xmr=ss0+rk.*xm; % this should be disceretize 
 % perform experiment at xmr
% Calculate the function value at xmr by User
 disp( sprintf( '%d function evaluation', ff ) );
 ff=ff+1;
       for ii=1:n
 disp( sprintf( 'Parameter %d: %e', ii,xmr(ii)) );
       end
prompt = 'What is Function Value? \n';
ymr = input(prompt);
clc
 % end of User Calculation 

ymr=fun(xmr);
yiT=[yiT ymr];

 
end

