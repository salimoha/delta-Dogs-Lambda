clear all
close all 
clc

load('1d-DG-data');
clearvars -except upTime xx

global n y0 Ain bin m kappa

 
% Dimension of problem
n=1;

% Get the upper and lower bound
ss0=0; ss1=20;

% Estimated value for the lower bound.
y0=0; delta_tol=0.01; k_max=1; 

% Calculations
kappa=1.2; m=2*n; inter_method=1;
rk=ss1-ss0;


% Initialization
Ain=[eye(n);-eye(n)]; bin=[ones(n,1);zeros(n,1)];
xie=[0 1];
yie=[];
ff=1;
% 
for ii=1:length(xie)
    xmr=ss0+rk.*xie(:,ii);
    [rt,ind]=min(abs(upTime-xmr));
     yie=[yie xx(ind)];
end

for k=1:11
 inter_par= interpolateparametarization(xie,yie,inter_method);
 [y,ind]=min(yie); xmin=xie(:,ind);  
 tri=delaunayn(xie.');
 [xm ym]=inter_min(xmin,inter_par);
 if ym>y0
 [xm cse]= tringulation_search_bound(inter_par,xie,tri);
 end
   delta=mindis(xm,xie);
   if delta<delta_tol
       break
   end
   xie=[xie xm];
 xmr=ss0+rk.*xm;
 [rt,ind]=min(abs(upTime-xmr));
 yie=[yie xx(ind)];
 
 plot(upTime,xx,'-',ss0+rk.*xie,yie,'*')
end





