% 1D Alpha_Dogs 
global n beta m
% n: Dimension, m:number of additional inequality constraints, meq=number of equality constraints 
n=1;  
% The upper and lower bound for variables
ss0=0; ss1=20;

% Scaling and changing the variables

% Initial calculations
bnd1=zeros(n,1); bnd2=ones(n,1); 
rk=ss1-ss0;

Ain=[]; bin=[];
Ain=[Ain; eye(n);-eye(n)];
bin=[bin; bnd2 ;-bnd1];

% Calculate the position of initial points
  xie=[0 1];
   
% Initial calculation ??????
for ii=1:size(xie,2)
    xmr=ss0+rk.*xie(:,ii);
   disp( sprintf( 'function evaluation') );
   disp( sprintf('Parameter: %e',xmr));
   prompt = 'What is Function value?\n';
   ymr = input(prompt); 
   prompt = 'What is uncertainty?\n';
   sigmamr = input(prompt);
end

% k_cr is the sigma multiplier
inter_method=0; beta=2; iter_max=100;
% Perform the iterations

for k=1:iter_max
    % Find the Delta-Dog point with regression
   [inter_par,yp]=regressionparametarization(xiT(:,n+2:end),yiT(n+2:end),sigmaT(:,n+2:end),inter_method);
   % Disceret search function
   sd=(min(yp,2*yiT(n+2:end)-yp)-y0)./sigmaT(n+2:end);
   [t,ind_min]=min(yiT); [t,ind_s]=min(sd);
   if ind_s~=ind_min
        % probably changed: Improve the measurement at xiT(:,ind_s)
        [yiT(:,ind_s), sigmaT(ind_s)]=improve(ind_s, transtime(ind_s));
        % Improvement 
   else
    [xie ,yie ,inde,tri_ind,tri] = neighberpoints_find(xiT,yiT,5);
    [y,ind]=min(yiT); xmin=xiT(:,ind); [xm ym]=inter_min(xmin,inter_par);
    deltam=mindis(xm,xiT);
    if (ym>y0) 
        [t,ind]=min(yiT);
        if glo==1
            [xm cse]= tringulation_search_bound(inter_par,xiT,tri);
        else
            [xm cse]=tringulation_search_local(inter_par,xiT,tri_ind ,tri);
        end
    end
   deltam=mindis(xm,xiT);
   if sigmaT(ind_min)<deltam*dxdt
    % scaled point
       xmr=ss0+rk.*(x0+V*xm);
   [yim, sigmam, transm]=initial_cal(xmr,size(xiT,2)+1); 
   %  Modify the data set
   sigmaT=[sigmaT sigmam]; yiT=[yiT yim]; transtime=[transtime transm];
   else
       xmr=ss0+rk.*(x0+V*xiT(:,ind_min));
   [yiT(:,ind_min), sigmaT(ind_min)]=improve(ind_min,transtime(ind_min));    
   end
   end
   % Find the most 
   figure(1)
    subplot(2,1,1)
    plot(1:length(yiT(n+2:end)),yiT(n+2:end),'-')
    ylim([0 2])
    subplot(2,1,2)
    plot(xiT(:,n+2:end)')
    drawnow
end