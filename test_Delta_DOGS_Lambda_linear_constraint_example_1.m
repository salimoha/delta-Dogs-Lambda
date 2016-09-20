% Second Order Delta-Dogs:
% Algorithm with second order convergence property.
% Constraints are simple bound constraints.
% Dimenson:
clear all; close all; clc
%   run('./FilesToPath/addTopath.m')
figure(1);clf;
qtype=1;
global n Ain bin y0 lattice  scale
n=5;
DELTA_MAX=1e-4*sqrt(n);

lattice = 'Zn '; [~,B,plane]=init_DOGS(n,lattice);
% xs = ones(n,1)*(-2.903);
% objective function:
% scaled to [0,1] domain
% xs =  (5+xs)/10;
xs = ones(n,1)*0.2097; 
%0<x<1 
fun=@(x)(sum((-5+10*x).^4)-16*sum((-5+10*x).^2)+5*sum(-5+10*x)+39.16616*2*n);
 %-5<x<5
% xs = ones(n,1)*(-2.903);
% fun=@(x) (sum((x).^4)-16*sum((x).^2)+5*sum(x)+39.16616*2*n);
% Constraints:
% upper nad lower bounds
% lob=-5*ones(n,1); upb=5*ones(n,1);
lob=zeros(n,1); upb=ones(n,1);
% Estimate of the solution:
y00=0; %y00 = -1e-5 ; % y00 = 0.01;
% Initial Calculation:
bnd1 = lob; bnd2 = upb;
ub = bnd2; lb = bnd1;
m=2*n;
% maximum number of iterations:
iter_max=100;
% interpolation strategy:
inter_method=8;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Input the equality constraints:?????????
% % Ain=[eye(n);-eye(n)]; bin=[bnd2 ;-bnd1];
% % Ain=[eye(n); -ones(1,n)/sqrt(n)]; bin=[5*ones(n,1); 2.8*n/sqrt(n)];
% % A=[ -ones(1,n)/sqrt(n)]; b=[ 2.8*n/sqrt(n)];
% A=([ -ones(1,n)/sqrt(n)]); b=[ 2.*n/sqrt(n)];
% Ain =[A; eye(n); -eye(n)]; 
% bin=[b ; ub; -lb];
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bs=0.7903; % active on constraint
bs = 0.81; % interior
A=[ones(1,n);-ones(1,n)]/sqrt(n); b=[ bs.*n/sqrt(n)];
Ain=[eye(n);-eye(n);A]
bin=[ones(n,1);zeros(n,1);b;b]
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% length_scale=10;
length_scale=1;
%length_scale=max(bnd2-bnd1);
% Calculate the Initial trinagulation points
% xU=bounds(bnd1,bnd2, n);
% xU=vertex_find(Ain,bin,[],[]);
xU=vertex_find(Ain,bin,lb,ub);
% Calculate initial delta
xE=(bnd1+0.5*(bnd2-bnd1));
% xE=(sum(xU')/(n+1))';
%xE=[0.350 ;0.4000 ;0.5000 ;0.5000 ;0.3000 ;0.3500; 0.3000];
% initial points: the midpoint point and its neighber
delta0=0.15*length_scale;
for ii=1:n
    e=zeros(n,1); e(ii)=1;
    xE(:,ii+1)=xE(:,1)+delta0*e;
end
% Calculate initial unevaluated delta
for ii=1:size(xU,2)
    deltaU(ii)=mindis(xU(:,ii),xE);
end
% Calculate the function at initial points
for ii=1:size(xE,2)
    yE(ii)=fun(xE(:,ii));
end
pE=yE;
interpolate_index=ones(1,length(yE));
%Nm=800;
Nm=32;
% Nm=100;
%y00=0;
%y00=0.02;
  y0=y00;
% Initialize constaints 
%tol=0.05;
rho=1;
for kk=1:8
for k=1:iter_max
%     keyboard
tic;
    stop=0;
    inter_par= interpolateparametarization(xE,yE,inter_method,interpolate_index);
    %y0=10;
    y0=y00;
    %y0=min(yE)-range(yE)/Nm;
    yup=[]; 
    yrp=[];
   % Calculate the metric of the unevaluated points.
            for ii=1:size(xU,2)
               yup(ii)=(interpolate_val(xU(:,ii),inter_par)-y0)/mindis(xU(:,ii),xE);
               yrp(ii)=interpolate_val(xU(:,ii),inter_par);
            end
          %  keyboard
   % check the minimum value of xU
            if size(xU,2)>0
                 [tup,ind_up]=min(yup);
            else
                 tup=inf;
            end
  if tup<0
      [trp,ind_rp]=min(yrp);
      [tt,ind]=min(yE); x0=xE(:,ind);
       x=min_decrease(x0, xU(:,ind_rp), inter_par);
       %x=xU(:,ind_rp);
       %x=round(x*Nm)/Nm;
       %keyboard
       x=General_linear_quatization(x, Ain, bin, Nm/length_scale, qtype);
       %keyboard
       if mindis(x,xU)<1e-4
      x=xU(:,ind_rp); xU(:,ind_rp)=[]; newadd=0;
       end
  else
    % tri=delaunayn([xE,xU]');
    % perform the trinagulation search
    while 1
    %    keyboard
     x=tringulation_search_bound(inter_par,[xE xU],yE);
     %keyboard
     %xq=round(x.*Nm)./Nm;
      [xq]=General_linear_quatization(x, Ain, bin, Nm/length_scale, qtype);
     if interpolate_val(x,inter_par)<y0
         break
     end
     [xq,xE,xU,newadd,success]=points_neighbers_find(xq,xE,xU);
     if success==1
         break
     else
         yup=[yup (interpolate_val(xq,inter_par)-y0)/mindis(xq,xE)];
         yrp=[yrp interpolate_val(xq,inter_par)]; 
     end
    end
        if (interpolate_val(x,inter_par)-y0)/mindis(x,xE)>min(yup)
            x=xU(:,ind_up); xU(:,ind_up)=[]; newadd=0;
       else
            x=xq;
        end
  end
      % Nm=50;
      [t,ind]=min(yE); x0=xE(:,ind);
      delta_mesh=1/Nm;
     % if (interpolate_val(x,inter_par)-y0)/(mindis(x,xE))>(interpolate_val(x0,inter_par)-y0)/(delta_mesh); 
     %     break
     % end
      ym=fun(x);
       
        %   true=check_add_point(x,ym,[xE xU],inter_par,newadd,size(xE,2));
           if mindis(x,xE)<DELTA_MAX
             %  keyboard
               break
           end
           %    keyboard
           %interpolate_index=[interpolate_index true];
           %if ym>2*(min(yE)-y0)+min(yE)
           %   interpolate_index(end)=0;
           %end
           %if ym<min(yE)
           %   interpolate_index(yE>2*(ym-y0)+ym)=0;
           %end  
           xE=[xE x]; yE=[yE ym];
           
          % interpolate_index=[interpolate_index true];
           pE=[pE interpolate_val(x,inter_par)];

   % illusteration
   
    % plotting
    figure(1)
    subplot(2,1,1)
    plot(1:length(yE(:)),yE,'-','linewidth',1.5)
%         ylim([0 200])
    subplot(2,1,2)
    plot(xE(:,1:end)','linewidth',1.5)
    drawnow
dt(k)=toc;
disp(strcat('finished iteration....', num2str(k*(kk-1)+k,'%02d'), '... ymin = ', num2str(min(yE)) ))
disp(strcat(' in ....  ', num2str(dt(k)), '... sec ' ))   

%      figure(gcf)
%     %subplot(2,1,1)
%     plot(1:length(yE),yE,'-','linewidth',2.5)
%     xlabel('number of fun eval')
%     ylabel('fun eval')
%     grid on
%     set(gca,'FontSize',18)
%     %, 1:length(yE),pE,'--')
%     %plot(1:length(yE),1./(exp(yE)-1),'-','linewidth',3)
%     ylim([0 200])
%     %ylim([-5 5])
%     %xlim([0 100])
%     %subplot(2,1,2)
%     %plot(xE')
%     %xlim([0 100])
%     drawnow
%       %plot(xi(1,:)0,xi(2,:),'rs')
%            
end
Nm=2*Nm;
end
%figure(2)
%plot(1:length(yE),yE(1:end),'-','linewidth',2)
%axis([0 12 -1 150])
%grid on
%set(gca,'fontsize',30)
%axes('Position',[0.7,0.3,0.28,0.28])
%plot(6:length(yE),yE(6:end),'-','linewidth',1.5)
%set(gca,'fontsize',15)
%grid on



% 
% 
%%%
figure(1);
subplot(2,1,1)
ylim([0,200])
% figure_to_publish(strcat('./results_Lambda/DD_Lambda_styb2_n_',num2str(n,'%02d'),'_',char(lattice)))
% figure_to_publish(strcat('./results_Lambda/DD_Lambda_styb2_n_',num2str(n,'%02d'),'_',char(lattice)))
[y_min,ind]=min(yE)
[xE(:,ind),xs]
% save(strcat('./results_Lambda/DD_Lambda_styb2_n_',num2str(n,'%02d'),'_','Dn_dual')) 

% save(strcat('./results_Lambda/DD_Lambda_styb2_n_',num2str(n,'%02d'),'_',char(lattice), 'spinterpolation_bug')) 
save(strcat('./results_Lambda/DD_Lambda_styb2_n_',num2str(n,'%02d'),'_',char(lattice), 'spinterpolation_bug')) 
