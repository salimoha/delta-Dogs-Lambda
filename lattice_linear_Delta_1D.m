% Second Order Delta-Dogs:
% Algorithm with second order convergence property.
% Constraints are simple bound constraints.
% Dimenson:
clear all
close all
clc

%setenvAVL

global n Ain bin y0 lattice
lattice ='Zn';

n=1;

% objective function:
% alpha=8;
% fun=@(x) 1./(9.*x(2,:)+0.1) +alpha*(x(1,:)-0.5).^2+(9.*x(2,:)+0.1); 
%fun=@ rastriginn;
%fun=@(x) rastriginn(3*(x-[0.3;0]));
%fun=@(x) sum(3*x)+sum(1./(3*x));
%fun=@(x) -20*exp(-0.2*sqrt(0.5*norm(x)^2))-exp(0.5*cos(2*pi*x(1))+0.5*cos(2*pi*x(2)))+exp(1)+20;
%fun= @(x) (1.5-x(1)+x(1)*x(2))^2+(2.25-x(1)+x(1)*x(2)^2)^2+(2.625-x(1)+x(1)*x(2)^3)^2;
%fun=@(x) log(1+(1.5-x(1)+x(1)*x(2))^2+(2.25-x(1)+x(1)*x(2)^2)^2+(2.625-x(1)+x(1)*x(2)^3)^2);
%fun=@(x)(x(1)-1)^2+10*(x(2)-x(1)^2)^2;
%fun=@ plygen_n;
% fun=@(x) sum(500*x.*sin(sqrt(abs(500*x))))/250;
%fun=@(x)
%(sum((-5+10*x).^4)-16*sum((-5+10*x).^2)+5*sum(-5+10*x)+39.16616*2*n);
%0<x<1
% fun=@(x) (sum((x).^4)-16*sum((x).^2)+5*sum(x)+39.16616*2*n);
%-5<x<5
%fun= @(x) sum((x-[0.3;0;-0.1;0.8;0.1;0;0.3]).^2);
%fun= @(x) -1/sum((x-[0.3;0;-0.1;0.8]).^2);
%fun= @(x) -1/sum((x-[0.3;-0.1]).^2)
%fun=@(x) sum(1./(0.2+10*x(1)))+sum(0.2+10*x(1))-2+(x(2)-0.1)^2;
% fun=@(x) sin(30.*(x-0.9).^4).*cos(2.*(x-0.9)) + (x-0.9)./2; y00=-0.625;
fun=@(x) sin(30.*((x-0.1)-0.9).^4).*cos(2.*((x-0.1)-0.9)) + ((x-0.1)-0.9)./2; y00=-0.8421;
% Constraints:
%lob=zeros(n,1); upb=ones(n,1);
%lob=[-2;-2]; upb=[2*pi;2*pi];
%funr=@(x) sum(500*x.*sin(sqrt(abs(500*x))))/250;
% upper nad lower bounds
%lob=zeros(n,1); upb=ones(n,1);
lob=0*ones(n,1); upb=1*ones(n,1);

%lob=[-4.5;-4.5]; upb=[4.5;4.5];
%lob=[-2;-2]; upb=[3.5;3.5];
%lob =0.1*ones(n,1); upb=2.2*ones(n,1);

% Estimate of the solution:
%y00=0.01;
%y00=4;
%y00=-100;
%y00=0.027;
%y00=1/37;
%y00=0;
%y00=-37;

% Initial Calculation:
bnd1 = lob; bnd2 = upb;
m=2*n;

% maximum number of iterations:
iter_max=50;

% interpolation strategy:
inter_method=7;
% inter_method=1;

% Input the equality constraints:?????????
Ain=[eye(n);-eye(n)]; bin=[bnd2 ;-bnd1];
% Ain=[eye(n); -ones(1,n)/sqrt(n)]; bin=[5*ones(n,1); 2.8*n/sqrt(n)];
length_scale=1;
% length_scale=10;
%length_scale=max(bnd2-bnd1);
% Calculate the Initial trinagulation points
%xU=bounds(bnd1,bnd2, n);
xU=vertex_find(Ain,bin,[],[]);
% Calculate initial delta
% xE=(bnd1+0.45*(bnd2-bnd1));
% xE=(sum(xU')/(n+1))';
%xE=[0.350 ;0.4000 ;0.5000 ;0.5000 ;0.3000 ;0.3500; 0.3000];
% initial points: the midpoint point and its neighber
xE=[0.5];
delta0=0.2*length_scale;
for ii=1:n
    e=zeros(n,1); e(ii)=1;
    xE(ii+1)=xE(1)+delta0*e;
end

% Calculate initial unevaluated delta
for ii=1:size(xU,2)
    deltaU(ii)=mindis(xU(:,ii),xE);
end
% Calculate the function at initial points
for ii=1:size(xE,2)
% % %     yE(ii)=fun(xE(:,ii));
     yE(ii)=fun(xE(ii));
end
pE=yE;
interpolate_index=ones(1,length(yE));
%Nm=800;
Nm=20;
%y00=0;
%y00=0.02;
% Initialize constaints 
%tol=0.05;
rho=1;
figure(3);clf;
xx=0:0.01:1;
plot(xx,fun(xx)); hold on

for kk=1:3

for k=1:iter_max
%    keyboard
    stop=0;
    inter_par= interpolateparametarization(xE,yE,inter_method,interpolate_index);
    %y0=10;
    y0=y00;
    %y0=min(yE)-range(yE)/Nm;
    yup=[]; 
    yrp=[];
   % Calculate the metric of the unevaluated points.
            for ii=1:size(xU,2)
               yup(ii)=(interpolate_val(xU(ii),inter_par)-y0)/mindis(xU(ii),xE);
               yrp(ii)=interpolate_val(xU(ii),inter_par);
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
      [tt,ind]=min(yE); x0=xE(ind); 
%       x=x0;
       x=min_decrease(x0, xU(ind_rp), inter_par);
       %keyboard
       %x=xU(:,ind_rp);
       
       %x=round(x*Nm)/Nm;
       %keyboard
       x=General_linear_quatization(x, Ain, bin, Nm/length_scale, 1);
       %keyboard
       if mindis(x,xU)<1e-4
      x=xU(ind_rp); xU(ind_rp)=[]; newadd=0;
       end
  else
    % tri=delaunayn([xE,xU]');
    % perform the trinagulation search
    while 1
    %    keyboard
     x=tringulation_search_bound(inter_par,[xE xU],yE);
%      keyboard
     %xq=round(x.*Nm)./Nm;
      [xq]=General_linear_quatization(x, Ain, bin, Nm/length_scale, 1);
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
        if (interpolate_val(x,inter_par)-y0)/mindis(x,xE)>tup
           % keyboard
            %x=round(x*Nm)/Nm;
           x=General_linear_quatization(x, Ain, bin, Nm/length_scale, 1);
           % keyboard
            x=xU(ind_up); xU(ind_up)=[]; newadd=0;
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
           if mindis(x,xE)<5e-16
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
%      figure(1)
    %subplot(2,1,1)
%     plot(1:length(yE),yE,'-')
    %, 1:length(yE),pE,'--')
    %plot(1:length(yE),1./(exp(yE)-1),'-','linewidth',3)
    %ylim([0 0.03])
    %ylim([-5 5])
    %xlim([0 100])
    %subplot(2,1,2)
    %plot(xE')
    %xlim([0 100])
pause(1)
    figure(3);
    
plot(xE,yE, 'k*')
    drawnow
      %plot(xi(1,:)0,xi(2,:),'rs')
           
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

