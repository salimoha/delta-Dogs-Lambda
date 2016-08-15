% 1D illustrative example for deltasearch
clear all
close all
clc

global n bnd1 bnd2
n=1; y0=2;
% The data points
xi=[0 20];
nn=20;
bnd1=0; bnd2=20;

% define the function
fun=@(x) cos(2*x)+sin(3*x)+x+3;
% interpolate method
inter_method=1;

% calculation
N=length(xi);
for i=1:N-1 xr(i)=(xi(i)+xi(i+1))/2; R2(i)=(-xi(i)+xi(i+1))^2/4; end
for i=1:N yi(i)=fun(xi(i)); end
inter_par= interpolateparametarization(xi,yi,inter_method);
xx=bnd1:(bnd2-bnd1)/(nn-1):bnd2;
for i=1:nn 
    yr(i)=fun(xx(i));
    yp(i)=interpolate_val(xx(i),inter_par);
    e(i)= max(R2- (xr-xx(i)).^2); 
    s(i)= (yp(i)-y0)/e(i);
end
[t,ind]=min(s);
% calculate the GPS
deltagps=0.2;
delta=mindis(xx(ind),xi);
ydelta=(interpolate_val(xx(ind),inter_par)-y0)/delta;
[y00,ind_min]=min(yi); x0=xi(:,ind_min);
xgps=repmat(x0,1,n+1)+deltagps*[-1 1];
for i=1:n+1
ygps(i)=interpolate_val(xgps(i),inter_par);
end
[yg,ind1]=min(ygps); yg=(yg-y0)/deltagps; xg=xgps(ind1);
xdelta=xx(ind);




% plot search-step figure
subplot(2,1,1)
plot(xx,yr,'k-','linewidth',3)
hold on
plot(xx,yp,'b--','linewidth',1.5)
hold on
plot(xi,yi,'ks','MarkerSize',10, 'MarkerFaceColor','k')
hold on
plot(xx,-1.7+3*e,'r-','linewidth',1)
hold on
plot(xx,xx*0+y0,'g-.','linewidth',2)
 set(gca,'YTick',[]);
  set(gca,'XTick',[]);
subplot(2,1,2)
plot(xx,s,'k-','linewidth',2)
hold on
plot(xx(ind),t,'ks','MarkerSize',10, 'MarkerFaceColor','k')
hold on
plot(xi,xi*0+2,'ks','MarkerSize',10, 'MarkerFaceColor','k')
hold on
text(xx(ind)-0.05,t-3.5,'$$\hat{x}$$','Interpreter','Latex','fontsize',30)
hold on
text((xx(ind)+0.8)/2,t+5,'$$\delta$$','Interpreter','Latex','fontsize',30)
ylim([2 20])
 set(gca,'YTick',[]);
  set(gca,'XTick',[]);
  
saveas(gcf,'deltasearch.eps','eps2c');

figure(2)
plot(xx,yr,'k-','linewidth',3)
hold on
plot(xx,yp,'b--','linewidth',1.5)
hold on
plot(xi,yi,'ks','MarkerSize',10, 'MarkerFaceColor','k')
hold on
plot(xdelta,1.2*ydelta,'ks','MarkerSize',10, 'MarkerFaceColor','b')
hold on
plot(xg,yg,'ks','MarkerSize',10, 'MarkerFaceColor','b')
%set(gca,'YTick',[]);
%  set(gca,'XTick',[]);
  text(xdelta+0.1,(ydelta+2)/2,'$$y_{\Delta}$$','Interpreter','Latex','fontsize',30)
  text(xg-0.6,(yg+2)/2,'$$y_{GPS}$$','Interpreter','Latex','fontsize',30)
  
  ylim([2 7])
