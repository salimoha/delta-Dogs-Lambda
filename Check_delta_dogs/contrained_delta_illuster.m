% Representation of subinterpolation 
clear all
close all
clc
warning('off','all')



fun=@(x) x.^2+sin(10*x);

gfun=@(x) x.^2+sin(10*(x+0.1));
% xi=[0.001 1 0.2 0.3 0.5 0.6 0.12 0.07 .17 .15];
xi=[0.1 0.7 0.2];
for ii=1:length(xi)
    yi(ii)=fun(xi(:,ii));
    gi(ii)=gfun(xi(:,ii));
end
%T=ones(1,3);
%K=2; L=1;
%[xi,ind]=sort(xi); yi=yi(ind);
xx=0:0.01:1;
for ii=1:length(xx)
yy(ii)=fun(xx(ii));
gg(ii)=gfun(xx(ii));
end
figure(1)
plot(xx,yy,'k-',xx,gg,'k--',xx,max(yy,gg)+2,'k-.')
xc=(xi(1:end-1)+xi(2:end))/2;
d=(xi(2:end)-xi(1:end-1))/2;
%[inter_par,yp]=interpolateparametarization(xi,yi,sigma0./sqrt(T),1);
[inter_par]=interpolateparametarization(xi,yi,1);
[inter_par_g]=interpolateparametarization(xi,gi,1);
%x=0.01; e=max(d.^2-(x-xc).^2);


%%
for ii=1:length(xx)
    x=xx(ii);
    %fr(ii)= funr(x);
    f(ii)=fun(x);
    p(ii) = interpolate_val(x,inter_par);
    pg(ii) = interpolate_val(x,inter_par_g);
    e(ii) = max(d.^2-(x-xc).^2);
    sc(ii) = max(p(ii)-;
    %sn(ii)= subinterpolate_search(x,e(ii),xi,yi,K);
    %sn(ii)=subquad_search(x,e(ii),xi,yi,K);
end
%[tc,indc]=min(sc);
figure(2)
plot(xx,pg,'k-',xx,p,'k--', xx, max(p,pg)+2, 'k-.') 
hold on
grid on
% errorbar(xi,yi,sigma0./sqrt(T),'k.','linewidth',2)
%    sd = min(yp,2*yi-yp)-L*sigma0./sqrt(T);
%[td,indd]=min(sd);
%hold on
plot(xi,yi,'ks', 'Markersize',10, 'MarkerFacecolor', 'k')
plot(xi,gi,'ks', 'Markersize',10, 'MarkerFacecolor', 'k')

%text(xx(indc),tc-.05,'x_k','fontsize',20)
plot(xx(indc),tc,'k*','Markersize',10)
%text(xi(indd),td-.05,'z_k','fontsize',20)
%set(gca,'YTickLabel',[])
%set(gca,'XTickLabel',[])
% ylim([-.3 10])
grid on
    


[t,ind]=min(sn);
ind
xx(ind)
x=xx(ind);
ii=ind;
[y,abc]=subquad_search(x,e(ii),xi,yi,K)
y_inter=abc(1)*xx.^2+ abc(2)*xx+abc(3);
hold on
plot(xx,y_inter,'r-')