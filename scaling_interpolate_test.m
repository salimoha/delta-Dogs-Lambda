% qudrtaic test for polyharmonic spline
clear all
close all
clc

n=2;
x0=0.5*ones(n,1);
xi=[x0 repmat(x0,1,n)+eye(n)*0.3 repmat(x0,1,n)-eye(n)*0.3];
%fun=@(x) norm(x-0.45*ones(n,1)).^2;
fun=@(x) (x(1)-0.45)^2+0.01*(x(2)-0.45)^2;
for ii=1:size(xi,2)
    yi(ii)=fun(xi(:,ii));
end

inter_par= interpolateparametarization(xi,yi,1);

xv=0:0.05:1;
for ii=1:length(xv)
   %  for jj=1:length(xv)
  yr(ii)=fun([xv(ii);0.5]);
  yp(ii)=interpolate_val([xv(ii);0.5],inter_par);
  yr1(ii)=fun([0.5;xv(ii)]);
  yp1(ii)=interpolate_val([0.5;xv(ii)],inter_par);
     for jj=1:length(xv)
  Ur(ii,jj)=fun([xv(ii);xv(jj)]);
  Up(ii,jj)=interpolate_val([xv(ii);xv(jj)].*[a;1],inter_par);
     end
end
sum(inter_par{2}.^2)
a=8;
inter_par= interpolateparametarization(xi.*repmat([a;1],1,size(xi,2)),yi,1);
%sum(inter_par{2}.^2)*norm([1;a])^6;
    for ii=1:length(xv)
    Ur11(ii)=interpolate_val([0.5*a;xv(ii)],inter_par);
    Ur12(ii)=interpolate_val([xv(ii)*a;0.5],inter_par);
    end


figure(1)
plot(xv,yr,'r-',xv,yp,'b--')

figure(2)
plot(xv,yr1,'r-',xv,yp1,'b--')


figure(3)
plot(xv,Ur11,'b--',xv,yr1,'r-')

figure(4)
plot(xv,yr,'r-',xv,Ur12,'b--')

for a=1:1:100
inter_par= interpolateparametarization(xi.*repmat([a;1],1,size(xi,2)),yi,1);
cost(a)=sum(inter_par{2}.^2)*norm([1;a])^3;
end

figure(5)
plot(1:1:100,cost)
%for ii=1:length(xv)
%    for jj=1:length(xv)
%Ur(ii,jj)=fun([xv(ii) ;xv(jj)]);
%Up(ii,jj)=interpolate_val([xv(ii) ;xv(jj)],inter_par);
%end
%end
