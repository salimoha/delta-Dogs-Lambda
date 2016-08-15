clear all 
close all
clc

% check the regression
 n=1; N=6;
 fun=@(x) x.^2;
 xi=[0 0.1 0.2 0.3 0.5 1];
 for ii=1:N
 yi_real(ii)=fun(xi(:,ii));
 end
 sigma=[0.01 0.2 0.01 0.2 0.01 0.2];
 V=randn(1,N).*sigma;
 yi=yi_real+V;
 
 
 rho=0:0.001:0.5;
 
 [inter_par,yp]= regressionparametarization(xi,yi,sigma,1);

 %figure(1)
 %plot(rho,weight)
 
 %[t,ind]=min(abs(weight-1));
 
 xx=0:0.1:1;
 for ii=1:length(xx)
     yp1(ii)=interpolate_val(xx(ii),inter_par);
  %    yp2(ii)=interpolate_val(xx(ii),inter_par2);
  %     yp3(ii)=interpolate_val(xx(ii),inter_par3);
 end
 
figure(2)
plot(xx,yp1,'b--',xx,xx.^2,'k-')
hold on
plot(xx,xx.^2,'k-')
 hold on
 errorbar(xi,yi,sigma,'.')
% xlim([0 1])