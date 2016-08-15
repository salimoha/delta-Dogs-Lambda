close all
tic
load('result_avl_noscale.mat')
inter_par = interpolateparametarization(xE,yE,7,interpolate_index);
toc
[ymin,ind] = min(yE); ymin = 1/ymin;
xm = xE(:,ind);
figure(17),
subplot(2,1,1)
tt=0:0.01:1;
for ii=1:n
%  figure,
   xx=xm;
   for jj=1:length(tt)
       xx(ii) = tt(jj); 
       ysp(ii,jj)=interpolate_val(xx,inter_par);
       
   end
xx = (tt+2*(ii-1))./2+0.7;
xmin = (xm(ii)+2*(ii-1))./2+0.7;
   hold on
  plot(xx,1./ysp(ii,:) )  
  plot(xmin,ymin, 'pk')
end

%% standard polyharmonic spline
inter_par01 = interpolateparametarization(xE,yE,1,interpolate_index);
% figure(18)
subplot(2,1,2)
for ii=1:n
%  figure,
   xx=xm;
   for jj=1:length(tt)
       xx(ii) = tt(jj); 
       yy(ii,jj)=interpolate_val(xx,inter_par01);
   end
   xx = (tt+2*(ii-1))./2+0.7;
   xmin = (xm(ii)+2*(ii-1))./2+0.7;
   hold on
  plot(xx,1./yy(ii,:) )   
   plot(xmin,ymin, 'pk')
end