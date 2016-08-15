for ii=n+2:size(xiT,2)
      p(ii)=interpolate_val(xiT(:,ii),inter_par);
end

x0=0.5*ones(n,1);
clear pp
xv=-4.5:0.2:4.5;
for ii=n+1:length(yiT)
  piE(ii)=interpolate_val(xiT(:,ii),inter_par);
end

w=0:0.1:1;
for ii=1:length(w)
    x1=w(ii)*xE(:,8)+(1-w(ii))*xE(:,10);
    yp(ii)=interpolate_val(x1,inter_par);
    U(ii)=fun(x1);
end


for ii=1:length(xv)
    for jj=1:length(xv)
        x=[xv(ii);xv(jj)];
     % pp(ii,jj)=interpolate_val(x,inter_par);
      U(ii,jj)=fun(x);
    end
end
%surf(pp)
contourf(-5+10*xv,-5+10*xv,U.',0:10:100)
colormap('bone')
hold on
 plot(-5+xE(1,:)*10,-5+10*xE(2,:),'ks','markersize',10,'markerfacecolor','w')
 plot(-5+10*xU(1,:),-5+10*xU(2,:),'k*','markersize',10)
 %plot(-5+xiT(1,4:end)*10,-5+10*xiT(2,4:end),'ks','markersize',10,'markerfacecolor','w')