% Visulatisation
xv=0:0.01:1; 
for ii=1:length(xv)
    for jj=1:length(xv)
xx=[xv(ii); xv(jj)];
yr(ii,jj)=fun(ss0+rk.*xx);
    end
end


contourf(xv,xv,yr,0:1:10)
hold on
plot(xiT(1,n+2:end),xiT(2,n+2:end),'ks','MarkerFacecolor','w','markersize',10)

