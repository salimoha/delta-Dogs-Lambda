% zerolevel curve test
clear all
close all
clc

x0=[0.5;0];
x1=[0;0.5];
R2=0.25;
xv=0:0.02:1;
for ii=1:length(xv)
    for jj=1:length(xv)
        x=[xv(ii);xv(jj)];
    con1(ii,jj)=R2-norm(x-x0)^2;
    con2(ii,jj)=R2-norm(x-x1)^2;
    con(ii,jj)=max(con2(ii,jj),con1(ii,jj));
    end
end

contourf(xv,xv,con',0:2:2)
colormap('bone')
brightness('0.2')