clear all
close all
clc

X=[];
N=4;

for ii=1:N
    for jj=1:N
        X=[X [ii;jj]];
    end
end

corners=[1 1 N N;1 N 1 N];
xi=[1 N 2; 2 2 N];
xk=[N 3.2];
yk=[N 3];
zk=[N N];
%plot(X(1,:),X(2,:),'ks','markersize',20)
plot(corners(1,:),corners(2,:),'k*','markersize',20)
hold on
plot(xi(1,:),xi(2,:),'ks','markersize',20)
plot(xk(1),xk(2),'ko','markersize',20,'markerfacecolor','k')
text(xk(1)+0.1,xk(2),'x_k','fontsize',20)
plot(yk(1),yk(2),'kv','markersize',20,'markerfacecolor','k')
text(yk(1)+0.1,yk(2),'y_k','fontsize',20)
text(zk(1)+0.1,zk(2)-0.1,'z_k','fontsize',20)
hold on
%plot(N,1.8,'k*','markersize',20)
%text(N+0.2,1.7,'x','fontsize',30)
%plot(1.5,1.5,'k*','markersize',10)
%text(N+0.2,2.1,'x_q','fontsize',30)
%text(N,2,'x_q','fontsize',10)
axis square
for ii=1:N
plot(1:N,(1:N)*0+ii,'k-');
hold on
plot((1:N)*0+ii,(1:N),'k-');
end
set(gca,'YTick',[])
set(gca,'XTick',[])
%h = text(1.6,1.8, '\delta_N','fontsize',50);
%set(h, 'rotation', 45);

%text(1.75,1.5,'\delta_{L^n}','fontsize',N0)
%annotation('doublearrow',[.N2 .55],[0.N2 0.55])