clear all
close all
clc

n=2;

% illusterate the issue of e(x)
xi=[0 0 1 1; 0 1 0 1];
xi=[xi rand(2,5)];
%xi=[xi [0.2 0;0.5 0.7]];
%xi=[0 1 0.2 0.7];
tri=delaunayn(xi.');
%A=A';
%[xc,R2]=circhyp2(A,2);
%plot(A(:,1),A(:,2),'ks')
%axis square
xlim([0 1])
ylim([0 1])

N=50;
for i=0:N
for j=0:N
xv=[i ;j]/N;
%xv=i/N;
e=[];
for ind=1:size(tri,1)
  [xc,R2]=circhyp2(xi(:,tri(ind,:)),1);  
  
  e(ind)=R2-norm(xv-xc)^2;
end
U(i+1,j+1)=max(e);
%keyboard
%U(i+1)=max(e);
end
end
if n==1
plot(0:1/N:1, U, 'k-', 'linewidth',3)
hold on
plot(xi,0*xi,'ks','markersize', 10, 'MarkerFaceColor','k')
%set(gca,'XTick',[]); 
%set(gca,'YTick',[]);
grid on
end
if n==2
 %contourf(0:1/N:1,0:1/N:1,U.')
 surf(0:1/N:1,0:1/N:1,U.')
 %colormap('bone')
 xlim([0 1])
 ylim([0 1])
 hold on
 %plot(A(1,:),A(2,:),'ks','markersize',10,'markerFacecolor','w')
 
% DT = delaunayTriangulation(xi.');
%triplot(DT,'k-','linewidth',2)
%plot(xi(1,:),xi(2,:), 'ks', 'linewidth',20)
end
%surf(0:1/N:1,0:1/N:1,U.')



