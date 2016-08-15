clear all
close all
clc

% illusterate the issue of e(x)
xi=[0 0 1 1; 0 1 0 1];
%xi=[xi rand(2,5)];
%xi=[xi [0.2;0.5] [0; 0.5]];
xi=[xi [0.2;0.5] [0; 0.7]];
%xi=[xi [0.5;0.5]];
tri=delaunayn(xi.');
%A=A';
%[xc,R2]=circhyp2(A,2);
%plot(A(:,1),A(:,2),'ks')
%axis square
xlim([0 1])
ylim([0 1])

N=100;
for i=0:N
for j=0:N
xv=[i ;j]/N;
e=0;
for ind=1:size(tri,1)
  [xc,R2]=circhyp2(xi(:,tri(ind,:)),2);  
  e(ind)=R2-norm(xv-xc)^2;
end
U(i+1,j+1)=max(e);
end
end

 contourf(0:1/N:1,0:1/N:1,U.')
 %colormap('bone')
 xlim([0 1])
 ylim([0 1])
 hold on
plot(xi(1,:),xi(2,:),'ks','markersize',10,'markerFacecolor','k')
 
%DT = delaunayTriangulation(A.');
DT = delaunayTriangulation(xi.');
triplot(DT,'w-','linewidth',2)
surf(0:1/N:1,0:1/N:1,U.')
colormap('bone')
view(140,40)
alpha(0.9)
box off

figure_to_publish('./err_plot_iluss')


