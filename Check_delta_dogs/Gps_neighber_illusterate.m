% Illusteration of GPS point in the 2D
clear all
close all
clc


N = 4;
x = linspace(0,1,N+1)
y = linspace(0,1,N+1)

% Horizontal grid 
for k = 1:length(y)
  line([x(1) x(end)], [y(k) y(k)])
end

% Vertical grid
for k = 1:length(y)
  line([x(k) x(k)], [y(1) y(end)])
end
axis square

hold on
plot(0.5,1,'k*','MarkerSize',10, 'MarkerFaceColor','b')
hold on
plot(0.75,1,'ks','MarkerSize',10, 'MarkerFaceColor','r')
hold on
plot(0.25,1,'ks','MarkerSize',10, 'MarkerFaceColor','r')
hold on
plot(0.5,0.75,'ks','MarkerSize',10, 'MarkerFaceColor','r')
hold on
text(0.52,1.02,'$$x_0$$','Interpreter','Latex','fontsize',20)

set(gca,'YTick',[]);
set(gca,'XTick',[]);