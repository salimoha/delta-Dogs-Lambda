% Illusteration of generalized grid
% Domain the uniform simplex
clear all; close all; clc
n=3;
% N=7;
N=6; M=N; Tol=1/N;
[V,a] = uniformsimplexrecursive(n);
V=V/a;
A=-V'; b=ones(n+1,1)/n;
% check
%V1=vertex_find(A,b,[],[]);
% 
% Defining B
B=eye(n);
% L_N,{}
datapoints=[];
for ii=-M:M
    for jj=-M:M
        if n==3
        for kk=-M:M
            x=[ii,jj,kk]';
            x=1/N*B*x;
            if b-A*x>Tol
            datapoints=[datapoints, x];
            end
        end
        else
            x=[ii,jj]';
            x=1/N*B*x;
            if b-A*x>Tol
            datapoints=[datapoints, x];
            end
        end
    end
end
%%
figure(1);clf;
if n==3
hold on
DT = delaunayTriangulation(V(1,:)',V(2,:)',V(3,:)');
tetramesh(DT);
alpha(0.3)
view(-107,30)
scatter3(datapoints(1,:),datapoints(2,:),datapoints(3,:),'MarkerFaceColor',[0 .75 .75])
patch('Vertices',V')
end
%
if n==2
hold on
plot(datapoints(1,:),datapoints(2,:),'ks')
f = [1 2 3];
patch('Faces',f,'Vertices',V');
alpha(0.3)
end