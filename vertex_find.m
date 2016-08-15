function [Vertex]=vertex_find(A,b,lb,ub)
% Find the vertice sof a simplex:

% Ax<b: Linear constraints:
% lb<x<ub: actual constraints:

%keyboard

if length(lb)~=0
Vertex=[];
[m,n]=size(A);
if m==0
Vertex=bounds(lb, ub, length(lb));  
else
for r=0:min(m,n)
%     keyboard
    C = nchoosek(1:m,r);
    D = nchoosek(1:n,n-r);
    if r==0
         F=bounds(lb, ub, n);
         for kk=1:size(F,2)
           x=F(:,kk);
           if min(A*x-b)<1e-6
           Vertex=[Vertex x];
           end
        end
         
    else
       for ii=1:size(C,1)
         index_A=C(ii,:);
         index_A_C=setdiff(1:m,index_A);
         A1=A(index_A,:);
         b1=b(index_A);
        for jj=1:size(D,1)
        index_B=D(jj,:); index_B_C=setdiff(1:n,index_B);
        F=bounds(lb(index_B), ub(index_B), n-r);
        A11=A1(:,index_B); A12=A1(:,index_B_C);
        for kk=1:size(F,2)
           A11=A1(:,index_B); A12=A1(:,index_B_C);
           xd=A12\(b1-A11*F(:,kk));
           x(index_B)=F(:,kk);
           x(index_B_C)=xd;
          % keyboard
           if (r==m || min(A(index_A_C,:)*x-b(index_A_C))<0)
               if (max(x-ub)<1e-6 && min(x-lb)>-1e-6 && mindis(x,Vertex)>1e-6)
               Vertex=[Vertex x];
               end
           end
         end
        end
        end
        end
end
end
else
    [m,n]=size(A);
    C = nchoosek(1:m,n);
    Vertex=[];
     for ii=1:size(C,1)
         index_A=C(ii,:);
         index_A_C=setdiff(1:m,index_A);
         A1=A(index_A,:); b1=b(index_A);
         A2=A(index_A_C,:); b2=b(index_A_C);
         x=A1\b1;
         if max(A2*x-b2)<1e-6
           Vertex=[Vertex x];
         end
         
     end
    
end
end



