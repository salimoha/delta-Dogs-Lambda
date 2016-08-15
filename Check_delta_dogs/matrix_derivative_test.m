%clear all
%close all


%A=rand(2,2); x=rand(2,1);
eps=1e-3;
y=error_cost(d,f,A,xi,yi,w);
[Dd,Df,DA_a]=error_grad(d,f,A,xi,yi,w);
for ii=1:2
    for jj=1:2
        %DA_a(ii,jj)=2*x(jj)*A(ii,:)*x; 
        A1=A; A1(ii,jj)=A1(ii,jj)+eps;
       % DA_n(ii,jj)= (norm(A1*x)^2-norm(A*x)^2)/eps;
   DA_n(ii,jj)= (error_cost(d,f,A1,xi,yi,w)-y)/eps;
    end
end

DA_a
DA_n 


