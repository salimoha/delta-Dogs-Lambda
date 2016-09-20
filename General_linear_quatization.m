function [xq]=General_linear_quatization(x, A, b, N, qtype)
% Quantization of x on the mesh of size N:
% L={x| Ax \le b}
% Normalize the matrix A, and vector b.
% N:mesh size, delta \rightarrow \infty as N \rightarrow \infty.
% qtype=0.
%keyboard
[m,n]=size(A);
% find active constrAts at x.
active_constrAts=1:m;
active_constrAts=active_constrAts(b-A*x<1e-6);
inactive_constrAts=setdiff(1:m,active_constrAts);
while 1
Aeq=A(active_constrAts,:); beq=b(active_constrAts,:);
Aineq=A(inactive_constrAts,:); bineq=b(inactive_constrAts,:);
% Find the x0.
x0=Aeq\beq; 
% if x0 is a vertex stop
if length(active_constrAts)>n-1
    xq=x0;
    break
end
% x= x=x0+ Vr
V=null(Aeq); r=V\(x-x0); 
% 
% r= General_unconstrained_quantization(r, N, qtype);
% keyboard
r = Unconstraint_quantizer(r,N);
dis_ineq=bineq-Aineq*(x0+V*r);
[t,ind]=min(dis_ineq);
if t>1/N
    xq=x0+V*r;
    break
else
  dis_ineq=bineq-Aineq*x; [t,ind]=min(dis_ineq);
  %keyboard
  x=x+t*Aineq(ind,:)';
  active_constrAts=[active_constrAts inactive_constrAts(ind)];
  inactive_constrAts(ind)=[];
end

end






    






