function [An,bn,Ain, bin]= normalize_Linear_constraints(A,b, lb, ub) 
% Normalize the constraints 
if (length(b)==0 && length(lb)==0)
  print('Error: The search doain is unbounded');
end
if length(b)~=0
B=diag(sqrt(A*A'));
bn=b./B;  An=diag(1./B)*A;
else
    bn=[]; An=[];
end

% Constrcut the inceremneted constraints
Ain=An; bin=bn;
%keyboard
   if length(lb)~=0
       n=length(lb);
       Ain=[Ain; eye(n); -eye(n)];
       bin=[bin; ub; -lb];
   end
end
