function g = interpolate_grad(x,inter_par)
% Calculate te interpolatated value at points x
% inter_par{1}=1 polyharmonic spline
% inter_par{1}=2 Quadratic interpolation

n=length(x);

% polyharmonoic spline
if inter_par{1}==1
    w=inter_par{2}; v=inter_par{3};
    xi=inter_par{4}; 
    N = size(xi, 2);
    g = zeros(n, 1);
for ii = 1 : N
    X = x - xi(:,ii);
    g = g + 3 *w(ii)* X*norm(X);
end
    g=g+v(2:end);
end

% Quadratic interpolation

if inter_par{1}==2
    k=inter_par{2}; v=inter_par{3};
    xmin=inter_par{4};
    g=2*k'*(x-xmin)+v;
end

if inter_par{1}==3
   % ymin=inter_par{2}; 
    xmin=inter_par{3};
    f=inter_par{4};
    H=inter_par{5};
    g=2*H*(x-xmin)+f;
end

% scaled polyharmonoic spline
if inter_par{1}==7 || inter_par{1}==8
% keyboard
    w=inter_par{2}; v=inter_par{3};
    xi=inter_par{4};  a = inter_par{7}; H = diag(a);
    N = size(xi, 2);
    g = zeros(n, 1);
%     keyboard
for ii = 1 : N
    X = x - xi(:,ii);
%     g = g + 3 *w(ii)* X'*norm(X);
        g = g + 3*w(ii)*H*X*(X'*H*X).^(1/2);
               term(ii,:)=  3*H*X*(X'*H*X).^(1/2);
%       dA(ii,jj,:) =3/2.* (xi(:,ii) - xi(:,jj)).^2 *  ((xi(:,ii) - xi(:,jj))' *H* (xi(:,ii) - xi(:,jj)))^(1/2)
end
    g=g+v(2:end);
end



end