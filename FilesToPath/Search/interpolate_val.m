function y = interpolate_val(x,inter_par)
% Calculate te interpolatated value at points x
% inter_par{1}=1 polyharmonic spline
% inter_par{1}=2 Quadratic interpolation
% inter_par{1}=7 Scaled Polyharmonic interpolation

if inter_par{1}==1
    w=inter_par{2};
    v=inter_par{3};
    xi=inter_par{4};
    S = bsxfun(@minus, xi, x);
    y = v' * [1; x] + w' * sqrt(diag(S' * S)) .^ 3;
end
if inter_par{1}==2
    k=inter_par{2};
    v=inter_par{3};
    xmin=inter_par{4}; ymin=inter_par{5};
    y=ymin+(x-xmin)'*v+k'*((x-xmin).^2);
end
if (inter_par{1}==3 || inter_par{1}==4)
    
   % keyboard
    ymin=inter_par{2};
    xmin=inter_par{3};
    f=inter_par{4};
    H=inter_par{5}; 
 %   keyboard
    y=(x-xmin)'*H*(x-xmin) +f'*(x-xmin)+ymin;
end
% keyboard
if inter_par{1}==7
    w=inter_par{2};
    v=inter_par{3};
    xi=inter_par{4};
    a = inter_par{7};
%      keyboard
    S = bsxfun(@minus, xi.*repmat(sqrt(a),1,size(xi,2)), x.*repmat(sqrt(a),1,size(x,2)));
    y = v' * [1; x] + w' * sqrt(diag(S' * S)) .^ 3;
end



end

