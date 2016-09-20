function [y,inter_par] = interpolate_val(x,inter_par)
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
    S = bsxfun(@minus, xi.*repmat(sqrt(a),1,size(xi,2)), x.*sqrt(a));
    y = v' * [1; x] + w' * sqrt(diag(S' * S)) .^ 3;
end

if inter_par{1}==8
    w=inter_par{2};
    v=inter_par{3};
    xi=inter_par{4};
    a = inter_par{7}; 
    a0=a; % to be used for evaluting the results from inter_method=7
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           if the distance from target value y0 was greater than the error
%           in the calculation of collapsed interpolation compare to inter_method
%           7 or original polyharmonic spline then use the collapsed
%           interpolation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%  Calculating the collapsed interpolation
    for ii=1:length(a)
       if a(ii)<0.1,
          a(ii) = 0; 
       end     
    end
    inter_par{7}=a;
%     keyboard
    S = bsxfun(@minus, xi.*repmat(sqrt(a),1,size(xi,2)), x.*sqrt(a));
    y = v' * [1; x] + w' * sqrt(diag(S' * S)) .^ 3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%  Calculating the scaled interpolation
    S0 = bsxfun(@minus, xi.*repmat(sqrt(a0),1,size(xi,2)), x.*sqrt(a0));
    y0 = v' * [1; x] + w' * sqrt(diag(S0' * S0)) .^ 3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

inter_err=max(abs(y-y0));

% do not collapse
if min(inter_par{8})< inter_err 
y=y0;
inter_par{7}=a0;
end


end

