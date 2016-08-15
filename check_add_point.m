function true=check_add_point(x,new_val,xi,inter_par,newpoint,evaluated_length)
%keyboard
global n 
e=0;
true=1;
%keyboard
if newpoint
    tri=delaunayn(xi.');
    % Find the uncertaity function at n_new
   for ii=1:size(tri)
   [xc,R2]=circhyp(xi(:,tri(ii,:)), n);
   if e<R2-norm(x-xc)^2
   e=R2-norm(x-xc)^2; tri_index=ii;
   end
   end
    % Find the linear estimator 
   X=xi(:,tri(tri_index,:));
    for ii=1:n+1
       Y(ii)=interpolate_val(X(:,ii),inter_par);   
    end
  %  keyboard
    wm=[ones(1,n+1);X]\[1 ;x]; linear_estimate=Y*wm;
    % check the new curvature
    L_hess=0.1;
    if (new_val-linear_estimate)/e>L_hess
       true=0;
    end
else
  %  keyboard
    [delta,ind]=mindis(x,xi(:,1:evaluated_length));
    L_linear=0.5;
    if  (new_val-interpolate_val(xi(:,ind),inter_par))/delta>L_linear
    true=0; 
    end
end
 

end