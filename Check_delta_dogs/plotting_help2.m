 for ii=n+2:size(xiT,2)
      if yiT(ii)<interpolate_val(xiT(:,ii),inter_par)
          error(ii)=yiT(ii)-interpolate_val(xiT(:,ii),inter_par);
      elseif interpolate_val(xiT(:,ii),inter_par)<min(yiT)
          error(ii)=min(yiT)-interpolate_val(xiT(:,ii),inter_par);
      end
 end
   
 
 
 for ii=1:n
     for jj=1:n
     e=zeros(1,n^2);
     e((jj-1)*n+ii)=1;
     pi((ii-1)*n+jj,:)=e;
     end
 end