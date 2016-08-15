function [xq]= General_unconstrained_quantization(x, N, qtype)
 % N: meshize
 % qtype: 1, cartesian, 2: Dn, 3: An.
 
 x=x*N;
 
 %keyboard
 
 % Cartesian
 if qtype==1
     xq=round(x);
 end
 
 % Checkerboard
 if qtype==2
     xq=round(x);
     if mod(sum(xq),2)==1
     dq=x-xq; [t,ind]=max(abs(dq));
     xq(ind)= xq(ind)+sign(dq(ind)+eps);
     end
 end
 
 % Zero sum
 if qtype==3
     
     n=length(x);
  %   keyboard
     V=null(ones(1,n+1));
     x=V*x;
     xq=round(x);
     delta=sum(xq);
     if delta~=0
     d=x-xq; [~,ind]=sort(d);
        if delta>0
           xq(ind(1:-delta))=xq(ind(1:-delta))-1;
        else
            xq(ind(n-delta+2:n+1))=xq(ind(n-delta+2:n+1))+1; 
        end
     end
     xq=V\xq;
 end
 
 xq=xq/N;  
end