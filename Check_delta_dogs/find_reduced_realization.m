lb=[0.2;0.5;0.5;-0.3;0.1;0.05;0.1]; 
ub=[0.5;1.5;1.5;0.3;0.3;0.5;2];

for ii=n+2:size(xiT,2)
 xi_real(ii-n-1,:)=((xiT(:,ii).*(ub-lb)+lb))';
 xi_real(ii-n-1,5)=1+1/xi_real(ii-n-1,5);
 xi_real(ii-n-1,7)=1+1/xi_real(ii-n-1,7);

end
tt=yiT(n+2:end);
xi_real(end,4)=-0.3;
ind=1;
for ii=1:length(tt)
if tt(ii)<min(tt(1:ii-1))
  ind=[ind,ii];
end
end
ind=ind(ind>10);
xi_real_reduced=xi_real(ind,:);