n=2;
% A= ones(1,n)*b;
a=5;b0=2.9;
b=2; % interor
% b=a-b0; %active on the constraint
% b=1.9
flin1=@(x) n*b-n*a-x;
flin2=@(x)(n*a-n*b-x);
xx=[-5:0.1:5]';
funplotting2d(@stybtang,xx,xx,n,0)
figure(2);
hold on 
plot([flin2(a),a],[a,flin2(a)],'k-','linewidth',2)
plot([-a,flin1(-a)],[flin1(-a),-a],'k-','linewidth',2)
rectangle('Position', [-a -a 2*a 2*a])
plot(-b0,-b0,'r*')
% figure_to_publish('./figures/example3_active_const')
% figure_to_publish('./figures/example3_interior')

%% scaled
% bs=2;
% Ain=[eye(n);-eye(n);ones(1,n);-ones(1,n)]
% bin=[ones(n,1)*5; ones(n,1)*5;bs;bs]
%%
n=2;
% A= ones(1,n)*b;
a=1;b0=-0.2097;
b=1-0.19; % interor
b=a+b0; %active on the constraint
% b=1.9
flin1=@(x) n*b-x;
flin2=@(x)(n*a-n*b-x);
xx=[0:0.01:a]';
fun=@(x)(sum((-5+10*x).^4)-16*sum((-5+10*x).^2)+5*sum(-5+10*x)+39.16616*2*n); % 0<x<1
funplotting2d(fun,xx,xx,n,0)
figure(2);
hold on 
plot([flin2(a),a],[a,flin2(a)],'k-','linewidth',2)
plot([0,flin1(0)],[flin1(0),0],'k-','linewidth',2)
rectangle('Position', [0 0 a a])
plot(-b0,-b0,'r*')
% figure_to_publish('./figures/example3_active_const')
% figure_to_publish('./figures/example3_interior')

