function sse = fit8(c)
%c=[10;c(1:5);100;c(6)];
c=c';
%c=[c(1:2) 60 c(3:7)];
%   C index 
%   1-c1_mc,2-c2_mc
%   3-c1_collagen1,4-c2_collagen1
%   5-c1_collagen2,6-c2_collagen2
%   7-c1_collagen3,8-c2_collagen3
%   9-c1_collagen4,10-c2_collagen4
%  [AV] = 7 - elastin 8 - angle
% load MT2_2.txt;
% 
%  h=(MT2_2(:,6)-MT2_2(:,7))*1e-3;
%  lz=MT2_2(:,4);
%  lt=MT2_2(:,3);
%  s_tt=MT2_2(:,1)*1e3;
%  s_zz=MT2_2(:,2)*1e3;
 %load c.txt;
 load MT3_2.txt;

 h=(MT3_2(:,6)-MT3_2(:,7))*1e-3;
 lz=MT3_2(:,4);
 lt=MT3_2(:,3);
 s_tt=MT3_2(:,1)*1e3;
 s_zz=MT3_2(:,2)*1e3;
 
 %s_tt_avg=mean(s_tt);
 %s_zz_avg=mean(s_zz);
 
%Fibre Angle
alpha_mp=pi/2;
alpha_co1=0;
alpha_co2=pi/2;
alpha_co3=c(8);
alpha_co4=-c(8);


% Prestretch
% G_e=1.4;
% G_c=1.08;
% G_mp=1.2;

G_e_cir=1.4;
G_e_axi=1.6;
G_e=1.4;
G_c=1.08;
G_mp=1.2;


%mass fractions
 m_c=.45; m_e=0.12;m_m=.43;
 m_c1=.1*m_c;m_c2=.1*m_c;m_c3=.4*m_c;m_c4=.4*m_c;


% For elastin work function
%c_e=382.4;%kPA

%For active SMC
Tmax=150;
lM=1.1;l0=0.4; CB=0.68;

% l - lambda, lz - lambda_z, lt - lambda_theta
l_mp=sqrt((lz*cos(alpha_mp)).^2+(lt*sin(alpha_mp)).^2);
l_co1=sqrt((lz*cos(alpha_co1)).^2+(lt*sin(alpha_co1)).^2);
l_co2=sqrt((lz*cos(alpha_co2)).^2+(lt*sin(alpha_co2)).^2);
l_co3=sqrt((lz*cos(alpha_co3)).^2+(lt*sin(alpha_co3)).^2);
l_co4=sqrt((lz*cos(alpha_co4)).^2+(lt*sin(alpha_co4)).^2);
%  -----------------------------------------Thetha Stress------------------

%Passive SMC
W_mp=c(1)*(exp(c(2)*(G_mp^2*l_mp.^2-1).^2));
dexpdlt_mp=c(2)*(G_mp^4*(4*(lt.^3).*(sin(alpha_mp)).^4+4*lt.*(lz.^2).*(sin(alpha_mp)).^2*(cos(alpha_mp)).^2)-G_mp^2*(4*lt.*(sin(alpha_mp)).^2));

dWdlt_mp=W_mp.*dexpdlt_mp;%

%Collagen 1 at 0 degrees
W_co1=c(3)*(exp(c(4)*(G_c^2*l_co1.^2-1).^2));
dexpdlt_co1=c(4)*(G_c^4*(4*(lt.^3).*(sin(alpha_co1)).^4+4*lt.*lz.^2*(sin(alpha_co1)).^2*(cos(alpha_co1)).^2)-G_c^2*(4*lt.*(sin(alpha_co1)).^2));

dWdlt_co1=W_co1.*dexpdlt_co1;

%Collagen 2 
W_co2=c(1)*(exp(c(2)*(G_c^2*l_co2.^2-1).^2));
dexpdlt_co2=c(2)*(G_c^4*(4*lt.^3*(sin(alpha_co2)).^4+4*lt.*lz.^2*(sin(alpha_co2)).^2*(cos(alpha_co2)).^2)-G_c^2*(4*lt.*(sin(alpha_co2)).^2));

dWdlt_co2=W_co2.*dexpdlt_co2;

%Collagen 3
W_co3=c(5)*(exp(c(6)*(G_c^2*l_co3.^2-1).^2));
dexpdlt_co3=c(6)*(G_c^4*(4*(lt.^3).*(sin(alpha_co3)).^4+4*lt.*lz.^2*(sin(alpha_co3)).^2*(cos(alpha_co3)).^2)-G_c^2*(4*lt.*(sin(alpha_co3)).^2));

dWdlt_co3=W_co3.*dexpdlt_co3;

%Collagen 4
W_co4=c(5)*(exp(c(6)*(G_c^2*l_co4.^2-1).^2));
dexpdlt_co4=c(6)*(G_c^4*(4*lt.^3*(sin(alpha_co4)).^4+4*lt.*lz.^2*(sin(alpha_co4)).^2*(cos(alpha_co4)).^2)-G_c.^2*(4*lt*(sin(alpha_co4)).^2));

dWdlt_co4=W_co4.*dexpdlt_co4;

% Elastin

dWdlt_e=c(7)*(2.*lt.*G_e_cir^2-2./(G_e_cir^4*lt.^3.*lz.^2));

% Active SMC
% l_act=1; 
% sact_theta=Tmax*m_m*(1-exp(-CB^2))*l_act*(1-(lM-l_act)^2/(lM-l0)^2);

%s_theta=lt.*(m_e*dWdlt_e+m_c1*dWdlt_co1+m_c2*dWdlt_co2+m_c3*dWdlt_co3+m_c4*dWdlt_co4+m_m*dWdlt_mp);
s_theta=lt.*(dWdlt_e+dWdlt_co1+dWdlt_co2+dWdlt_co3+dWdlt_co4+dWdlt_mp);


%------------------------------------ Z Stress-----------------------------

%Passive SMC
W_mp=c(1)*(exp(c(2)*(G_mp^2*l_mp.^2-1).^2));
dexpdlz_mp=c(2)*(G_mp^4*(4*(lz.^3).*(cos(alpha_mp)).^4+4*lz.*(lt.^2).*(sin(alpha_mp)).^2*(cos(alpha_mp)).^2)-G_mp^2*(4*lz.*(cos(alpha_mp)).^2));

dWdlz_mp=W_mp.*dexpdlz_mp;%

%Collagen 1 at 0 degrees
W_co1=c(3)*(exp(c(4)*(G_c^2*l_co1.^2-1).^2));
dexpdlz_co1=c(4)*(G_c^4*(4*(lz.^3).*(cos(alpha_co1)).^4+4*lz.*lt.^2*(sin(alpha_co1)).^2*(cos(alpha_co1)).^2)-G_c^2*(4*lz.*(cos(alpha_co1)).^2));

dWdlz_co1=W_co1.*dexpdlz_co1;

%Collagen 2 
W_co2=c(1)*(exp(c(2)*(G_c^2*l_co2.^2-1).^2));
dexpdlz_co2=c(2)*(G_c^4*(4*(lz.^3)*(cos(alpha_co2)).^4+4*lz.*lt.^2*(sin(alpha_co2)).^2*(cos(alpha_co2)).^2)-G_c^2*(4*lz.*(cos(alpha_co2)).^2));

dWdlz_co2=W_co2.*dexpdlz_co2;

%Collagen 3
W_co3=c(5)*(exp(c(6)*(G_c^2*l_co3.^2-1).^2));
dexpdlz_co3=c(6)*(G_c^4*(4*(lz.^3).*(cos(alpha_co3)).^4+4*lz.*lt.^2*(sin(alpha_co3)).^2*(cos(alpha_co3)).^2)-G_c^2*(4*lz.*(cos(alpha_co3)).^2));

dWdlz_co3=W_co3.*dexpdlz_co3;

%Collagen 4
W_co4=c(5)*(exp(c(6)*(G_c^2*l_co4.^2-1).^2));
dexpdlz_co4=c(6)*(G_c^4*(4*(lz.^3).*(cos(alpha_co4)).^4+4*lz.*lt.^2*(sin(alpha_co4)).^2*(cos(alpha_co4)).^2)-G_c^2*(4*lz.*(cos(alpha_co4)).^2));

dWdlz_co4=W_co4.*dexpdlz_co4;

% Elastin

dWdlz_e=c(7)*(2.*lz.*G_e_axi^2-2./(G_e_axi^4*lt.^2.*lz.^3));

%s_z=1./lt.*(m_e*dWdlz_e+m_c1*dWdlz_co1+m_c2*dWdlz_co2+m_c3*dWdlz_co3+m_c4*dWdlz_co4+m_m*dWdlz_mp);
s_z=lz.*(dWdlz_e+dWdlz_co1+dWdlz_co2+dWdlz_co3+dWdlz_co4+dWdlz_mp);

%== [AV]

%figure(1);plot(lt,s_theta,lt,s_tt,'r','LineWidth',2);
%xlabel('lambda-theta');ylabel('sigma-theta(Pa)')
%legend('fit','experiment')
%figure(2);plot(lt,s_z,lt,s_zz,'r','LineWidth',2);
%xlabel('lambda_z');ylabel('sigma_z(Pa)')
% ==[/AV]

% legend('fit','experiment')
% figure;plot(lt,dWdlz_e,'*',lt,dWdlz_mp,'^',lt,dWdlz_co1,'v',lt,dWdlz_co2,'o',lt,dWdlz_co3,'-',lt,dWdlz_co4);
% Error_Vector=(((s_theta-s_tt).^2)+((s_zz-s_z).^2));
%%sse=1/length(s_theta)*sqrt(sum(Error_Vector));
Error_Vector=(s_z-s_zz).^2./s_zz.^2+(s_tt-s_theta).^2./s_tt.^2;

%Error_Vector=((s_theta-s_tt).^2)./(s_tt).^2;
%Error_Vector=(((s_theta-s_tt).^2)/s_tt_avg^2+((s_zz-s_z).^2)/s_zz_avg^2);
%sse=sqrt(sum(Error_Vector));
sse=sqrt(sum(Error_Vector)/length(s_theta));
sse=min(sse,1000);

dlmwrite('sse.txt',sse);

end
