% Second Order Delta-Dogs:
% Algorithm with second order convergence property.
% Constraints are simple bound constraints.
% Dimenson:
clear all; close all; clc
load Lambda
global n Ain bin y0 lattice
% for kk=[2,3,4,5,8]
for kk=[4]
clear yE xE Nm 
n=kk;
figure(1);clf;
qtype=0; iPlot=1; iSave=1;
run('./FilesToPath/addTopath.m')
% minCoveringRadius = 5e-6;
minCoveringRadius = 5e-3;
% EXAMPLE = 'interior'; bs = 0.81; % interior
EXAMPLE = 'activ'; bs=0.7903; % active on constraint
%  n=3;
%  lattice = 'Zn '; % cartesian Z
 lattice = lambda{kk}.lattice
%%%%%%%%%best quantizers and densest %%%%%%%%%%%%%%%
% n=1; lattice = 'Zn ';
% n=2;lattice = 'An ';
%  n=3; lattice = 'An*'; %%%BETTER %lattice = 'An*';
% n=3; lattice = 'Dn ';
% n=5; lattice = 'An*'; %lattice = 'Dn ';
% n=5; lattice = 'Dn*';
% n=6; lattice = 'An '; % lattice = 'E6 ';

% n=6; lattice = 'E6*'; % lattice = 'E6 ';
% n=7; lattice = 'E7*'; % lattice = 'E7 ';
% n=8; lattice = 'E8';
meshNAme=char(lattice);
if meshNAme(end)=='*', Dual =1; else, Dual =0; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% xs = ones(n,1)*(-2.903);
% objective function:
% scaled to [0,1] domain
% xs =  (5+xs)/10;
xs = ones(n,1)*0.2097;
%0<x<1
fun=@(x)(sum((-5+10*x).^4)-16*sum((-5+10*x).^2)+5*sum(-5+10*x)+39.16616*2*n);
%-5<x<5
% xs = ones(n,1)*(-2.903);
% fun=@(x) (sum((x).^4)-16*sum((x).^2)+5*sum(x)+39.16616*2*n);
% Constraints:
% upper nad lower bounds
% lob=-5*ones(n,1); upb=5*ones(n,1);
lob=zeros(n,1); upb=ones(n,1);
% Estimate of the solution:
y00=0; %y00 = -1e-5 ; % y00 = 0.01;
% Initial Calculation:
bnd1 = lob; bnd2 = upb;
ub = bnd2; lb = bnd1;
m=2*n;
% maximum number of iterations:
iter_max=100;
% interpolation strategy:
inter_method=1;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Input the equality constraints:?????????
% % Ain=[eye(n);-eye(n)]; bin=[bnd2 ;-bnd1];
% % Ain=[eye(n); -ones(1,n)/sqrt(n)]; bin=[5*ones(n,1); 2.8*n/sqrt(n)];
% % A=[ -ones(1,n)/sqrt(n)]; b=[ 2.8*n/sqrt(n)];
% A=([ -ones(1,n)/sqrt(n)]); b=[ 2.*n/sqrt(n)];
% Ain =[A; eye(n); -eye(n)];
% bin=[b ; ub; -lb];
%%%%%%%%%%%%%%%%%%%%%%%%%%%
A=[-ones(1,n)]/sqrt(n); b=[-xs(1).*n/sqrt(n)];
Ain=[eye(n);-eye(n);A]; bin=[ones(n,1);zeros(n,1);b];
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% length_scale=10;
% length_scale=1;
length_scale=max(bnd2-bnd1);
% Calculate the Initial trinagulation points
% xU=bounds(bnd1,bnd2, n);
% xU=vertex_find(Ain,bin,[],[]);
xU=vertex_find(Ain,bin,lb,ub);
% Calculate initial delta
xE=(bnd1+0.45*(bnd2-bnd1));
% xE=(sum(xU')/(n+1))';
% initial points: the midpoint point and its neighber
delta0=0.15*length_scale;
for ii=1:n; e=zeros(n,1); e(ii)=1;xE(:,ii+1)=xE(:,1)+delta0*e; end
% Calculate initial unevaluated delta
for ii=1:size(xU,2),  deltaU(ii)=mindis(xU(:,ii),xE); end


% Calculate the function at initial points
for ii=1:size(xE,2),   yE(ii)=fun(xE(:,ii)); end
pE=yE;
interpolate_index=ones(1,length(yE));
%Nm=800;
scale = lambda{kk}.scale;
Nm=4; Nm= Nm/scale;
% Nm=100;
rho=1;
for kk=1:16
    for k=1:iter_max
        %     keyboard
        stop=0;
        inter_par= interpolateparametarization(xE,yE,inter_method,interpolate_index);
        %y0=10;
        y0=y00;
        %y0=min(yE)-range(yE)/Nm;
        ys_disc=[];
        yp_U=[];
        % Calculate the metric of the unevaluated points.
        for ii=1:size(xU,2)
            ys_disc(ii)=(interpolate_val(xU(:,ii),inter_par)-y0)/mindis(xU(:,ii),xE);
            yp_U(ii)=interpolate_val(xU(:,ii),inter_par);
        end
        %  keyboard
        % check the minimum value of xU
        if size(xU,2)>0
            [tup,ind_up]=min(ys_disc);
        else
            tup=inf;
        end
        if tup<0
            % if there is a point that is ys_disc < 0 then
            [trp,ind_rp]=min(yp_U);
            
            % trust region decreasing step
            [tt,ind]=min(yE); x0=xE(:,ind); %x0=xU(:,ind_rp);
            x=min_decrease(x0, xU(:,ind_rp), inter_par);
            
            %keyboard
            x=General_linear_quatization(x, Ain, bin, Nm/length_scale, qtype);
            %keyboard
            
            %      when a new point is close to xU prject x to xU. ????
            if mindis(x,xU)<1e-4, x=xU(:,ind_rp); xU(:,ind_rp)=[]; newadd=0; end
            
        else
            % perform the trinagulation search
            while 1
                %    keyboard
                %  minimizer of the continuous search function. x_k step 4
                x=tringulation_search_bound(inter_par,[xE xU],yE);
                [xq]=General_linear_quatization(x, Ain, bin, Nm/length_scale, qtype);
                
                %it the interpolation was less then y0
                if interpolate_val(x,inter_par)<y0, break, end
                
                %      keyboard
%                 [xq,xE,xU,newadd,success]=points_neighbers_find(xq,xE,xU);
                
                [xq,xE,xU,newadd,success]=points_neighbers_find_lambda(xq,xE,xU);
                
                
                if success==1, break,
                    %
                else
                    ys_disc=[ys_disc, (interpolate_val(xq,inter_par)-y0)/mindis(xq,xE)];
                    yp_U=[yp_U, interpolate_val(xq,inter_par)];
                end
                
                
            end
            
            %      minimizer s_cont    >   minimizer s_disc
            if (interpolate_val(x,inter_par)-y0)/mindis(x,xE) > min(ys_disc)
                % consider the minimizer of discrete seach function
                x = xU(:,ind_up); xU(:,ind_up)=[]; newadd=0;
            else
                % consider the minimizer of quantized continious seach function
                x = xq;
            end
            
            
            
        end
        
        
        %   true=check_add_point(x,ym,[xE xU],inter_par,newadd,size(xE,2));
        if mindis(x,xE) < minCoveringRadius/scale
            %  keyboard
            break
        end
        
        
        %  perform function evalution
        ym=fun(x); xE=[xE x]; yE=[yE ym];
        
        % interpolate_index=[interpolate_index true];
        pE=[pE interpolate_val(x,inter_par)];
        
        %       [t,ind]=min(yE); x0=xE(:,ind);delta_mesh=1/Nm;
        % if (interpolate_val(x,inter_par)-y0)/(mindis(x,xE))>(interpolate_val(x0,inter_par)-y0)/(delta_mesh);
        %     break, end
        
        if iPlot ==1
            % illusteration
            figure(gcf)
            %subplot(2,1,1)
            plot(1:length(yE),yE,'-','linewidth',2.5)
            xlabel('number of fun eval'); ylabel('fun eval')
            grid off;
            set(gca,'FontSize',18)
            %, 1:length(yE),pE,'--')
            %plot(1:length(yE),1./(exp(yE)-1),'-','linewidth',3)
            ylim([0 200])
            %ylim([-5 5])
            %xlim([0 100])
            %subplot(2,1,2)
            %plot(xE')
            %xlim([0 100])
            drawnow
            %plot(xi(1,:)0,xi(2,:),'rs')
        end
    end
    Nm=2*Nm/scale;
end
%
% figure(2);clf
% plot(1:length(yE),yE(1:end),'-','linewidth',2)
% axis([0 length(yE)*1.2 -1 150])
% grid on
% set(gca,'fontsize',24)
% axes('Position',[0.7,0.3,0.28,0.28])
% t=floor(.8*length(yE))
% plot(t:length(yE),yE(t:end),'-','linewidth',1.5)
% set(gca,'fontsize',14)
% grid on
figure(2);clf; plotting_optimal_point(yE,'semilog')
%
POINT='BODCENTR';

fpath = '/Users/shahrouz/Desktop';
% POINT='GooD';
if iPlot ==1
    switch EXAMPLE
        case('interior')
            if Dual==1
                CC=char(lattice); CC=CC(1:end-1);
                figure_to_publish(strcat('./results_Lambda/example_interior_bs081/ex1_', char(POINT), 'interior_Lambda_',  char(CC) ,'_dual', '_n_', num2str(n,'%02d')))
                
            else
                figure_to_publish(strcat('./results_Lambda/example_interior_bs081/ex1_', char(POINT), 'interior_Lambda_', char(lattice), '_n_', num2str(n,'%02d')))
            end
            
        case('activ')
            if Dual==1
                CC=char(lattice); CC=CC(1:end-1);
                figure_to_publish(strcat('./results_Lambda/example_3_activ_bs07903/ex3_activ_Lambda_', char(CC) ,'_dual', '_n_', num2str(n,'%02d')))
                %         save(strcat('./results_Lambda/example_3_activ_bs07903/ex3_activ_Lambda_', char(CC) ,'_dual', '_n_', num2str(n,'%02d')));
                
            else
                figure_to_publish(strcat(fpath,'/results_Lambda/example_3_activ_bs07903/ex3_activ_Lambda_',char(lattice),'_n_', num2str(n,'%02d')))
            end
    end
    
end
% saving data
if iSave ==1
    switch EXAMPLE
        case('interior')
            if Dual==1
                CC=char(lattice); CC=CC(1:end-1);
                save(strcat('./results_Lambda/example_interior_bs081/ex1_', char(POINT), 'interior_Lambda_', char(CC) ,'_dual', '_n_', num2str(n,'%02d')))
                
            else
                save(strcat('./results_Lambda/example_interior_bs081/ex1_', char(POINT), 'interior_Lambda_', char(lattice), '_n_', num2str(n,'%02d')))
            end
            
        case('activ')
            if Dual==1
                CC=char(lattice); CC=CC(1:end-1);
                save(strcat('./results_Lambda/example_3_activ_bs07903/ex3_activ_Lambda_', char(CC) ,'_dual', '_n_', num2str(n,'%02d')))
                %         save(strcat('./results_Lambda/example_3_activ_bs07903/ex3_activ_Lambda_', char(CC) ,'_dual', '_n_', num2str(n,'%02d')));
                
            else
                save(strcat(fpath, '/results_Lambda/example_3_activ_bs07903/ex3_activ_Lambda_',char(lattice),'_n_', num2str(n,'%02d')))
            end
    end
    
end

end