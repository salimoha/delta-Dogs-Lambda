%% Restriction Algorithm Library
%  \Delta-DOGS (\Lambda) Library
%
%  Title: Quantization process generation for unconstraint lattices
%
%  Author: Shahrouz Alimo
%
%  Description:
%    Domain setup.
%    This function will restrict a point x onto the lattice defined by the variable 'lattice'
%defined by the coarseness variable 'scale = 1/MeshSize'.
%
%
%
function [Sq,Errq,Sz, S]=Unconstraint_quantizer(S,scale) %S Y
% Input:
% S: Intial set of points, note that: S = [x1;x2;...;xn]. Each column of
% the point matrix S is a one point. (like mathematical vectors).
%
% scale: the coarseness variable 'scale = 1/MeshSize'.
%
% lattice: type of desired lattice:
%   lattice=1: Zn or Cartesian lattice
%              Z2: square
%              Z3: cube
%   lattice=2: An or Zero-sum lattice. see Conway&Sloane
%              A2: hexagonal
%              A3: face-centered cubic (FCC)
%              A8: zero-sum
%   lattice=3: An* or Dual of zero-sum lattice.
%
%   lattice=4: Dn or Checkerboard lattice. %page 445 conway and sloane
%              D2: rotated squares
%   lattice=5: Dn* or Dual of Checkerboard lattice.
%
%   lattice=6: En or Gosset lattice.
%
%   lattice=7: En* or Dual of Gosset lattice.
%

%% check the inputs
global plane  lattice 

%% quantization process
% keyboard
for k=1:length(S(1,:))                     %iterate through matrix of pts
    x=S(:,k);
    Sz(:,k)=[x;nan];
    N=length(x);
    %pt of interest
    switch(lattice)
        case('Zn ')                        %cartesian lattice
            x = x.*scale;                  %scale x
            x = round(x);                  %round to integer
            x=x./scale;                    %scale back to grid
            %% Zero sum,
            %         case(2)
            %
            %             % Zero sum, Shahrouz
            %             n=length(x);
            % %             x=x./scale.*sqrt(2);                    %scale to integers
            %                     x=x./scale;                    %scale to integers
            % %                 keyboard
            %             V=null(ones(1,n+1));
            %             x=V*x;
            %             xq=round(x);
            %             delta=sum(xq);
            %             if delta~=0
            %                 d=x-xq; [~,ind]=sort(d);
            %                 if delta>0
            %                     xq(ind(1:-delta))=xq(ind(1:-delta))-1;
            %                 else
            % %                     xq(ind(n-delta+2:n+1))=xq(ind(n-delta+2:n+1))+1;
            %                        xq(ind(n-delta+1:n+1))=xq(ind(n-delta+1:n+1))+1;
            %                 end
            %             end
            % %             keyboard
            %             xq=V\xq;
            % %             x=(xq.*scale./sqrt(2));
            %                x=(xq.*scale);
            %             % Paul
            % %         case(3)
        case('An ')                        %An lattice, see Conway&Sloane
            x=x.*scale.*sqrt(2);                    %scale to integers
            x=(plane'*x)';
            %             keyboard
            %algorithm from Conway & Sloane
            s=sum(x); x_p = x - s/(N+1)*ones(1,N+1);  %project onto plane
            f_x = round(x_p); DELTA = sum(f_x);  %calculate deficiency
            delta = x_p - f_x;  [delta, i] = sort(delta);
            if(DELTA)>0, for j=1:DELTA, f_x(i(j)) = f_x(i(j))-1;end, end
            if(DELTA)<0, for j=N-DELTA:N+1, f_x(i(j)) = f_x(i(j))+1;end,end
            %back into N
            x=(plane'\f_x'./scale./sqrt(2));
            %this is all straight out of C&S except for dividing by sqrt(2)
            %in the last step. This is necessary because otherwise the
            %algorithm will restrict to the lattice An which has a base
            %vector length of sqrt(2).
            Sz(:,k)=f_x./scale./sqrt(2);
            %% Checkerboard
        case('Dn ') %Dn
            %             % Checkerboard, Shahrouz
            %             xq=round(x./scale);
            %             %         xq=round(x);
            % %             xq=round(x);
            %             if mod(sum(xq),2)==1
            %                 dq=x-xq; [~,ind]=max(abs(dq));
            %                 xq(ind)= xq(ind)+sign(dq(ind)+eps);
            %             end
            %             x=xq.*scale;
            %         case(5) % Dn
            %page 445 conway and sloane
            x=x.*scale;
            f_x = round(x);
            g_x=f_x; d=x-g_x; [a,j]=max(abs(d));
            if(x(j)-g_x(j)>0), g_x(j)=g_x(j)+1;
            elseif(x(j)-g_x(j)<=0), g_x(j)=g_x(j)-1; end
            a=sum(f_x); b=sum(g_x);
            x=f_x;
            if(rem(a,2)==0),x=f_x; else, x=g_x; end,
            x=x./scale;
            %% Dual of An
        case('An*');
            % case(6)
%             TODO
            %% Dual of Dn*
            %             keyboard
        case('Dn*');
            x=x';
            x=x.*scale; x1=round(x); x2=(round(x-0.5.*ones(1,N))+0.5.*ones(1,N));
            if(norm(x-x2)<norm(x-x1)), x=x2; else x=x1; end, x=x./scale;
            x=x';
            %% Gosset lattice
            %             case(7)
        case('E6 ');
            
        case('E7 ');
            
            
        case('E8 ');
%             keyboard
            x=x';
            %first quantize to Dn*
            x=x.*scale; xorg=x;
            xb=x;
            x1=round(x);
            x2=(round(x-0.5.*ones(1,N))+0.5.*ones(1,N));
            if(norm(x-x2)<norm(x-x1)), x=x2; else x=x1; end
            xa=x;  %save x as xa
%   keyboard           
            %now quantize to Dn+[1]
            x=xb+.5.*ones(1,N);
            x1=round(x);
            x2=(round(x-0.5.*ones(1,N))+0.5.*ones(1,N));
            if(norm(x-x2)<norm(x-x1)), x=x2; else x=x1; end
            xb=x;  %save x as xb
            
            %figure out which the original x is closer to
            d1=norm(xorg-xa);
            d2=norm(xorg-xb);
            
            if(d1<d2)
                x=xa;
            else
                x=xb;
            end
            x=x./scale;
            x=x';
            
    end
    
    Sq(:,k)=x;
    Errq(k) = norm(x-S(:,k));
    
end
