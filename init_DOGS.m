%calculate the cloud of nearest neighbors, the basis matrix, and the plane
%that the lattice lies on, if applicable. Optioins for the variable lattice
%are 'An', 'An*', 'Dn', 'Dn*', 'E8 '

function [neigh,matrix,plane]=init_DOGS(N,lattice);

disp(strcat('initializing  ',' lattice ' ,'....  ',  '   ', char(lattice),'  ....'))

[matrix,v,len] = make_matrix(lattice,N);   %build the lattice basis matrix
plane          = [];
iter           = 2;

%if(lattice=='An '), iter=1; end

%find neighbors in matrix space
indx           = -iter*ones(1,length(matrix(1,:))); 
indxorg        = indx; 
c              = 0; 
b              = 0;



neigh=[];
indx(end)=indx(end)-1;


for j=1:(2*iter+1)^(N), 
    
    indx(end)=indx(end)+1; 

    for i=length(indx):-1:2
        if(indx(i)>-indxorg(i))
            indx(i)=indxorg(i);
            indx(i-1)=indx(i-1)+1;
        end
        
    end
    


    pt=matrix*indx';
    v=norm(pt);
    if(v>0 && v <1.001*len)
        neigh=[neigh; pt'];
    end
end


%done with neighbors in matrix space
plane=[];
switch(lattice)
    case('An ')
        plane = QRHouseholder(ones(N+1,1)); plane=plane(2:end,:);  
        %QRHouseholder returns 2:end orthogonal vectors - plane basis
        for i=1:length(neigh)
            neigh2(i,:)=(plane'\neigh(i,:)')';
        end
        neigh=neigh2;
    case('An*')
        v=ones(N+1,1); v(end)=-N;
        plane = QRHouseholder(v); plane=plane(2:end,:);
        for i=1:length(neigh)
            neigh2(i,:)=(plane'\neigh(i,:)')';
        end
        neigh=neigh2;
end

disp('initialization complete, DELTA DOGS Lambda starting...')

function [Q,A] = QRHouseholder(A)
% Compute a QR decomposition A=QR by applying a sequence of Householder reflections
% to any MxN matrix A to reduce it to upper triangular form.

[M,N]=size(A); Q=eye(M,M);
for i=1:min(N,M-1)
    [A(i:M,i:N),sigma,w] = Reflect(A(i:M,i:N));
    wdot=w';
    Q(:,i:M)=Q(:,i:M)-(Q(:,i:M)*w)*(sigma*wdot);
end % Eqn (1.9b)
%end
% end function QRHouseholder.m


function [X,sig,w] = Reflect(X)
% Apply a Householder reflector matrix H?H to a MxN matrix X (i.e., calculate H?H*X),
% with [sigma,w] arranged to give zeros in the (2:end,1) locations of the result.
x=X(:,1); if real(x(1))<0, s=-1; else, s=1; end, nu=s*norm(x); % Eqn (1.7b)
if nu==0, sig=0; w=0; else, sig=(x(1)+nu)/nu; w=[x(1)+nu; x(2:end)]/(x(1)+nu); end
X(1,1)=-nu; X(2:end,1)=0; % Eqn (1.8)
wdot=w';
X(:,2:end)=X(:,2:end)-(conj(sig)*w)*(wdot*X(:,2:end)); % Eqn (1.9a)
% end function Reflect.m


%%%%%% Make_Matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [matrix,v,len] = make_matrix(lattice,N)
plane=[];
switch(lattice)
    case('Zn ')
        matrix = eye(N);
        len=1;
        v=[];
  
    case('An ')
        v=ones(N,1);
        M=diag(v,0);
        M(:,length(M)+1)=zeros(1,N);
        P=diag(v,1);
        P=P(1:N,:);
        matrix = (-M+P)';
        len=1;%sqrt(2);
        matrix=matrix./sqrt(2); 
         
    case('Dn ')
        v=-1.*ones(N,1); 
        matrix=diag(v,0);
        P=diag(-v,1);  P(end,:)=[]; P(:,end)=[];
        matrix=matrix+P;
        matrix(2,1)=-1;
        v=[];
        len=sqrt(2);
        
    case('Dn*')
        v=ones(N,1); 
        matrix=diag(v,0);
        matrix(:,end)=0.5.*ones(1,N);
        v=[];
        len1=norm(matrix(:,1));
        len2=norm(matrix(:,end));
        len=min(len1,len2);
    case('An*')
        P=diag(-ones(N,1),-1);
        P(:,end)=[];
        P(1,1:end)=1; 
        P(:,end)=(1/(N+1))*ones(N+1,1);P(1,end)=-N/(N+1);
        matrix=P;
        v=ones(N+1,1);
        len=norm(matrix(:,end));
    case('E8 ')
        matrix = diag(ones(1,8));
        matrix = matrix + diag(-1.*ones(7,1),1);
        matrix(1:4,end)=.5;
        matrix(4:8,end)=-.5;
        matrix(1,1)=2;
        v=[];
        len=sqrt(2);
        
    otherwise
        err('no viable lattice entered')
end
