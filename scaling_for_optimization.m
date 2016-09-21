%% finding the scaling factor for Z and Lambda 2<N<8

clear all
close all
clc
global neigh lattice N options
N=2; 
RL=[];
for N=2:1:8
    
    switch(N)
        
        case(2)
            lattice1='An ';
            Dz=0.78540;
            Da=0.90690;
        case(3)
            lattice1='An ';
            Dz=0.52360;
            Da=0.74048;
        case(4)
            lattice1='Dn ';
            Dz=0.30843;
            Da=0.61685;
        case(5)
            lattice1='Dn ';
            Dz=0.16449;
            Da=0.46526;
        case(6)
            lattice1='E6 ';
            Dz=0.08075;
            Da=0.37295;
        case(7)
            lattice1='E7 ';
            Dz=0.03691;
            Da=0.29530;
        case(8)
            lattice1='E8 ';
            Dz=0.01585;
            Da=0.25367;
    end

    
    da=(Da/Dz)^(1/N); 
    lambda{N}.Da = Da;
    
    lambda{N}.Dz = Dz;
    RL=[RL , da];
    end
    Scale = [RL; [2:8]];

 
 

Sn=1;
dz=1;
%%%%%%%%%%%%%%%%%%%%%%%%%
da=(Da/Dz)^(1/N)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% de=(De/Dz)^(1/N)
% 
%     oscale=1;
%     
%     lattice='An ';
%     
%     [neigh,matrix,plane]=init_DOGS(N,lattice);