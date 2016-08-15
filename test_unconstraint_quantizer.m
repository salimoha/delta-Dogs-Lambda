    % test the unconstraint quantizer
    % shahrouz alimo August 11 2016
    clear all;clc;
    %%
    %lattice= 'Zn ' , 'An ', 'An*?!','Dn ','Dn*', 'E8 '
    global lattice plane 
    % initialization 
    ii=0
    while ii<20
    n=2;  S=rand(n,1) % data point
    %%
    lattice='An ';   [~,B,plane]=init_DOGS(n,lattice);
    %scale factor. Different for each lattice!
    scale=1;
    %
    %quantization 
    [Sq,Errq,Sz]=Unconstraint_quantizer(S,scale);
    
    X = [S, Sq]
    e = Errq
    ii=ii+1;
    saveas(gcf, strcat('./pics2D/example_',num2str(ii,'%02d')))
    end