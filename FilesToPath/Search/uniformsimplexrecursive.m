function [V,a] = uniformsimplexrecursive(N)                                                                   
%unifrom simplex around origin with radius 1.
%keyboard
if N==1
    V=[-1/2 1/2];
    a=0.5;
end
if N>1
    [V,a]=uniformsimplexrecursive(N-1);
    x=(1-2*a^2)/(2*sqrt(1-a^2));
    a=sqrt(a^2+x^2);
    V=[V;-x*ones(1,N)];
    V=[V,[zeros(N-1,1);a]];
end


